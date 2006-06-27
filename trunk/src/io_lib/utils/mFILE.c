#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdarg.h>

#include "os.h"
#include "mFILE.h"
#include "vlen.h"

/*
 * This file contains memory-based versions of the most commonly used
 * (by io_lib) stdio functions.
 *
 * Actual file IO takes place either on opening or closing an mFILE.
 *
 * Coupled to this are a bunch of rather scary macros which can be obtained
 * by including stdio_hack.h. It is recommended though that you use mFILE.h
 * instead and replace fopen with mfopen (etc). This is more or less
 * mandatory if you wish to use both FILE and mFILE structs in a single file.
 */

static mFILE *m_channel[3];  /* stdin, stdout and stderr fakes */

/*
 * Reads the entirety of fp into memory. If 'fn' exists it is the filename
 * associated with fp. This will be used for more optimal reading (via a
 * stat to identify the size and a single read). Otherwise we use successive
 * reads until EOF.
 *
 * Returns a malloced buffer on success of length *size
 *         NULL on failure
 */
static char *mfload(FILE *fp, const char *fn, size_t *size, int binary) {
    struct stat sb;
    char *data = NULL;
    size_t allocated = 0, used = 0;
    int bufsize = 8192;

#ifdef _WIN32
    if (binary)
	_setmode(_fileno(fp), _O_BINARY);
    else 
	_setmode(_fileno(fp), _O_TEXT);
#endif

    if (fn && -1 != stat(fn, &sb)) {
	data = malloc(allocated = sb.st_size);
	bufsize = sb.st_size;
    } else {
	fn = NULL;
    }

    do {
	size_t len;
	if (used + bufsize > allocated) {
	    allocated += bufsize;
	    data = realloc(data, allocated);
	}
	len = fread(data + used, 1, allocated - used, fp);
	if (len > 0)
	    used += len;
    } while (!feof(fp) && (fn == NULL || used < sb.st_size));

    *size = used;

    return data;
}

/*
 * Creates and returns m_channel[0].
 * We initialise this on the first attempted read, which then slurps in
 * all of stdin until EOF is met.
 */
mFILE *mstdin(void) {
    if (m_channel[0])
	return m_channel[0];

    m_channel[0] = mfcreate(NULL, 0);
    m_channel[0]->fp = stdin;
    return m_channel[0];
}

static void init_mstdin(void) {
    static int done_stdin = 0;
    if (done_stdin)
	return;

    m_channel[0]->data = mfload(stdin, NULL, &m_channel[0]->size, 1);
    m_channel[0]->mode = MF_READ;
    done_stdin = 1;
}

/*
 * Creates and returns m_channel[1]. This is the fake for stdout. It starts as
 * an empty buffer which is physically written out only when mfflush or
 * mfclose are called.
 */
mFILE *mstdout(void) {
    if (m_channel[1])
	return m_channel[1];

    m_channel[1] = mfcreate(NULL, 0);
    m_channel[1]->fp = stdout;
    m_channel[1]->mode = MF_WRITE;
    return m_channel[1];
}

/*
 * Stderr as an mFILE.
 * The code handles stderr by returning m_channel[2], but also checking
 * for stderr in fprintf (the common usage of it) to auto-flush.
 */
mFILE *mstderr(void) {
    if (m_channel[2])
	return m_channel[2];

    m_channel[2] = mfcreate(NULL, 0);
    m_channel[2]->fp = stderr;
    m_channel[2]->mode = MF_WRITE;
    return m_channel[2];
}


/*
 * For creating existing mFILE pointers directly from memory buffers.
 */
mFILE *mfcreate(char *data, int size) {
    mFILE *mf = (mFILE *)malloc(sizeof(*mf));
    mf->fp = NULL;
    mf->data = data;
    mf->alloced = size;
    mf->size = size;
    mf->eof = 0;
    mf->offset = 0;
    mf->flush_pos = 0;
    mf->mode = MF_READ | MF_WRITE;
    return mf;
}

/*
 * Recreate an existing mFILE to house new data/size.
 * It also rewinds the file.
 */
void mfrecreate(mFILE *mf, char *data, int size) {
    if (mf->data)
	free(mf->data);
    mf->data = data;
    mf->size = size;
    mf->alloced = size;
    mf->eof = 0;
    mf->offset = 0;
    mf->flush_pos = 0;
}


/*
 * Converts a FILE * to an mFILE *.
 * Use this for wrapper functions to turn external prototypes requring
 * FILE * as an argument into internal code using mFILE *.
 */
mFILE *mfreopen(const char *path, const char *mode_str, FILE *fp) {
    mFILE *mf;
    int r = 0, w = 0, a = 0, b = 0, mode = 0;

    /* Parse mode:
     * r = read file contents (if truncated => don't read)
     * w = write on close
     * a = position at end of buffer
     */
    if (strchr(mode_str, 'r'))
	r = 1, mode |= MF_READ;
    if (strchr(mode_str, 'w'))
	w = 1, mode |= MF_WRITE | MF_TRUNC;
    if (strchr(mode_str, 'a'))
	w = a = 1, mode |= MF_WRITE | MF_APPEND;
    if (strchr(mode_str, 'b'))
	b = 1, mode |= MF_BINARY;
    if (strchr(mode_str, '+')) {
        w = 1, mode |= MF_READ | MF_WRITE;
	if (a)
	    r = 1;
    }

    if (r) {
	mf = mfcreate(NULL, 0);
	if (!(mode & MF_TRUNC)) {
	    mf->data = mfload(fp, path, &mf->size, b);
	    mf->alloced = mf->size;
	    if (!a)
		fseek(fp, 0, SEEK_SET);
	}
    } else {
	/* Write - initialise the data structures */
	mf = mfcreate(NULL, 0);
    }
    mf->fp = fp;
    mf->mode = mode;
    
    if (a) {
	mf->flush_pos = mf->size;
	fseek(fp, 0, SEEK_END);
    }

    return mf;
}

/*
 * Opens a file. If we have read access (r or a+) then it loads the entire
 * file into memory. If We have write access then the pathname is stored.
 * We do not actually write until an mfclose, which then checks this pathname.
 */
mFILE *mfopen(const char *path, const char *mode) {
    FILE *fp;

    if (NULL == (fp = fopen(path, mode)))
	return NULL;
    return mfreopen(path, mode, fp);
}

/*
 * Closes an mFILE. If the filename is known (implying write access) then this
 * also writes the data to disk.
 *
 * Stdout is handled by calling mfflush which writes to stdout if appropriate.
 */
int mfclose(mFILE *mf) {
    if (!mf)
	return -1;

    mfflush(mf);

    if (mf->fp)
	fclose(mf->fp);

    mfdestroy(mf);

    return 0;
}

/*
 * Closes the file pointer contained within the mFILE without destroying
 * the in-memory data.
 */
int mfdetatch(mFILE *mf) {
    if (!mf)
	return -1;

    mfflush(mf);

    if (mf->fp) {
	fclose(mf->fp);
	mf->fp = NULL;
    }

    return 0;
}

/*
 * Destroys an mFILE structure but does not flush or close it
 */
int mfdestroy(mFILE *mf) {
    if (!mf)
	return -1;

    if (mf->data)
	free(mf->data);
    free(mf);

    return 0;
}

/*
 * Seek/tell functions. Nothing more than updating and reporting an
 * in-memory index. NB we can seek on stdin or stdout even provided we
 * haven't been flushing.
 */
int mfseek(mFILE *mf, long offset, int whence) {
    switch (whence) {
    case SEEK_SET:
	mf->offset = offset;
	break;
    case SEEK_CUR:
	mf->offset += offset;
	break;
    case SEEK_END:
	mf->offset = mf->size + offset;
	break;
    default:
	errno = EINVAL;
	return -1;
    }

    mf->eof = 0;
    return 0;
}

long mftell(mFILE *mf) {
    return mf->offset;
}

void mrewind(mFILE *mf) {
    mf->offset = 0;
    mf->eof = 0;
}

/*
 * mftruncate is not directly a translation of ftruncate as the latter
 * takes a file descriptor instead of a FILE *. It performs the analogous
 * role though.
 *
 * If offset is -1 then the file is truncated to be the current file
 * offset.
 */
void mftruncate(mFILE *mf, long offset) {
    mf->size = offset != -1 ? offset : mf->offset;
    if (mf->offset > mf->size)
	mf->offset = mf->size;
}

int mfeof(mFILE *mf) {
    return mf->eof;
}

/*
 * mFILE read/write functions. Basically these turn fread/fwrite syntax
 * into memcpy statements, with appropriate memory handling for writing.
 */
size_t mfread(void *ptr, size_t size, size_t nmemb, mFILE *mf) {
    size_t len;
    char *cptr = (char *)ptr;
    
    if (mf == m_channel[0]) init_mstdin();

    if (mf->size <= mf->offset)
	return 0;

    len = size * nmemb <= mf->size - mf->offset
	? size * nmemb
	: mf->size - mf->offset;
    if (!size)
	return 0;

    memcpy(cptr, &mf->data[mf->offset], len);
    mf->offset += len;
    cptr += len;
    
    if (len != size * nmemb) {
	mf->eof = 1;
    }

    return len / size;
}

size_t mfwrite(void *ptr, size_t size, size_t nmemb, mFILE *mf) {
    if (!(mf->mode & MF_WRITE))
	return 0;

    /* Append mode => forced all writes to end of file */
    if (mf->mode & MF_APPEND)
	mf->offset = mf->size;

    /* Make sure we have enough room */
    while (size * nmemb + mf->offset > mf->alloced) {
	mf->alloced = mf->alloced ? mf->alloced * 2 : 1024;
	mf->data = (void *)realloc(mf->data, mf->alloced);
    }

    /* Record where we need to reflush from */
    if (mf->offset < mf->flush_pos)
	mf->flush_pos = mf->offset;

    /* Copy the data over */
    memcpy(&mf->data[mf->offset], ptr, size * nmemb);
    mf->offset += size * nmemb;
    if (mf->size < mf->offset)
	mf->size = mf->offset;

    return nmemb;
}

int mfgetc(mFILE *mf) {
    if (mf == m_channel[0]) init_mstdin();
    if (mf->offset < mf->size) {
	return (unsigned char)mf->data[mf->offset++];
    }

    mf->eof = 1;
    return -1;
}

int mungetc(int c, mFILE *mf) {
    if (mf->offset > 0) {
	mf->data[--mf->offset] = c;
	return c;
    }
    
    mf->eof = 1;
    return -1;
}

char *mfgets(char *s, int size, mFILE *mf) {
    int i;

    if (mf == m_channel[0]) init_mstdin();
    *s = 0;
    for (i = 0; i < size-1;) {
	if (mf->offset < mf->size) {
	    s[i] = mf->data[mf->offset++];
	    if (s[i++] == '\n')
		break;
	} else {
	    mf->eof = 1;
	    break;
	}
    }

    s[i] = 0;
    return i ? s : NULL;
}

/*
 * Flushes an mFILE. If this is a real open of a file in write mode then
 * mFILE->fp will be set. We then write out any new data in mFILE since the
 * last flush. We cannot tell what may have been modified as we don't keep
 * track of that, so we typically rewrite out the entire file contents between
 * the last flush_pos and the end of file.
 *
 * For stderr/stdout we also reset the offsets so we cannot modify things
 * we've already output.
 */
int mfflush(mFILE *mf) {
    if (!mf->fp)
	return 0;

    /* FIXME: only do this when opened in write mode */
    if (mf == m_channel[1] || mf == m_channel[2]) {
	fwrite(mf->data + mf->flush_pos, 1, mf->size - mf->flush_pos, mf->fp);
	fflush(mf->fp);

	/* Stdout & stderr are non-seekable streams so throw away the data */
	mf->offset = mf->size = mf->flush_pos = 0;
    }

    /* only flush when opened in write mode */
    if (mf->mode & MF_WRITE) {
	if (mf->flush_pos < mf->size) {
	    fseek(mf->fp, mf->flush_pos, SEEK_SET);
	    fwrite(mf->data + mf->flush_pos, 1,
		   mf->size - mf->flush_pos, mf->fp);
	    fflush(mf->fp);
	}
	ftruncate(fileno(mf->fp), ftell(mf->fp));
	mf->flush_pos = mf->size;
    }

    return 0;
}

/*
 * A wrapper around vsprintf() to write to an mFILE. This also uses vflen() to
 * estimate how many additional bytes of storage will be required for the
 * vsprintf to work.
 */
int mfprintf(mFILE *mf, char *fmt, ...) {
    int ret;
    size_t est_length;
    va_list args;

    va_start(args, fmt);
    est_length = vflen(fmt, args);
    va_end(args);
    while (est_length + mf->offset > mf->alloced) {
	mf->alloced = mf->alloced ? mf->alloced * 2 : 1024;
	mf->data = (void *)realloc(mf->data, mf->alloced);
    }

    va_start(args, fmt);
    ret = vsprintf(&mf->data[mf->offset], fmt, args);
    va_end(args);

    if (ret > 0) {
	mf->offset += ret;
	if (mf->size < mf->offset)
	    mf->size = mf->offset;
    }

    if (mf->fp == stderr) {
	/* Auto-flush for stderr */
	mfflush(mf);
    }

    return ret;
}

/*
 * Converts an mFILE from binary to ascii mode by replacing all
 * cr-nl with nl.
 *
 * Primarily used on windows when we've uncompressed a binary file which
 * happens to be a text file (eg Experiment File). Previously we would have
 * seeked back to the start and used _setmode(fileno(fp), _O_TEXT).
 *
 * Side effect: resets offset and flush_pos back to the start.
 */
void mfascii(mFILE *mf) {
    size_t p1, p2;

    for (p1 = p2 = 1; p1 < mf->size; p1++, p2++) {
	if (mf->data[p1] == '\n' && mf->data[p1-1] == '\r') {
	    p2--; /* delete the \r */
	}
	mf->data[p2] = mf->data[p1];
    }
    mf->size = p2;

    mf->offset = mf->flush_pos = 0;
}
