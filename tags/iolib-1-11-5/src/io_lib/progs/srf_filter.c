/*
 * ======================================================================
 * This software has been created by Genome Research Limited (GRL).
 *
 * GRL hereby grants permission to use, copy, modify and distribute
 * this software and its documentation for non-commercial purposes
 * without fee at the user's own risk on the basis set out below.
 *
 * GRL neither undertakes nor accepts any duty whether contractual or
 * otherwise in connection with the software, its use or the use of
 * any derivative, and makes no representations or warranties, express
 * or implied, concerning the software, its suitability, fitness for
 * a particular purpose or non-infringement.
 *
 * In no event shall the authors of the software or GRL be responsible
 * or liable for any loss or damage whatsoever arising in any way
 * directly or indirectly out of the use of this software or its
 * derivatives, even if advised of the possibility of such damage.
 *
 * Our software can be freely distributed under the conditions set out
 * above, and must contain this copyright notice.
 * ======================================================================
 */

/*
 * This performs a linear (non-indexed) search for a trace in an SRF archive.
 *
 * It's not intended as a suitable production program or as a library of code
 * to use, but as a test and benchmark statistic.
 */

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/srf.h>

#define CR 13            /* Decimal code of Carriage Return char */
#define LF 10            /* Decimal code of Line Feed char */
#define EOF_MARKER 26    /* Decimal code of DOS end-of-file marker */
#define MAX_REC_LEN 1024 /* Maximum size of input buffer */

#define CHUNK_BASE (1 << 0)
#define CHUNK_CNF1 (1 << 1)
#define CHUNK_CNF4 (1 << 2)
#define CHUNK_SAMP (1 << 3)
#define CHUNK_SMP4 (1 << 4)
#define CHUNK_ALL (CHUNK_BASE | CHUNK_CNF1 | CHUNK_CNF4 | CHUNK_SAMP | CHUNK_SMP4)

#define NCHUNKS 6

#define CHUNK_BASE_STR "BASE"
#define CHUNK_CNF1_STR "CNF1"
#define CHUNK_CNF4_STR "CNF4"
#define CHUNK_SAMP_STR "SAMP"
#define CHUNK_SMP4_STR "SMP4"
#define CHUNK_ALL_STR  "ALL"

#define TYPE_PROC (1 << 0)
#define TYPE_SLXI (1 << 1)
#define TYPE_SLXN (1 << 2)
#define TYPE_0FAM (1 << 3)
#define TYPE_1CY3 (1 << 4)
#define TYPE_2TXR (1 << 5)
#define TYPE_3CY5 (1 << 6)
#define TYPE_ALL (TYPE_PROC | TYPE_SLXI | TYPE_SLXN | TYPE_0FAM | TYPE_1CY3 | TYPE_2TXR | TYPE_3CY5)

#define NTYPES 8

#define TYPE_PROC_STR "PROC"
#define TYPE_SLXI_STR "SLXI"
#define TYPE_SLXN_STR "SLXN"
#define TYPE_0FAM_STR "0FAM"
#define TYPE_1CY3_STR "1CY3"
#define TYPE_2TXR_STR "2TXR"
#define TYPE_3CY5_STR "3CY5"
#define TYPE_ALL_STR  "ALL"

/* ------------------------------------------------------------------------ */

/*
 * All the reads and prefixes will be used to filter reads in an srf archive.
 * For reads the match with the reads in the archive must be exact.  Prefixes
 * need only match to the beginning of the read, so they are a way to filter
 * several reads that, for example, belong to the same lane and tile.  For
 * example, if the reads are <Center>:<Run>:<Lane>:<Tile>:<X>:<Y>, then a
 * prefix can progress from center all the way to y coordinate in terms of
 * read specificity.
 */
typedef struct {
    int prefixes_size;
    int reads_size;
    char** prefixes;  /* Prefixes to filter on. */
    char** reads;     /* Reads to filter on. */
} read_filter_t;

/* ------------------------------------------------------------------------ */

/*
 * Print usage message to stderr and exit with the given \"code\".
 */
void usage(int code) {
    fprintf(stderr, "Usage: srf_filter [-c chunk_types] [-f read_filter] [-C] [-o] [-v] input(s) output\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -c chunk_types\n");
    fprintf(stderr, "              Chunk types to output given as a comma delimited list of types.\n");
    fprintf(stderr, "              The following types are allowed: \"ALL\", \"BASE\", \"CNF1\", \"CNF4\"\n");
    fprintf(stderr, "              \"SAMP\", \"SMP4\".\n");
    fprintf(stderr, "              The default is \"ALL\".\n");
    fprintf(stderr, "\n    -m mdata_types\n");
    fprintf(stderr, "              SAMP/SMP4 mdata types to output given as a comma delimited list of types.\n");
    fprintf(stderr, "              The following types are allowed: \"ALL\", \"PROC\", \"SLXI\", \"SLXN\"\n");
    fprintf(stderr, "              \"0FAM\", \"1CY3\", \"2TXR\", \"3CY5\".\n");
    fprintf(stderr, "              The default is \"ALL\".\n");
    fprintf(stderr, "\n    -f read_filter\n");
    fprintf(stderr, "              The filter to apply to reads in the archive.  If reads match the\n");
    fprintf(stderr, "              filter they are dumped.\n");
    fprintf(stderr, "              The filter takes the form of <name>=<value>, where <name> can be\n");
    fprintf(stderr, "              \"read\", \"prefix\", \"file\".\n");
    fprintf(stderr, "              If the <name> is \"read\" the value is represents the name of a\n");
    fprintf(stderr, "              read to dump.  Only the matching read will be dumped.\n");
    fprintf(stderr, "              If the <name> is \"prefix\" the value is represents the prefix of\n");
    fprintf(stderr, "              the name of the reads to dump.  Only the matching reads will be\n");
    fprintf(stderr, "              dumped.\n");
    fprintf(stderr, "              If the <name> is \"file\" the value is a file name where any\n");
    fprintf(stderr, "              number of \"read\" and \"prefix\" name value pairs can be included,\n");
    fprintf(stderr, "              one per line.\n");
    fprintf(stderr, "              The default is no filter, which means all reads are dumped.\n");
    fprintf(stderr, "\n    -b      exclude bad reads using readsFlags bitmask in data block header.\n");
    fprintf(stderr, "\n    -v      Print verbose messages.\n");
    fprintf(stderr, "\n");

    exit(code);
}

/*
 * Reads the lines from the fiven \"input\" file and creates the reads and prefixes
 * for the read filter.  The \"read_filter\" given is allocated.
 *
 * Returns 0 on success.  Errors cause the usage message and exit. 
 */
int read_filter_from_file(FILE *input, read_filter_t *read_filter)
{
    int   isNewline;              /* Boolean indicating we've read a CR or LF */
    long  lFileLen;               /* Length of file */
    long  lIndex;                 /* Index into cThisLine array */
    long  lLineCount;             /* Current line number */
    long  lLineLen;               /* Current line length */
    long  lStartPos;              /* Offset of start of current line */
    long  lTotalChars;            /* Total characters read */
    char  cThisLine[MAX_REC_LEN]; /* Contents of current line */
    char *cFile;                  /* Dynamically allocated buffer (entire file) */
    char *cThisPtr;               /* Pointer to current position in cFile */

    char *filter_type;
    char *prefix;
    char *read;

    fseek(input, 0L, SEEK_END);  /* Position to end of file */
    lFileLen = ftell(input);     /* Get file length */
    rewind(input);               /* Back to start of file */

    cFile = calloc(lFileLen + 1, sizeof(char));

    if(cFile == NULL )
	{
	    fprintf(stderr, "\nInsufficient memory to read file.\n");
	    return 0;
	}

    fread(cFile, lFileLen, 1, input); /* Read the entire file into cFile */

    lLineCount  = 0L;
    lTotalChars = 0L;

    cThisPtr    = cFile;              /* Point to beginning of array */

    while (*cThisPtr)                 /* Read until reaching null char */
	{
	    lIndex    = 0L;                 /* Reset counters and flags */
	    isNewline = 0;
	    lStartPos = lTotalChars;

	    while (*cThisPtr)               /* Read until reaching null char */
		{
		    if (!isNewline)               /* Haven't read a CR or LF yet */
			{
			    if (*cThisPtr == CR || *cThisPtr == LF) /* This char IS a CR or LF */
				isNewline = 1;                        /* Set flag */
			}

		    else if (*cThisPtr != CR && *cThisPtr != LF) /* Already found CR or LF */
			break;                                     /* Done with line */

		    /* Don't copy LS or CR */
		    if (*cThisPtr != CR && *cThisPtr != LF) {
			cThisLine[lIndex++] = *cThisPtr++; /* Add char to output and increment */
			++lTotalChars;
		    } else {
			cThisPtr++;
		    }

		} /* end while (*cThisPtr) */

	    cThisLine[lIndex] = '\0';     /* Terminate the string */
	    ++lLineCount;                 /* Increment the line counter */
	    lLineLen = strlen(cThisLine); /* Get length of line */

	    /* Find the one and only = in the string. */
	    if(strchr(cThisLine,'=') != NULL && (strchr(cThisLine,'=') == strrchr(cThisLine,'='))) {
		filter_type = strtok (cThisLine, "=");
	    } else {
		fprintf(stderr, "Baddly formatted read filter \"%s\".  Expected an \"=\" character in middle of filter.\n", cThisLine);
		usage(1);
	    }

	    if (!strcmp(filter_type, "prefix")) {
		prefix = strtok (NULL, "=");
		if(prefix == NULL) {
		    fprintf(stderr, "Bad prefix \"%s\" in read filter \"%s\".\n", prefix, cThisLine);
		    usage(1);
		} else {
		    ++(read_filter->prefixes_size);
		    read_filter->prefixes = (char**) realloc (read_filter->prefixes, read_filter->prefixes_size * sizeof(char *));
		    read_filter->prefixes[read_filter->prefixes_size - 1] =  (char*) calloc (strlen(prefix) + 1,sizeof(char));
		    strcpy(read_filter->prefixes[read_filter->prefixes_size - 1], prefix);
		}
	    } else if (!strcmp(filter_type, "read")) {
		read = strtok (NULL, "=");
		if(read == NULL) {
		    fprintf(stderr, "Bad read \"%s\" in read filter \"%s\".\n", read, cThisLine);
		    usage(1);
		} else {
		    ++(read_filter->reads_size);
		    read_filter->reads = (char**) realloc (read_filter->reads, read_filter->reads_size * sizeof(char *));
		    read_filter->reads[read_filter->reads_size - 1] =  (char*) calloc (strlen(read) + 1,sizeof(char));
		    strcpy(read_filter->reads[read_filter->reads_size - 1], read);
		}
	    } else {
		fprintf(stderr, "Unrecognized filter type \"%s\" given as part of read filter \"%s\".  The valid filter types are \"%s\".\n", filter_type, cThisLine, "prefix or read");
		usage(1);
	    }

	}

    free(cFile);
    return 0;
}

/*
 * Parse the given \"filter_value\" string for the filter value.
 *
 * Returns a read filter allocated on the heap which needs to be freed by the
 * calling function. Errors usually cause the usage message and exit.
 */
read_filter_t *get_read_filter(char *filter_value)
{
    char *filter_type = NULL;
    char *file_name = NULL;
    FILE *fp = NULL;
    char *prefix = NULL;
    char *read = NULL;

    /* Create read filter. */
    read_filter_t *read_filter = (read_filter_t *)calloc(1, sizeof(read_filter_t));
    if(read_filter == NULL) {
	return NULL;
    }
    read_filter->prefixes_size = 0;
    read_filter->reads_size = 0;

    /* Find the one and only = in the string. */
    if( strchr(filter_value,'=') != NULL &&
	(strchr(filter_value,'=') == strrchr(filter_value,'=')) ) {
	filter_type = strtok (filter_value,"=");
    } else {
        fprintf(stderr, "Baddly formatted read filter \"%s\".  Expected an \"=\" character in middle of filter.\n", filter_value);
        usage(1);
    }

    /* Check the string before the = is a valid filter type. */
    if(!strcmp("file", filter_type)) {

        /* Read the file. */
	file_name = strtok (NULL, "=");
	if(file_name == NULL) {
	    fprintf(stderr, "Bad file name \"%s\" in read filter \"%s\".\n", file_name, filter_value);
	    usage(1);
	}
	fp = fopen(file_name, "r");
	if(fp == NULL) {
            fprintf(stderr, "Bad file name \"%s\" in read filter \"%s\".\n", file_name, filter_value);
	    usage(1);
	}

	/* Read line by line. */
	if(read_filter_from_file(fp, read_filter)) {
	    fprintf(stderr, "Bad contents of file %s.\n", file_name);
	    usage(1);
	}

    } else if (!strcmp("prefix", filter_type)) {
	prefix = strtok (NULL, "=");
	if(prefix == NULL) {
	    fprintf(stderr, "Bad prefix \"%s\" in read filter \"%s\".\n", prefix, filter_value);
	    usage(1);
	} else {
	    ++(read_filter->prefixes_size);
	    read_filter->prefixes = (char**) malloc (read_filter->prefixes_size * sizeof(char *));
	    read_filter->prefixes[read_filter->prefixes_size - 1] =  (char*) calloc (strlen(prefix) + 1,sizeof(char));
	    strcpy(read_filter->prefixes[read_filter->prefixes_size - 1], prefix);
	}
    } else if (!strcmp("read", filter_type)) {
	read = strtok (NULL, "=");
	if(read == NULL) {
	    fprintf(stderr, "Bad read \"%s\" in read filter \"%s\".\n", read, filter_value);
	    usage(1);
	} else {
	    ++(read_filter->reads_size);
	    read_filter->reads = (char**) malloc (read_filter->reads_size * sizeof(char *));
	    read_filter->reads[read_filter->reads_size - 1] =  (char*) calloc (strlen(read) + 1, sizeof(char));
	    strcpy(read_filter->reads[read_filter->reads_size - 1], read);
	}
    } else {
	fprintf(stderr, "Unrecognized filter type \"%s\" given as part of read filter \"%s\".  The valid filter types are \"%s\".\n", filter_type, filter_value, "prefix, read, or file");
	usage(1);
    }

    return read_filter;
}

/*
 * Parse the comma delimited list of chunk types and put them in the single character \"mode\".
 *
 * Returns 0 on success.
 */
int get_chunk_types(char *arg, char *mode) {
    int num_allowed_types = NCHUNKS;
    char *allowed_str_types[] = {CHUNK_BASE_STR,CHUNK_CNF1_STR,CHUNK_CNF4_STR,CHUNK_SAMP_STR,CHUNK_SMP4_STR,CHUNK_ALL_STR};
    char allowed_types[] = {CHUNK_BASE,CHUNK_CNF1,CHUNK_CNF4,CHUNK_SAMP,CHUNK_SMP4,CHUNK_ALL};
    char *type;
    int i = 0;

    type = strtok (arg,",");
    while(type) {
	for(i = 0; i < num_allowed_types; i++) {
	    if(!strcmp(type, allowed_str_types[i]) && !(*mode & allowed_types[i])) {
		*mode += allowed_types[i];
		break;
	    }
	}
        type = strtok (NULL, ",");
    }

    return 0;
}

/*
 * Parse the comma delimited list of mdata types and put them in the single character \"mode\".
 *
 * Returns 0 on success.
 */
int get_mdata_types(char *arg, char *mode) {
    int num_allowed_types = NTYPES;
    char *allowed_str_types[] = {TYPE_PROC_STR, TYPE_SLXI_STR, TYPE_SLXN_STR, TYPE_0FAM_STR, TYPE_1CY3_STR, TYPE_2TXR_STR, TYPE_3CY5_STR, TYPE_ALL_STR};
    char allowed_types[] = {TYPE_PROC, TYPE_SLXI, TYPE_SLXN, TYPE_0FAM, TYPE_1CY3, TYPE_2TXR, TYPE_3CY5, TYPE_ALL};
    char *type;
    int i = 0;

    type = strtok (arg,",");
    while(type) {
	for(i = 0; i < num_allowed_types; i++) {
	    if(!strcmp(type, allowed_str_types[i]) && !(*mode & allowed_types[i])) {
		*mode += allowed_types[i];
		break;
	    }
	}
        type = strtok (NULL, ",");
    }

    return 0;
}

/*
 * Returns 1 is the read \"name\" matches any of the reads or prefixes in the \"read_filter\".
 */
int check_read_name(read_filter_t *read_filter, char *name) {
    int i;

    for(i = 0; i < read_filter->prefixes_size; i++) {
	if(name == strstr(name, read_filter->prefixes[i])) {
	    return 1;
	}
    }

    for(i = 0; i < read_filter->reads_size; i++) {
	if(!strcmp(name, read_filter->reads[i])) {
	    free(read_filter->reads[i]);
	    read_filter->reads[i] = read_filter->reads[read_filter->reads_size - 1];
	    read_filter->reads_size--;
	    return 1;
	}
    }

    return 0;
}

void dump_read_filter(read_filter_t *read_filter) {
    int i;

    printf("Read filter used:\n");

    for(i = 0; i < read_filter->prefixes_size; i++) {
	printf("\tprefix[%d] = %s\n", i, read_filter->prefixes[i]);
    }

    for(i = 0; i < read_filter->reads_size; i++) {
	printf("\tread[%d] = %s\n", i, read_filter->reads[i]);
    }
}

void dump_chunk_mode(char mode) {
    printf("mode: %d.\n", mode);

    if(mode & CHUNK_BASE) {
	printf("BASE chunk required.\n");
    }
    if(mode & CHUNK_CNF1) {
	printf("CNF1 chunk required.\n");
    }
    if(mode & CHUNK_CNF4) {
	printf("CNF4 chunk required.\n");
    }
    if(mode & CHUNK_SAMP) {
	printf("SAMP chunk required.\n");
    }
    if(mode & CHUNK_SMP4) {
	printf("SMP4 chunk required.\n");
    }
}

void dump_mdata_mode(char mode) {
    printf("mode: %d.\n", mode);

    if(mode & TYPE_PROC) {
	printf("Illumina PROC SAMP/SMP4 chunk required.\n");
    }
    if(mode & TYPE_SLXI) {
	printf("Illumina SLXI SAMP/SMP4 chunk required.\n");
    }
    if(mode & TYPE_SLXN) {
	printf("Illumina SLXN SAMP/SMP4 chunk required.\n");
    }
    if(mode & TYPE_0FAM) {
	printf("Solid 0FAM SAMP/SMP4 chunk required.\n");
    }
    if(mode & TYPE_1CY3) {
	printf("Solid 1CY3 SAMP/SMP4 chunk required.\n");
    }
    if(mode & TYPE_2TXR) {
	printf("Solid 2TXR SAMP/SMP4 chunk required.\n");
    }
    if(mode & TYPE_3CY5) {
	printf("Solid 3CY3 SAMP/SMP4 chunk required.\n");
    }
}

/*
 * Given the input archive name (input), the output archive name (output),
 * the chunk types to output (chunk_mode) and some other parameters such as
 * the read filter generates a filtered srf file.
 *
 * Note the generated srf file is NOT indexed
 *
 * Returns 0 on success.
 */
int srf_filter(char *input, srf_t *out_srf, char chunk_mode, char mdata_mode, int filter_mode, read_filter_t *read_filter, int read_mask) {
    srf_t *in_srf;
    char name[1024];

    if (NULL == (in_srf = srf_open(input, "rb"))) {
	perror(input);
        return 1;
    }

    do {
      int type;
      ztr_chunk_t *chunk;

      switch(type = srf_next_block_type(in_srf)) {
	case SRFB_CONTAINER:
          if (0 != srf_read_cont_hdr(in_srf, &in_srf->ch)) {
            fprintf(stderr, "Error reading container header.\nExiting.\n");
	    exit(1);
          }
          if (0 != srf_write_cont_hdr(out_srf, &in_srf->ch)) {
            fprintf(stderr, "Error writing container header.\nExiting.\n");
	    exit(1);
          }
          break;

        case SRFB_XML:
          if (0 != srf_read_xml(in_srf, &in_srf->xml)) {
            fprintf(stderr, "Error reading XML.\nExiting.\n");
	    exit(1);
          }
          if (0 != srf_write_xml(out_srf, &in_srf->xml)) {
            fprintf(stderr, "Error writing XML.\nExiting.\n");
	    exit(1);
          }
          break;

	case SRFB_TRACE_HEADER:
          if (0 != srf_read_trace_hdr(in_srf, &in_srf->th)) {
            fprintf(stderr, "Error reading trace header.\nExiting.\n");
	    exit(1);
          }

#if 1
          if(chunk_mode == CHUNK_ALL && mdata_mode == TYPE_ALL ){
            if (0 != srf_write_trace_hdr(out_srf, &in_srf->th)) {
              fprintf(stderr, "Error writing trace header.\nExiting.\n");
              exit(1);
            }
            break;
          }
#endif          
      
          /* Decode ZTR chunks in the header */
          if (in_srf->mf)
              mfdestroy(in_srf->mf);

          in_srf->mf = mfcreate(NULL, 0);
          if (in_srf->th.trace_hdr_size)
            mfwrite(in_srf->th.trace_hdr, 1, in_srf->th.trace_hdr_size, in_srf->mf);
          if (in_srf->ztr)
            delete_ztr(in_srf->ztr);
          mrewind(in_srf->mf);

          /* create the trace header data */
          mFILE *mf = mfcreate(NULL, 0);

          if (NULL != (in_srf->ztr = partial_decode_ztr(in_srf, in_srf->mf, NULL))) {
            in_srf->mf_pos = mftell(in_srf->mf);
            mfseek(in_srf->mf, 0, SEEK_END);
            in_srf->mf_end = mftell(in_srf->mf);

            mfseek(in_srf->mf, 0, SEEK_SET);
            mfwrite(in_srf->mf->data, 1, sizeof(ztr_header_t), mf);
            mfseek(in_srf->mf, sizeof(ztr_header_t), SEEK_CUR);

            int pos = mftell(in_srf->mf);
            while (chunk = ztr_read_chunk_hdr(in_srf->mf)) {
              char *key = ztr_lookup_mdata_value(in_srf->ztr, chunk, "TYPE");
              int flag = 0;

              /* filter on chunk type */
              switch (chunk->type) {
	      case ZTR_TYPE_BASE:
                if (chunk_mode & CHUNK_BASE)
                  flag = 1;
                break;
	      case ZTR_TYPE_CNF1:
                if (chunk_mode & CHUNK_CNF1)
                  flag = 1;
                break;
	      case ZTR_TYPE_CNF4:
                if (chunk_mode & CHUNK_CNF4)
                  flag = 1;
                break;
	      case ZTR_TYPE_SAMP:
                if (chunk_mode & CHUNK_SAMP) {
                  if (mdata_mode == TYPE_ALL)
                    flag = 1;
                  if ((mdata_mode & TYPE_0FAM) && (key && 0 == strcmp(key, "0FAM")))
                    flag = 1;
                  if ((mdata_mode & TYPE_1CY3) && (key && 0 == strcmp(key, "1CY3")))
                    flag = 1;
                  if ((mdata_mode & TYPE_2TXR) && (key && 0 == strcmp(key, "2TXR")))
                    flag = 1;
                  if ((mdata_mode & TYPE_3CY5) && (key && 0 == strcmp(key, "3CY5")))
                    flag = 1;
                }
	        break;
              case ZTR_TYPE_SMP4:
                if (chunk_mode & CHUNK_SMP4) {
                  if (mdata_mode == TYPE_ALL)
                    flag = 1;
                  if ((mdata_mode & TYPE_PROC) && (NULL == key || 0 == strcmp(key, "PROC")))
                    flag = 1;
                  if ((mdata_mode & TYPE_SLXI) && (key && 0 == strcmp(key, "SLXI")))
                    flag = 1;
                  if ((mdata_mode & TYPE_SLXN) && (key && 0 == strcmp(key, "SLXN")))
                    flag = 1;
                }
                break;
              default:
                flag = 1;
                break;
              }

              if (flag)
                mfwrite(in_srf->mf->data+pos, 1, (4+4+chunk->mdlength+4+chunk->dlength), mf);
              mfseek(in_srf->mf, chunk->dlength, SEEK_CUR);
              pos = mftell(in_srf->mf);

   	      if (chunk->mdata)
	        xfree(chunk->mdata);
  	      xfree(chunk);
            }

          } else {
            /* Maybe not enough to decode or no headerBlob. */
            /* So delay until decoding the body. */
            in_srf->mf_pos = in_srf->mf_end = 0;
          }

          /* construct the new trace header */
          srf_trace_hdr_t th;
          srf_construct_trace_hdr(&th, in_srf->th.id_prefix, (unsigned char *)mf->data, mftell(mf));
          if (0 != srf_write_trace_hdr(out_srf, &th)) {
            fprintf(stderr, "Error writing trace header.\nExiting.\n");
            exit(1);
          }

	  mfdestroy(mf);

          break;

	case SRFB_TRACE_BODY: {
          srf_trace_body_t old_tb;
          ztr_t *ztr_tmp;

          if (0 != srf_read_trace_body(in_srf, &old_tb, 0)) {
            fprintf(stderr, "Error reading trace body.\nExiting.\n");
	    exit(1);
          }
          
          if (-1 == construct_trace_name(in_srf->th.id_prefix,
                                         (unsigned char *)old_tb.read_id,
                                         old_tb.read_id_length,
                                         name, 512)) {
            fprintf(stderr, "Error constructing trace name.\nExiting.\n");
	    exit(1);
          }

          if (old_tb.flags & read_mask)
   	    break;

          if(filter_mode && !check_read_name(read_filter, name))
	    break;
          
#if 1
          if(chunk_mode == CHUNK_ALL && mdata_mode == TYPE_ALL ){
            if (0 != srf_write_trace_body(out_srf, &old_tb)) {
              fprintf(stderr, "Error writing trace body.\nExiting.\n");
              exit(1);
            }
            break;
          }
#endif          

          if (!in_srf->mf) {
            fprintf(stderr, "Error reading trace body.\nExiting.\n");
	    exit(1);
          }

          mfseek(in_srf->mf, in_srf->mf_end, SEEK_SET);
          if (old_tb.trace_size) {
              mfwrite(old_tb.trace, 1, old_tb.trace_size, in_srf->mf);
              free(old_tb.trace);
              old_tb.trace = NULL;
          }
          mftruncate(in_srf->mf, mftell(in_srf->mf));
          mfseek(in_srf->mf, in_srf->mf_pos, SEEK_SET);

          if (in_srf->ztr)
              ztr_tmp = ztr_dup(in_srf->ztr); /* inefficient, but simple */
          else
              ztr_tmp = NULL;

          if (NULL != partial_decode_ztr(in_srf, in_srf->mf, ztr_tmp)) {

            /* create the trace body data */
            mFILE *mf = mfcreate(NULL, 0);

            /* include the ztr header if it wasn't in the trace header block */
            if( !in_srf->mf_pos ){
              mfseek(in_srf->mf, 0, SEEK_SET);
              mfwrite(in_srf->mf->data, 1, sizeof(ztr_header_t), mf);
              mfseek(in_srf->mf, sizeof(ztr_header_t), SEEK_CUR);
            }else{
              mfseek(in_srf->mf, in_srf->mf_pos, SEEK_SET);
            }

            int pos = mftell(in_srf->mf);
            while (chunk = ztr_read_chunk_hdr(in_srf->mf)) {
              char *key = ztr_lookup_mdata_value(in_srf->ztr, chunk, "TYPE");
              int flag = 0;

              /* filter on chunk type */
              switch (chunk->type) {
	      case ZTR_TYPE_BASE:
                if (chunk_mode & CHUNK_BASE)
                  flag = 1;
                break;
	      case ZTR_TYPE_CNF1:
                if (chunk_mode & CHUNK_CNF1)
                  flag = 1;
                break;
	      case ZTR_TYPE_CNF4:
                if (chunk_mode & CHUNK_CNF4)
                  flag = 1;
                break;
	      case ZTR_TYPE_SAMP:
                if (chunk_mode & CHUNK_SAMP) {
                  if (mdata_mode == TYPE_ALL)
                    flag = 1;
                  if ((mdata_mode & TYPE_0FAM) && (key && 0 == strcmp(key, "0FAM")))
                    flag = 1;
                  if ((mdata_mode & TYPE_1CY3) && (key && 0 == strcmp(key, "1CY3")))
                    flag = 1;
                  if ((mdata_mode & TYPE_2TXR) && (key && 0 == strcmp(key, "2TXR")))
                    flag = 1;
                  if ((mdata_mode & TYPE_3CY5) && (key && 0 == strcmp(key, "3CY5")))
                    flag = 1;
                }
	        break;
              case ZTR_TYPE_SMP4:
                if (chunk_mode & CHUNK_SMP4) {
                  if (mdata_mode == TYPE_ALL)
                    flag = 1;
                  if ((mdata_mode & TYPE_PROC) && (NULL == key || 0 == strcmp(key, "PROC")))
                    flag = 1;
                  if ((mdata_mode & TYPE_SLXI) && (key && 0 == strcmp(key, "SLXI")))
                    flag = 1;
                  if ((mdata_mode & TYPE_SLXN) && (key && 0 == strcmp(key, "SLXN")))
                    flag = 1;
                }
                break;
              default:
                flag = 1;
                break;
              }

              if (flag)
                mfwrite(in_srf->mf->data+pos, 1, (4+4+chunk->mdlength+4+chunk->dlength), mf);
              mfseek(in_srf->mf, chunk->dlength, SEEK_CUR);
              pos = mftell(in_srf->mf);

   	      if (chunk->mdata)
	        xfree(chunk->mdata);
  	      xfree(chunk);
            }

            /* construct the new trace body */
            srf_trace_body_t new_tb;
            srf_construct_trace_body(&new_tb, name+strlen(in_srf->th.id_prefix), -1, mf->data, mf->size, old_tb.flags);

            if (0 != srf_write_trace_body(out_srf, &new_tb)) {
              fprintf(stderr, "Error writing trace body.\nExiting.\n");
              exit(1);
            }

	    mfdestroy(mf);
          }

	  if( ztr_tmp )
	      delete_ztr(ztr_tmp);

          break;
        }
      }

      if( type == -1 || type == SRFB_INDEX || type == SRFB_NULL_INDEX )
          break;

    } while (1);

    srf_destroy(in_srf, 1);
    return 0;
}

/* ------------------------------------------------------------------------ */

/*
 * Main method.
 */
int main(int argc, char **argv) {
  int nfiles, ifile;
    int filter_mode = 0;
    char *input = NULL;
    char *output = NULL;
    read_filter_t *read_filter = NULL;
    char *filter_value = NULL;

    int c;
    int errflg = 0;
    extern char *optarg;
    extern int optind, optopt;

    char chunk_mode = CHUNK_ALL;
    char mdata_mode = TYPE_ALL;
    int verbose = 0;
    int read_mask = 0;

    srf_t *srf = NULL;

    while ((c = getopt(argc, argv, ":c:m:f:vb")) != -1) {
        switch (c) {
        case 'c':
	    chunk_mode = 0;
	    if(get_chunk_types(optarg, &chunk_mode) || !chunk_mode) {
                fprintf(stderr,
			"Invalid value \"%s\" given to option -%c.\n", optarg, c);
		errflg++;
	    }
	    break;
        case 'm':
	    mdata_mode = 0;
	    if(get_mdata_types(optarg, &mdata_mode) || !mdata_mode) {
                fprintf(stderr,
			"Invalid value \"%s\" given to option -%c.\n", optarg, c);
		errflg++;
	    }
	    break;
        case 'f':
  	    filter_mode++;
	    filter_value = optarg;
            break;
        case 'v':
	    verbose++;
            break;
        case 'b':
	    read_mask = SRF_READ_FLAG_BAD_MASK;
            break;
        case ':':       /* -? without operand */
            fprintf(stderr,
                    "Option -%c requires an operand\n", optopt);
            errflg++;
            break;
        case '?':
            fprintf(stderr,
                    "Unrecognised option: -%c\n", optopt);
            errflg++;
        }
    }

    if (errflg) {
	usage(1);
    }

    nfiles = (argc-optind);
    if( nfiles < 2 ){
        fprintf(stderr, "Please specify input archive name(s) and an output archive name.\n");
        usage(1);
    }
    output = argv[optind+nfiles-1];
    nfiles--;
    
    if(filter_mode) {
	read_filter = get_read_filter(filter_value);
	if(verbose) {
	    dump_read_filter(read_filter);
	}
    }

    if(chunk_mode && verbose) {
	dump_chunk_mode(chunk_mode);
    }

    if(mdata_mode && verbose) {
	dump_mdata_mode(mdata_mode);
    }

    if (NULL == (srf = srf_open(output, "wb"))) {
        perror(output);
        return 1;
    }
    
    for (ifile=0; ifile<nfiles; ifile++) {
        input = argv[optind+ifile];
        printf("Reading archive %s.\n", input);

        srf_filter(input, srf, chunk_mode, mdata_mode, filter_mode, read_filter, read_mask);
    }

    if(NULL != srf)
        srf_destroy(srf, 1);

    return 0;
}
