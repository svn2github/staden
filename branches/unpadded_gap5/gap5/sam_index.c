#include <staden_config.h>
#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tg_gio.h"
#include "tg_struct.h"
#include "tg_index_common.h"
#include "sam_index.h"

#include <staden_config.h>
#ifdef HAVE_SAMTOOLS

#define _IOLIB 2
#include "bam.h"
#include "sam.h"
#include "sam_header.h"
#include "faidx.h"
#include "bam_maqcns.h"
#include "depad_seq_tree.h"

/*
 * Uncomment this if you want sam auxillary tags to be added as tag in
 * gap5.
 */
//#define SAM_AUX_AS_TAG

typedef struct {
    bam1_t *b;
    char *seq;
    char *conf;
    int seq_len;
    int alloc_len;
    int mqual;
    int pos;
    int left;
    tg_rec rec;
} bam_seq_t;

typedef struct rec_list {
    tg_rec rec;
    struct rec_list *next;
} rec_list_t;

typedef struct {
    GapIO *io;
    const char *fn;
    bam_seq_t *seqs;
    int nseq;
    int max_seq;
    rec_list_t *rec_head;
    rec_list_t *rec_tail;
    HacheTable *pair;
    HacheTable *libs;
    contig_t *c;
    int n_inserts;
    int npads;
    int count;
    int skip;
    bam_header_t *header;
    tg_args *a;
    struct PAD_COUNT *tree; /* re-padding */
    int last_tid;
    void *rg2pl_hash;
} bam_io_t;


/*
 * Converts a string from the read-group PL field to a constant STECH_*
 * value for sequencing technology.
 * str may be null in which case STECH_UNKNOWN is returned.
 */
int stech_str2int(const char *str) {
    if (!str)
	return STECH_UNKNOWN;

    if (strcasecmp(str, "ILLUMINA") == 0)
	return STECH_SOLEXA;
    if (strcasecmp(str, "SOLEXA") == 0)
	return STECH_SOLEXA;

    if (strcasecmp(str, "ABI") == 0)
	return STECH_SANGER;
    if (strcasecmp(str, "CAPILLARY") == 0)
	return STECH_SANGER;
    if (strcasecmp(str, "SANGER") == 0)
	return STECH_SANGER;

    if (strcasecmp(str, "454") == 0)
	return STECH_454;
    if (strcasecmp(str, "LS454") == 0)
	return STECH_454;

    if (strcasecmp(str, "SOLID") == 0)
	return STECH_SOLID;

    return STECH_UNKNOWN;
}

/*
 * Attempts to guess a sequencing technology based on sequence name.
 * This is far from ideal and is also a bit Sanger Institute centric.
 */
int stech_guess_by_name(char *str) {
    size_t l;
    char *cp;
    int colons;

    if (!str || !*str)
	return STECH_UNKNOWN;

    l = strlen(str);

    /* 454 follow a rigid pattern, [0-9A-Z]{14}, first 7 being the plate */
    if (l == 14) {
	int i;
	for (i = 0; i < l; i++) {
	    if (!isalnum(str[i]))
		break;
	}
	if (i == l)
	    return STECH_454;
    }

    /* SOLID appear to start VAB_ (vendor = AB) */
    if (strncmp(str, "VAB_", 4) == 0)
	return STECH_SOLID;

    /*
     * Illumina tend to start with machine-name followed by lane, run,
     * tile, x, y and possibly #index. We look for 4 colons or for
     * machine name IL[0-9]+.
     */
    if (strncmp(str, "IL", 2) == 0 && isdigit(str[2]))
	return STECH_SOLEXA;

    colons = 0;
    cp = str;
    do {
	cp = strchr(cp, ':');
	if (cp) {
	    colons++;
	    cp++;
	}
    } while(cp);

    if (colons == 4) {
	return STECH_SOLEXA;
    }

    
    /*
     * Sanger capillary sequences tend to be template_name.[pq][0-9]k.
     * Very sanger specific, but there's just too much variation.
     */
    cp = strchr(str, '.');
    if (cp) {
	if (cp[1] == 'p' || cp[1] == 'q') {
	    if (isdigit(cp[2])) {
		if (cp[3] == 'k') {
		    return STECH_SANGER;
		}
	    }
	}
    }

    return STECH_UNKNOWN;
}


int bio_extend_seq(bam_io_t *bio, int snum, char base, int conf);

/*
 * Allocates a new bam_seq_t entry in bam_io_t struct.
 * We use a single indexed array counting from 0 to N-1 representing the
 * N active sequences in a bam_pileup1_t struct. This may mean we have
 * O(D^2) complexity for depth D, so if we find this becomes a bottleneck
 * then we can replace the bam_seq_t array with an ordered linked list
 * instead.
 *
 * Returns the sequence index on success
 *         -1 on failure
 */
int bio_new_seq(bam_io_t *bio, const bam_pileup1_t *p, int pos) {
    int i;
    bam_seq_t *s;
    rec_list_t *tmp;

    //printf("New seq %d at %d %s\n", bio->nseq, pos, bam1_qname(p->b));

    if (bio->nseq == bio->max_seq) {
	bio->max_seq = bio->max_seq ? bio->max_seq * 2 : 256;
	bio->seqs = (bam_seq_t *)realloc(bio->seqs,
					 bio->max_seq * sizeof(bam_seq_t));
	if (!bio->seqs)
	    return -1;
    }
    
    s = &bio->seqs[bio->nseq];
    s->alloc_len = p->b->core.l_qseq+10;
    s->seq = (char *)malloc((int)(s->alloc_len * 1.2));
    s->conf = (char *)malloc((int)(s->alloc_len * 1.2));
    for (i = 0; i < p->qpos; i++) {
	s->seq[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), i)];
	s->conf[i] = bam1_qual(p->b)[i];
    }
    s->seq_len = i;
    s->pos = pos - i;
    s->b = p->b;
    s->left = i;

    s->rec = bio->rec_head->rec;
    tmp = bio->rec_head;
    bio->rec_head = bio->rec_head->next;
    xfree(tmp);

    return bio->nseq++;
}

typedef union {
    char  *s;
    int    i;
    float  f;
    double d;
} bam_aux_t;

/*
 * An iterator on bam aux fields. NB: This code is not reentrant or multi-
 * thread capable. The values returned are valid until the next call to
 * this function.
 * key:  points to an array of 2 characters (eg "RG", "NM")
 * type: points to an address of 1 character (eg 'Z', 'i')
 * val:  points to an address of a bam_aux_t union.
 *
 * Pass in *iter_handle as NULL to initialise the search and then
 * pass in the modified value on each subsequent call to continue the search.
 *
 * Returns 0 if the next value is valid, setting key, type and val.
 *        -1 when no more found.
 */
int bam_aux_iter(bam1_t *b, char **iter_handle,
		 char *key, char *type, bam_aux_t *val) {
    char *s;

    if (!iter_handle || !*iter_handle) {
	s = (char *)bam1_aux(b);
    } else {
	s = *iter_handle;
    }

    if ((uint8_t *)s >= b->data + b->data_len)
	return -1;

    key[0] = s[0];
    key[1] = s[1];
    
    switch (s[2]) {
    case 'A':
	if (type) *type = 'A';
	if (val) val->i = *(s+3);
	s+=4;
	break;

    case 'C':
	if (type) *type = 'i';
	if (val) val->i = *(uint8_t *)(s+3);
	s+=4;
	break;

    case 'c':
	if (type) *type = 'i';
	if (val) val->i = *(int8_t *)(s+3);
	s+=4;
	break;

    case 'S':
	if (type) *type = 'i';
	if (val) {
	    char tmp[2]; /* word aligned data */
	    tmp[0] = s[3]; tmp[1] = s[4];
	    val->i = *(uint16_t *)tmp;
	}
	s+=5;
	break;

    case 's':
	if (type) *type = 'i';
	if (val) {
	    char tmp[2]; /* word aligned data */
	    tmp[0] = s[3]; tmp[1] = s[4];
	    val->i = *(int16_t *)tmp;
	}
	s+=5;
	break;

    case 'I':
	if (type) *type = 'i';
	if (val) {
	    char tmp[4]; /* word aligned data */
	    tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
	    val->i = *(uint32_t *)tmp;
	}
	s+=7;
	break;

    case 'i':
	if (type) *type = 'i';
	if (val) {
	    char tmp[4]; /* word aligned data */
	    tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
	    val->i = *(int32_t *)tmp;
	}
	s+=7;
	break;

    case 'f':
	if (type) *type = 'f';
	if (val) memcpy(&val->f, s+3, 4);
	s+=7;
	break;

    case 'd':
	if (type) *type = 'd';
	if (val) memcpy(&val->d, s+3, 8);
	s+=11;
	break;

    case 'Z': case 'H':
	if (type) *type = s[2];
	s+=3;
	if (val) val->s = s;
	while (*s++);
	break;

    default:
	fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	return -1;
    }

    if (iter_handle)
	*iter_handle = s;

    return 0;
}

/*
 * Searches for 'key' in the bam auxillary tags.
 *
 * If found, key type and value are filled out in the supplied
 * 'type' and 'val' pointers. These may be supplied as NULL if the
 * caller simply wishes to test for the presence of a key instead.
 *
 * Returns 0 if found
 *        -1 if not
 */
int bam_aux_find(bam1_t *b, char *key, char *type, bam_aux_t *val) {
    char *h = NULL;
    char k[2];

    while (0 == bam_aux_iter(b, &h, k, type, val)) {
	if (k[0] == key[0] && k[1] == key[1])
	    return 0;
    }

    return -1;
}

#endif /* HAVE_SAMTOOLS */

char *sam_aux_stringify_old(char *s, int len) {
    static char str[8192];
    char *cp = str, *s_end = s+len;
    int first = 1;

    //write(2, s, (int)(b->data + b->data_len - (uint8_t *)s));

    while (s < s_end) {
	if (first)
	    first = 0;
	else
	    *cp++ = '\t';

	switch (s[2]) {
	case 'A':
	    cp += sprintf(cp, "%c%c:A:%c", s[0], s[1], *(s+3));
	    s+=4;
	    break;

	case 'C':
	    cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint8_t *)(s+3));
	    s+=4;
	    break;

	case 'c':
	    cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int8_t *)(s+3));
	    s+=4;
	    break;

	case 'S':
	    {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 's':
	    {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 'I':
	    {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'i':
	    {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'f':
	    {
		float f;
		memcpy(&f, s+3, 4);
		cp += sprintf(cp, "%c%c:f:%f", s[0], s[1], f);
	    }
	    s+=7;
	    break;

	case 'd':
	    {
		double d;
		memcpy(&d, s+3, 8);
		cp += sprintf(cp, "%c%c:d:%f", s[0], s[1], d);
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    cp += sprintf(cp, "%c%c:%c:%s", s[0], s[1], s[2], s+3);
	    s+=3;
	    while (*s++);
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *cp = 0;

    return str;
}

static char *append_int(char *cp, int i) {
    int j, k = 0;

    if (i < 0) {
	*cp++ = '-';
	i = -i;
    } else if (i == 0) {
	*cp++ = '0';
	return cp;
    }

    if (i < 1000)
	goto b1;
    if (i < 100000)
	goto b2;
    if (i < 100000000)
	goto b3;

    j = i / 1000000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 1000000000;

    j = i / 100000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 100000000;
    
 b3:
    j = i / 10000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 10000000;
    
    j = i / 1000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 1000000;
    
    j = i / 100000;
    if (j || k) *cp++ = j + '0', k=1, i %= 100000;
    
 b2:
    j = i / 10000;
    if (j || k) *cp++ = j + '0', k=1, i %= 10000;

    j = i / 1000;
    if (j || k) *cp++ = j + '0', k=1, i %= 1000;

 b1:
    j = i / 100;
    if (j || k) *cp++ = j + '0', k=1, i %= 100;

    j = i / 10;
    if (j || k) *cp++ = j + '0', k=1, i %= 10;

    if (i || k) *cp++ = i + '0';

    return cp;
}

#define APPEND_FMT(fmt) \
    *cp++ = s[0]; \
    *cp++ = s[1]; \
    *cp++ = ':'; \
    *cp++ = (fmt); \
    *cp++ = ':';

/* Experiments to speed up printf */
#if 0

/* A basic printf without a lot of the %-10 or %5.5 etc field width
 * specifiers.
 * It is, however, much faster.
 */
#define EVEN_SIMPLER
int fast_fprintf(FILE *fp, char *fmt, ...) {
    char *cp, c;
    long l;
    int i;
    double d; 
    char tmp[128], *ct;
    va_list ap;

    va_start(ap, fmt);

    for (cp = fmt; *cp; cp++) {
	switch(*cp) {

	/* A format specifier */
	case '%': {
	    char *endp;
	    long conv_len1=0, conv_len2=0, conv_len=0;
	    signed int arg_size;

	    cp++;
#ifndef EVEN_SIMPLER
	    /* Firstly, strip the modifier flags (+-#0 and [space]) */
	    for(; c=*cp;) {
		if ('#' == c)
		    ;
		else if ('-' == c || '+' == c || ' ' == c)
		    ;
		else
		    break;
	    }

	    /* Width specifier */
	    if (*cp >= '0' && *cp <= '9') {
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len = conv_len1 = l;
	    } else if (*cp == '*') {
		conv_len = conv_len1 = (int)va_arg(ap, int);
		cp++;
	    }
#endif

	    /* Precision specifier */
	    if ('.' == *cp) {
		cp++;
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len2 = l;
		if (*cp == '*') {
		    conv_len2 = (int)va_arg(ap, int);
		    cp++;
		}
		conv_len = MAX(conv_len1, conv_len2);
	    }

#ifndef EVEN_SIMPLER
	    /* Short/long identifier */
	    if ('h' == *cp) {
		arg_size = -1; /* short */
		cp++;
	    } else if ('l' == *cp) {
		arg_size = 1; /* long */
		cp++;
	    } else {
		arg_size = 0; /* int */
	    }
#endif

	    /* The actual type */
	    switch (*cp) {
	    case '%':
		/*
		 * Not real ANSI I suspect, but we'll allow for the
		 * completely daft "%10%" example.
		 */
		putc('%', fp);
		break;

	    case 'd':
	    case 'i':
	    case 'u':
		/* Remember: char and short are sent as int on the stack */
		if (arg_size == -1)
		    l = (long)va_arg(ap, int);
		else if (arg_size == 1)
		    l = va_arg(ap, long); 
		else 
		    l = (long)va_arg(ap, int);

		ct = append_int(tmp, (int)l);
		*ct++ = 0;
		fputs(tmp, fp);
		break;

	    case 'c':
		i = va_arg(ap, int);
		putc(i, fp);
		break;

	    case 'f':
		d = va_arg(ap, double);
		fprintf(fp, "%f", d);
		break;

	    case 'e':
	    case 'E':
	    case 'g':
	    case 'G':
		d = va_arg(ap, double);
		fprintf(fp, "%g", d);
		break;

	    case 'p':
	    case 'x':
	    case 'X':
		l = (long)va_arg(ap, void *);
		puts("TODO");
		break;

	    case 'n':
		/* produces no output */
		break;

	    case 's': {
		char *s = (char *)va_arg(ap, char *);
		if (conv_len2) {
		    fwrite(s, conv_len2, 1, fp);
		} else
		    fputs(s, fp);
		break;
	    }

	    default:
		/* wchar_t types of 'C' and 'S' aren't supported */
		puts("TODO");
	    }
	    
	}

	case '\0':
	    break;

	default:
	    putc(*cp, fp);
	}
    }

    va_end(ap);
    
    return 0;
}

/* Max 8k of output */
int faster_fprintf(FILE *fp, char *fmt, ...) {
    char *cp, c;
    long l;
    int i;
    double d; 
    char tmp[128], *ct;
    va_list ap;
    char max_line[8192], *out = max_line;

    va_start(ap, fmt);

    for (cp = fmt; *cp; cp++) {
	switch(*cp) {

	/* A format specifier */
	case '%': {
	    char *endp;
	    long conv_len1=0, conv_len2=0, conv_len=0;
	    signed int arg_size;

	    cp++;

#ifndef EVEN_SIMPLER
	    /* Firstly, strip the modifier flags (+-#0 and [space]) */
	    for(; c=*cp;) {
		if ('#' == c)
		    ;
		else if ('-' == c || '+' == c || ' ' == c)
		    ;
		else
		    break;
	    }

	    /* Width specifier */
	    if (*cp >= '0' && *cp <= '9') {
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len = conv_len1 = l;
	    } else if (*cp == '*') {
		conv_len = conv_len1 = (int)va_arg(ap, int);
		cp++;
	    }
#endif

	    /* Precision specifier */
	    if ('.' == *cp) {
		cp++;
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len2 = l;
		if (*cp == '*') {
		    conv_len2 = (int)va_arg(ap, int);
		    cp++;
		}
		conv_len = MAX(conv_len1, conv_len2);
	    }

#ifndef EVEN_SIMPLER
	    /* Short/long identifier */
	    if ('h' == *cp) {
		arg_size = -1; /* short */
		cp++;
	    } else if ('l' == *cp) {
		arg_size = 1; /* long */
		cp++;
	    } else {
		arg_size = 0; /* int */
	    }
#endif

	    /* The actual type */
	    switch (*cp) {
	    case '%':
		/*
		 * Not real ANSI I suspect, but we'll allow for the
		 * completely daft "%10%" example.
		 */
		*out++ = '%';
		break;

	    case 'd':
	    case 'i':
	    case 'u':
		/* Remember: char and short are sent as int on the stack */
		if (arg_size == -1)
		    l = (long)va_arg(ap, int);
		else if (arg_size == 1)
		    l = va_arg(ap, long); 
		else 
		    l = (long)va_arg(ap, int);

		out = append_int(out, (int)l);
		break;

	    case 'c':
		i = va_arg(ap, int);
		*out++ = i;
		break;

	    case 'f':
		d = va_arg(ap, double);
		out += sprintf(out, "%f", d);
		break;

	    case 'e':
	    case 'E':
	    case 'g':
	    case 'G':
		d = va_arg(ap, double);
		out += sprintf(out, "%g", d);
		break;

	    case 'p':
	    case 'x':
	    case 'X':
		l = (long)va_arg(ap, void *);
		puts("TODO");
		break;

	    case 'n':
		/* produces no output */
		break;

	    case 's': {
		char *s = (char *)va_arg(ap, char *);
		if (conv_len2) {
		    //memcpy(out, s, conv_len2);
		    //out += conv_len2;

		    //strncpy(out, s, conv_len2);
		    //out += MIN(conv_len2, strnlen(s));

		    while(conv_len2 && *s) {
		        *out++ = *s++;
		        conv_len2--;
		    }
		} else {
		    //strcpy(out, s);
		    //out += strlen(s);

		    while(*s) *out++ = *s++;
		}
		//		if (conv_len2) {
		//		    fwrite(s, conv_len2, 1, fp);
		//		} else
		//		    fputs(s, fp);
		break;
	    }

	    default:
		/* wchar_t types of 'C' and 'S' aren't supported */
		puts("TODO");
	    }
	    
	}

	case '\0':
	    break;

	default:
	    *out++ = *cp;
	}
    }

    *out = 0;
    fwrite(max_line, 1, out-max_line, fp);

    va_end(ap);
    
    return 0;
}
#endif

char *sam_aux_stringify(char *s, int len) {
    static char str[8192];
    char *cp = str, *s_end = s+len;
    int first = 1;

    //write(2, s, (int)(b->data + b->data_len - (uint8_t *)s));

    while (s < s_end) {
	if (first)
	    first = 0;
	else
	    *cp++ = '\t';

	switch (s[2]) {
	case 'A':
	    APPEND_FMT('A');
	    *cp++ = s[3];
	    s+=4;
	    break;

	case 'C':
	    APPEND_FMT('i');
	    cp = append_int(cp, *(uint8_t *)(s+3));
	    s+=4;
	    break;

	case 'c':
	    APPEND_FMT('i');
	    cp = append_int(cp, *(int8_t *)(s+3));
	    s+=4;
	    break;

	case 'S':
	    {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		APPEND_FMT('i');
		cp = append_int(cp, *(uint16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 's':
	    {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		APPEND_FMT('i');
		cp = append_int(cp, *(int16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 'I':
	    {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		APPEND_FMT('i');
		cp = append_int(cp, *(uint32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'i':
	    {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		APPEND_FMT('i');
		cp = append_int(cp, *(int32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'f':
	    {
		float f;
		memcpy(&f, s+3, 4);
		APPEND_FMT('f');
		cp += sprintf(cp, "%f", f);
	    }
	    s+=7;
	    break;

	case 'd':
	    {
		double d;
		memcpy(&d, s+3, 8);
		APPEND_FMT('f');
		cp += sprintf(cp, "%f", d);
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    APPEND_FMT(s[2]);
	    s+=3;
	    while (*s)
		*cp++ = *s++;
	    s++;
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *cp = 0;

    return str;
}

#ifdef HAVE_SAMTOOLS
char *bam_aux_stringify_old(bam1_t *b, int no_RG) {
    static char str[8192];
    char *s = (char *)bam1_aux(b), *cp = str;
    int first = 1;
    int keep;

    no_RG = 1;

    //write(2, s, (int)(b->data + b->data_len - (uint8_t *)s));

    while ((uint8_t *)s < b->data + b->data_len) {

	keep = (no_RG && s[0] == 'R' && s[1] == 'G') ? 0 : 1;
	if (keep) {
	    if (first)
		first = 0;
	    else
		*cp++ = '\t';
	}

	switch (s[2]) {
	case 'A':
	    if (keep)
		cp += sprintf(cp, "%c%c:A:%c", s[0], s[1], *(s+3));
	    s+=4;
	    break;

	case 'C':
	    if (keep)
		cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint8_t *)(s+3));
	    s+=4;
	    break;

	case 'c':
	    if (keep)
		cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int8_t *)(s+3));
	    s+=4;
	    break;

	case 'S':
	    if (keep) {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 's':
	    if (keep) {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 'I':
	    if (keep) {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'i':
	    if (keep) {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'f':
	    if (keep) {
		float f;
		memcpy(&f, s+3, 4);
		cp += sprintf(cp, "%c%c:f:%f", s[0], s[1], f);
	    }
	    s+=7;
	    break;

	case 'd':
	    if (keep) {
		double d;
		memcpy(&d, s+3, 8);
		cp += sprintf(cp, "%c%c:d:%f", s[0], s[1], d);
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    if (keep)
		cp += sprintf(cp, "%c%c:%c:%s", s[0], s[1], s[2], s+3);
	    s+=3;
	    while (*s++);
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *cp = 0;

    return str;
}

char *bam_aux_stringify(bam1_t *b, int no_RG) {
    static char str[8192];
    char *s = (char *)bam1_aux(b), *cp = str;
    int first = 1;
    int keep;

    no_RG = 1;

    //write(2, s, (int)(b->data + b->data_len - (uint8_t *)s));

    while ((uint8_t *)s < b->data + b->data_len) {

	keep = (no_RG && s[0] == 'R' && s[1] == 'G') ? 0 : 1;
	if (keep) {
	    if (first)
		first = 0;
	    else
		*cp++ = '\t';
	}

	switch (s[2]) {
	case 'A':
	    if (keep) {
		APPEND_FMT('A');
		*cp++ = s[3];
	    }
	    s+=4;
	    break;

	case 'C':
	    if (keep) {
		APPEND_FMT('i');
		cp = append_int(cp, *(uint8_t *)(s+3));
	    }
	    s+=4;
	    break;

	case 'c':
	    if (keep) {
		APPEND_FMT('i');
		cp = append_int(cp, *(int8_t *)(s+3));
	    }
	    s+=4;
	    break;

	case 'S':
	    if (keep) {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		APPEND_FMT('i');
		cp = append_int(cp, *(uint16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 's':
	    if (keep) {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		APPEND_FMT('i');
		cp = append_int(cp, *(int16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 'I':
	    if (keep) {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		APPEND_FMT('i');
		cp = append_int(cp, *(uint32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'i':
	    if (keep) {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		APPEND_FMT('i');
		cp = append_int(cp, *(int32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'f':
	    if (keep) {
		float f;
		memcpy(&f, s+3, 4);
		APPEND_FMT('f');
		cp += sprintf(cp, "%f", f);
	    }
	    s+=7;
	    break;

	case 'd':
	    if (keep) {
		double d;
		memcpy(&d, s+3, 8);
		APPEND_FMT('f');
		cp += sprintf(cp, "%f", d);
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    if (keep) {
		APPEND_FMT(s[2]);
		s+=3;
		while (*s)
		    *cp++ = *s++;
		s -= 3;
	    }
	    s+=3;
	    if (!keep)
		while (*s++);
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *cp = 0;

    return str;
}

/*
 * Filters out specific types from the sam aux records.
 * Returns a new aux record, also in the compact binary form.
 * 'len' holds the size of the returned string.
 *
 * Value returned is statically allocated. Do not free.
 */
char *bam_aux_filter(bam1_t *b, char **types, int ntypes, int *len) {
    static char str[8192];
    char *s = (char *)bam1_aux(b), *cp = str;
    int keep, i;

    while ((uint8_t *)s < b->data + b->data_len) {
	keep = 1;
	for (i = 0; i < ntypes; i++) {
	    if (s[0] == types[i][0] &&
		s[1] == types[i][1]) {
		keep = 0;
		break;
	    }
	}

	if (keep) {
	    *cp++ = s[0];
	    *cp++ = s[1];
	    *cp++ = s[2];
	}

	switch (s[2]) {
	case 'A':
	case 'C':
	case 'c':
	    if (keep)
		*cp++ = s[3];
	    s+=4;
	    break;

	case 'S':
	case 's':
	    if (keep) {
		*cp++ = s[3];
		*cp++ = s[4];
	    }
	    s+=5;
	    break;

	case 'I':
	case 'i':
	case 'f':
	    if (keep) {
		*cp++ = s[3];
		*cp++ = s[4];
		*cp++ = s[5];
		*cp++ = s[6];
	    }
	    s+=7;
	    break;

	case 'd':
	    if (keep) {
		*cp++ = s[3];
		*cp++ = s[4];
		*cp++ = s[5];
		*cp++ = s[6];
		*cp++ = s[7];
		*cp++ = s[8];
		*cp++ = s[9];
		*cp++ = s[10];
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    s+=3;
	    if (keep)
		while ((*cp++ = *s++));
	    else
		while (*s++);
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *len = cp - str;

    return str;
}

/*
 * Creates a new contig and updates the bam_io_t struct.
 */
void bio_new_contig(bam_io_t *bio, int tid) {
    char *cname = bio->header->target_name[tid];

    /* header->target_name[b.core.tid] */
    printf("\n++Processing contig %d / %s\n", tid, cname);
	
    create_new_contig(bio->io, &(bio->c), cname, bio->a->merge_contigs);
    bio->n_inserts = 0;
    bio->npads = 0;
    bio->skip = 0;

    if (bio->a->repad) {
	bio->tree = depad_consensus(bio->io, bio->c->rec);
	//padtree_dump(bio->tree);
    }
	
    bio->last_tid = tid;
}

/*
 * Samtools pileup won't iterate over unmapped reads. Therefore we have a
 * separate function to add these to the database - this one.
 * Although it shares much of the same code so is a candidate for merging
 * at some stage.
 */
int bio_add_unmapped(bam_io_t *bio, bam1_t *b) {
    char type;
    const char *LB;
    HacheItem *hi;
    HacheData hd;
    seq_t s;
    char tname[1024];
    library_t *lib = NULL;
    bam_aux_t val;
    int new = 0;
    char *name;
    int name_len;
    char *aux;
    int i, flags;
    tg_rec recno;
    int paired, is_pair = 0;
    char *filter[] = {"RG"};
    int stech;

    bio->count++;

    /* Check if it's a new contig, create if so */
    if (b->core.tid != bio->last_tid) {
	bio_new_contig(bio, b->core.tid);
    }

    /* Fetch read-group and pretend it's a library for now */
    if (0 == bam_aux_find(b, "RG", &type, &val) && type == 'Z') {
	LB = val.s;
	stech = stech_str2int(sam_tbl_get(bio->rg2pl_hash, LB));
    } else {
	LB = bio->fn;
	stech = STECH_UNKNOWN;
    }

    hd.p = NULL;
    hi = HacheTableAdd(bio->libs, (char *)LB, strlen(LB), hd, &new);
    if (new) {
	tg_rec lrec;
	printf("New library %s\n", LB);

	lrec = library_new(bio->io, (char *)LB);
	lib = get_lib(bio->io, lrec);
	lib = cache_rw(bio->io, lib);
	lib->machine = stech;
	hi->data.p = lib;
	cache_incr(bio->io, lib);
    }
    lib = hi->data.p;

    /* Construct a seq_t struct */
    name = bam1_qname(b);
    name_len = strlen(name);

    /*
     * Uncomment one of the two sections below if you wish to allow storing
     * sam auxillary records in gap5.
     * This is experimental and the format for these hasn't been finalised
     * yet.
     */
    aux = NULL;
    s.aux_len = 0;

    if (bio->a->sam_aux)
	aux = bam_aux_filter(b, filter, 1, &s.aux_len);

    //aux = bam_aux_stringify(b, 1);
    //s.aux_len = strlen(aux);
    
    //aux = bam1_aux(b);
    //s.aux_len = (int)(b->data + b->data_len - (uint8_t *)aux);

    s.pos = bio->npads +
	get_padded_coord(bio->tree, b->core.pos + 1 + bio->n_inserts
			 - bio->npads);
    //s.pos = b->core.pos+1;
    s.len = b->core.l_qseq;
    s.rec = 0;
    s.seq_tech = stech != STECH_UNKNOWN ? stech : stech_guess_by_name(name);
    s.flags = 0;
    s.left  = 1;
    s.right = s.len;
    s.parent_type = 0;
    s.parent_rec = 0;
    if (bio->a->data_type & DATA_NAME) {
	s.name_len = name_len;
	s.name = s.data = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, name);
    } else {
	char *n = "";
	s.name_len = 0;
	s.name = s.data = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, n);
    }
    s.trace_name = s.name + s.name_len + 1;
    *s.trace_name = 0;
    s.trace_name_len = 0;
    s.alignment = (unsigned char *)s.trace_name + s.trace_name_len + 1;
    *s.alignment = 0;
    s.alignment_len = 0;
    s.seq = (char *)s.alignment + s.alignment_len+1;
    s.conf = s.seq+s.len;
    s.mapping_qual = b->core.qual;
    s.format = SEQ_FORMAT_MAQ; /* pack bytes */
    s.anno = NULL;
    s.sam_aux = s.conf + s.len;

    for (i = 0; i < b->core.l_qseq; i++) {
	s.seq[i] = bio->a->data_type & DATA_SEQ
	    ? bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)]
	    : 'N';
	s.conf[i] = bio->a->data_type & DATA_QUAL
	    ? bam1_qual(b)[i]
	    : 0;
    }

    if (bam1_strand(b)) {
	complement_seq_t(&s);
	s.flags |= SEQ_COMPLEMENTED;
    }

    memcpy(s.sam_aux, aux, s.aux_len);

    /* Create the range, save the sequence */
    paired = (b->core.flag & BAM_FPAIRED) ? 1 : 0;
    flags = paired ? GRANGE_FLAG_TYPE_PAIRED : GRANGE_FLAG_TYPE_SINGLE;

    if (b->core.flag & BAM_FREAD1)
	s.flags |= SEQ_END_FWD;

    if (b->core.flag & BAM_FREAD2)
	s.flags |= SEQ_END_REV;

    strcpy(tname, name);

    if (name_len >= 2 && name[name_len-2] == '/') {
	tname[name_len-2] = 0;

	/* Check validity of name vs bit-fields */
	if ((name[name_len-1] == '1' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_FWD) ||
	    (name[name_len-1] == '2' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_REV)) {
	    fprintf(stderr, "Inconsistent read name vs flags: %s vs 0x%02x\n",
		    name, b->core.flag);
	}
    }

    if (paired)
	flags |= (s.flags & SEQ_END_MASK) == SEQ_END_FWD
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;
    else
	/* Guess work here. For now all <--- are rev, all ---> are fwd */
	flags |= bam1_strand(b)
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;

    if (bam1_strand(b)) {
	flags |= GRANGE_FLAG_COMP1;
    }

    if (bio->pair) is_pair = 1;

    flags |= GRANGE_FLAG_ISUMSEQ;

    recno = save_range_sequence(bio->io, &s, 0, 0, s.mapping_qual, bio->pair,
    				is_pair, tname, bio->c, bio->a, flags, lib);


#ifdef SAM_AUX_AS_TAG
    /* Make an annotation out of the sam auxillary data */
    aux = bam_aux_stringify(b, 1);
    if (aux && *aux) {
	anno_ele_t *e;
	bin_index_t *bin;
	range_t r;

	r.mqual = str2type("SAMX");
	r.start = s.pos;
	r.end = s.pos;
	r.pair_rec = recno;
	r.flags = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;
	r.rec = anno_ele_new(bio->io, 0, GT_Seq, recno, 0, r.mqual, aux);
	e = (anno_ele_t *)cache_search(bio->io, GT_AnnoEle, r.rec);
	e = cache_rw(bio->io, e);
	
	bin = bin_add_range(bio->io, &bio->c, &r, NULL, NULL);
	e->bin = bin->rec;
    }
#endif

    return 0;
}

/*
 * Removes a sequence from the bam_io_t struct.
 * This actually performs the main work of adding a sequence to the gap5
 * database.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bio_del_seq(bam_io_t *bio, const bam_pileup1_t *p, int snum) {
    bam_seq_t *bs;
    bam1_t *b;
    seq_t s;
    HacheItem *hi;
    tg_rec recno;
    int i, paired;
    int is_pair = 0;
    int flags;
    char tname[1024];
    library_t *lib = NULL;
    bam_aux_t val;
    char type;
    const char *LB;
    HacheData hd;
    int new = 0;
    char *name;
    int name_len;
    char *aux;
    char *filter[] = {"RG"};
    char *handle, aux_key[2];
    int stech;

    if (snum < 0 || snum >= bio->nseq)
	return -1;

    bio->count++;

    bs = &bio->seqs[snum];
    b = bs->b;

    /* Fetch read-group and pretend it's a library for now */
    if (0 == bam_aux_find(b, "RG", &type, &val) && type == 'Z') {
	LB = val.s;
	stech = stech_str2int(sam_tbl_get(bio->rg2pl_hash, LB));
    } else {
	LB = bio->fn;
	stech = STECH_UNKNOWN;
    }

    hd.p = NULL;
    hi = HacheTableAdd(bio->libs, (char *)LB, strlen(LB), hd, &new);
    if (new) {
	tg_rec lrec;
	printf("New library %s\n", LB);

	lrec = library_new(bio->io, (char *)LB);
	lib = get_lib(bio->io, lrec);
	lib = cache_rw(bio->io, lib);
	lib->machine = stech;
	hi->data.p = lib;
	cache_incr(bio->io, lib);
    }
    lib = hi->data.p;

    /*
    printf("\nSeq %d @ %6d: '%.*s' '%.*s' => nseq=%d\n",
	   snum, bs->pos, bs->seq_len, bs->seq, bs->seq_len, bs->conf,
	   bio->nseq-1);
    */

    /* Construct a seq_t struct */
    s.right = bs->seq_len;
    for (i = p->qpos+1; i < b->core.l_qseq; i++) {
	int base = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), i)];
	int qual = bam1_qual(p->b)[i];
	bio_extend_seq(bio, snum, base, qual);
    }
    
    name = bam1_qname(b);
    name_len = strlen(name);

    /*
     * Uncomment one of the two sections below if you wish to allow storing
     * sam auxillary records in gap5.
     * This is experimental and the format for these hasn't been finalised
     * yet.
     */
    aux = NULL;
    s.aux_len = 0;

    if (bio->a->sam_aux)
	aux = bam_aux_filter(b, filter, 1, &s.aux_len);
    
    //aux = bam_aux_stringify(b, 1);
    //s.aux_len = strlen(aux);

    //aux = bam1_aux(b);
    //s.aux_len = (int)(b->data + b->data_len - (uint8_t *)aux);
    
    s.rec = bs->rec;
    s.pos = bs->pos;
    s.len = bs->seq_len;
    s.seq_tech = stech != STECH_UNKNOWN ? stech : stech_guess_by_name(name);
    s.flags = 0;
    s.left = bs->left+1;
    s.parent_type = 0;
    s.parent_rec = 0;
    if (bio->a->data_type & DATA_NAME) {
	s.name_len = name_len;
	s.name = s.data = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, name);
    } else {
	char *n = "";
	s.name_len = 0;
	s.name = s.data = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, n);
    }
    s.trace_name = s.name + s.name_len + 1;
    *s.trace_name = 0;
    s.trace_name_len = 0;
    s.alignment = (unsigned char *)s.trace_name + s.trace_name_len + 1;
    *s.alignment = 0;
    s.alignment_len = 0;
    s.seq = (char *)s.alignment + s.alignment_len+1;
    s.conf = s.seq+s.len;
    s.mapping_qual = b->core.qual;
    s.format = SEQ_FORMAT_MAQ; /* pack bytes */
    s.anno = NULL;
    s.sam_aux = s.conf + s.len;

    if (bio->a->data_type & DATA_SEQ)
	memcpy(s.seq,  bs->seq,  s.len);
    else
	memset(s.seq, 'N', s.len);

    if (bio->a->data_type & DATA_QUAL)
	memcpy(s.conf, bs->conf, s.len);
    else
	memset(s.conf, 0, s.len);
    
    if (bam1_strand(b)) {
	complement_seq_t(&s);
	s.flags |= SEQ_COMPLEMENTED;
    }

    memcpy(s.sam_aux, aux, s.aux_len);

    /* Create the range, save the sequence */
    paired = (b->core.flag & BAM_FPAIRED) ? 1 : 0;
    flags = paired ? GRANGE_FLAG_TYPE_PAIRED : GRANGE_FLAG_TYPE_SINGLE;

    if (b->core.flag & BAM_FREAD1)
	s.flags |= SEQ_END_FWD;

    if (b->core.flag & BAM_FREAD2)
	s.flags |= SEQ_END_REV;

    strcpy(tname, name);

    if (name_len >= 2 && name[name_len-2] == '/') {
	tname[name_len-2] = 0;

	/* Check validity of name vs bit-fields */
	if ((name[name_len-1] == '1' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_FWD) ||
	    (name[name_len-1] == '2' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_REV)) {
	    fprintf(stderr, "Inconsistent read name vs flags: %s vs 0x%02x\n",
		    name, b->core.flag);
	}
    }

    if (paired)
	flags |= (s.flags & SEQ_END_MASK) == SEQ_END_FWD
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;
    else
	/* Guess work here. For now all <--- are rev, all ---> are fwd */
	flags |= bam1_strand(b)
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;

    if (bam1_strand(b)) {
	flags |= GRANGE_FLAG_COMP1;
    }

    if (bio->pair) is_pair = 1;

    recno = save_range_sequence(bio->io, &s, 0, 0, s.mapping_qual, bio->pair,
    				is_pair, tname, bio->c, bio->a, flags, lib);


#ifdef SAM_AUX_AS_TAG
    /* Make an annotation out of the sam auxillary data */
    aux = bam_aux_stringify(b, 1);
    if (aux && *aux) {
	anno_ele_t *e;
	bin_index_t *bin;
	range_t r;

	r.mqual = str2type("SAMX");
	r.start = s.pos;
	r.end = s.pos;
	r.pair_rec = recno;
	r.flags = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;
	r.rec = anno_ele_new(bio->io, 0, GT_Seq, recno, 0, r.mqual, aux);
	e = (anno_ele_t *)cache_search(bio->io, GT_AnnoEle, r.rec);
	e = cache_rw(bio->io, e);
	
	bin = bin_add_range(bio->io, &bio->c, &r, NULL, NULL);
	e->bin = bin->rec;
    }
#endif

    /* Add tags */
    handle = NULL;
    while (0 == bam_aux_iter(b, &handle, aux_key, &type, &val)) {
	range_t r;
	anno_ele_t *e;
	bin_index_t *bin;
	char *tokens[4], *cp, *tag_text, tag_type[5];
	int ntok;
	int tag_pos, tag_len;

	if (!(aux_key[0] == 'Z' && (aux_key[1] == 's' || aux_key[1] == 'c')))
	    continue;

	tokens[0] = val.s;
	for (ntok = 1, cp = val.s; *cp && ntok < 4; cp++) {
	    if (*cp == '|') {
		*cp = 0;
		tokens[ntok++] = cp+1;
	    }
	}

	/* Parse it */
	tag_type[0] = tag_type[1] = tag_type[2] = tag_type[3] = '-';
	tag_type[4] = 0;
	strncpy(tag_type, tokens[0], 4);
	tag_pos  = ntok >= 2 ? atoi(tokens[1]) : 0;
	tag_len  = ntok >= 3 ? atoi(tokens[2]) : 0;
	tag_text = ntok >= 4 ? (unescape_line(tokens[3]),tokens[3]) : NULL;

	/* Create the tag */
	r.mqual    = str2type(tag_type);
	r.start    = tag_pos-1 + s.pos;
	r.end      = tag_pos-1 + s.pos + tag_len-1;
	r.pair_rec = (aux_key[1] == 'c') ? bio->c->rec : recno;
	r.flags    = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;
	r.rec      = anno_ele_new(bio->io, 0, GT_Seq, recno, 0, r.mqual,
				  tag_text);

	/* Link it to a bin */
	e = (anno_ele_t *)cache_search(bio->io, GT_AnnoEle, r.rec);
	e = cache_rw(bio->io, e);

	bin = bin_add_range(bio->io, &bio->c, &r, NULL, NULL);
	e->bin = bin->rec;
    }

    /* Tidy up */
    if (bs->seq)  free(bs->seq);
    if (bs->conf) free(bs->conf);
    
    return 0;
}

/*
 * Adds a base to a sequence in the bam_io_t struct.
 * Returns 0 on success
 *        -1 on failure
 */
int bio_extend_seq(bam_io_t *bio, int snum, char base, int conf) {
    bam_seq_t *s;

    if (snum < 0 || snum >= bio->nseq)
	return -1;

    s = &bio->seqs[snum];

    /* Extend if appropriate */
    if (s->seq_len >= s->alloc_len) {
	s->alloc_len = (s->alloc_len + 100)*1.5;
	if (NULL == (s->seq  = (char *)realloc(s->seq,  s->alloc_len)))
	    return -1;
	if (NULL == (s->conf = (char *)realloc(s->conf, s->alloc_len)))
	    return -1;
    }
    
    /* And append */
    s->seq [s->seq_len] = base;
    s->conf[s->seq_len] = conf;
    s->seq_len++;

    return 0;
}

/*
 * Called once for every position (pos) in every contig (tid), with an 
 * array of n sequences (pl).
 *
 * We produce our own stack of reads that we append to, spotting new
 * sequences (p->is_head) and terminating sequences (p->is_tail) flushing
 * out our data to disk only on sequence removal. It's possible we could
 * hijack the pl->b bam1_t struct with our own pointer, putting it back to
 * the the bam1_t only on removal.
 */
int bio_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl,
		void *data) {
    bam_io_t *bio = (bam_io_t *)data;
    int i, j, k, insertions = 0;
    int np;

    //printf("\nCallback at pos=%d n=%d tid=%d\n", pos, n, tid);

    /*
     * tid is reference id - aka contig
     * n is number of sequences at this point (pos).
     * pl is the pileup info for these sequences.
     */
    /* Create new contig if appropriate */
    if (tid != bio->last_tid) {
	bio_new_contig(bio, tid);
    }
    
    np = 0;
    if (bio->a->repad) {
	if ((np=padtree_pad_at(bio->tree, pos+1+bio->n_inserts-bio->npads))) {
	    /* Add pads to match existing consensus gaps */
	    for (j = 0; j < np; j++) {
		if (bio->skip) {
		    //printf("Skipping import of pads\n");
		    bio->skip--;
		    continue;
		}

		//printf("Import pads from existing consensus at %d\n",
		//       pos+1+bio->n_inserts-bio->npads);
		for (i = 0; i < n; i++) {
		    bio_extend_seq(bio, i, '*', 0);
		}
	    }
	}
    }

    for (j = 0; j <= insertions; j++) {
	//printf("%d pos: %6d.%d ", tid, pos+1+bio->n_inserts, j);

//	printf("%5d: ", pos+1+bio->n_inserts);
//	for (i = 0; i < n; i++) {
//	    const bam_pileup1_t *p = &pl[i];
//	    printf("%d%d ", p->indel, j);
//	}
//	printf("\n");

	if (j && bio->a->repad /*&& j > np*/) {
	    /* FIXME: And not already in the originally padded version */
	    np=padtree_pad_at(bio->tree, pos+2);

	    //printf("Insert at %d: j=%d np=%d\n", pos+2+bio->n_inserts, j, np);
	    //	    contig_insert_base(bio->io, &bio->c, pos+2+bio->n_inserts,
	    //			       '*', 0);
	    if (j > np) {
		contig_insert_base(bio->io, &bio->c,
				   get_padded_coord(bio->tree, pos+2)
				   +bio->n_inserts,
				   '*', 0);
		bio->npads++;
	    } else {
		bio->n_inserts--;
		bio->skip++;
	    }
	}

	for (i = 0; i < n; i++) {
	    const bam_pileup1_t *p = &pl[i];
	    
	    if (j == 0 && p->is_head) {
		int i2, ppos;
		/* New sequence */
		//printf("^%c", p->b->core.qual+33);
		ppos = bio->npads
		     + get_padded_coord(bio->tree,
					pos + 1 + bio->n_inserts - bio->npads);
		i2 = bio_new_seq(bio, p, ppos);

		/*
		 * The following fails if we have a previous sequence that
		 * ended on a deletion. Eg CIGAR string 36M1D. This causes
		 * samtools tview to break too.
		 *
		 * 2/7/2010: It also fails when dealing with split reads
		 * using [0-9]+N CIGAR entries. We omit this assertion for
		 * now until we can implement spliced alignments better.
		 */
		//assert(i == i2);
	    }

	    if (p->is_del) {
		/* Undercall, aka deleted base compared to ref */
		int q = (bam1_qual(p->b)[p->qpos] + 
			 (p->qpos < p->b->core.l_qseq-1
			  ? bam1_qual(p->b)[p->qpos+1]
			  : bam1_qual(p->b)[p->qpos])) / 2;
		//printf("Undercall in seq %d\n", i);
		bio_extend_seq(bio, i, '*', q);
	    } else {
		if (p->indel >= 0) {
		    int c, q;

		    /* Overcall - edit all other seqs */
		    if (insertions < p->indel)
			insertions = p->indel;

		    if (j <= p->indel) {
			c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b),
							 p->qpos + j)];
			q = bam1_qual(p->b)[p->qpos+j];
		    } else {
			c = '*';
			q = (bam1_qual(p->b)[p->qpos] + 
			     (p->qpos < p->b->core.l_qseq-1
			      ? bam1_qual(p->b)[p->qpos+1]
			      : bam1_qual(p->b)[p->qpos])) / 2;
			if (p->is_tail)
			    continue;
		    }
		    bio_extend_seq(bio, i, c, q);
		} else if (p->indel < 0 /*&& j == 0*/) {
		    /*
		     * This occurs when the next base is an undercall.
		     * However THIS base should still exist, I think.
		     * It's all a little bit confusing if truth be known.
		     */
		    int c = j == 0
		      ? bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)]
		      : '*';
		    int q = bam1_qual(p->b)[p->qpos];
		    //printf("%s extend seq %d base %c\n",
		    //       bam1_qname(p->b), i, c);
		    bio_extend_seq(bio, i, c, q);
		}
	    }
	}
    }

    for (i = j = k = 0; i < n; i++, j++, k++) {
	const bam_pileup1_t *p = &pl[j];
	
	bio->seqs[i] = bio->seqs[j];
	if (p->is_tail) {
	    bio_del_seq(bio, p, i);
	    i--;
	    n--;
	    bio->nseq--;
	}
    }

    bio->n_inserts += insertions;
    //    if (insertions) {
    //	printf("%d insertions (%d) at pos=%d\n", insertions, bio->n_inserts, pos);
    //    }

    return 0;
}

int parse_sam_or_bam(GapIO *io, const char *fn, tg_args *a, char *mode) {
    bam_io_t *bio = (bam_io_t*)calloc(1, sizeof(*bio));
    bam1_t *b;
    int count = 0;
    bam_plbuf_t *plbuf;
    samfile_t *fp;
    rec_list_t *tmp;

    /* for pair data */
    open_tmp_file();

    /* Setup bam_io_t object and create our pileup interface */
    bio->io = io;
    bio->seqs = NULL;
    bio->nseq = 0;
    bio->max_seq = 0;
    bio->a = a;
    bio->c = NULL;
    bio->count = 0;
    bio->fn = fn;
    bio->libs = HacheTableCreate(256, HASH_DYNAMIC_SIZE);
    bio->libs->name = "libs";
    bio->last_tid = -1;
    bio->tree = NULL;
    bio->rec_head = bio->rec_tail = NULL;

    if (a->pair_reads) {
	bio->pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
	bio->pair->name = "pair";
    } else {
	bio->pair = NULL;
    }

    fp = samopen(fn, mode, NULL);
    assert(fp);
    bio->header = fp->header;
    if (!bio->header)
	return -1;

    if (!bio->header->dict) {
	bio->header->dict = sam_header_parse2(bio->header->text);
    }
    bio->rg2pl_hash = sam_header2tbl(bio->header->dict, "RG", "ID", "PL");
    plbuf = bam_plbuf_init(bio_callback, bio);
    //bam_plbuf_set_mask(plbuf, BAM_DEF_MASK /* or BAM_FUNMAP? */);
    bam_plbuf_set_mask(plbuf, 0 /* or BAM_DEF_MASK, or BAM_FUNMAP? */);

    /*
     * Loop through reads - the bulk of the work here is done in the
     * bio_callback function.
     */
    b = (bam1_t*)calloc(1, sizeof(bam1_t));
    while (samread(fp, b) >= 0) {
	int k;

	if (a->store_unmapped && b->core.flag & BAM_FUNMAP) {
	    bio_add_unmapped(bio, b);
	    if ((++count & 0xffff) == 0) {
	        putchar('.');
	        fflush(stdout);
	        cache_flush(io);
	    }
	    continue;
	} else if (b->core.flag & BAM_FUNMAP)
	    continue;

	//printf("push %s\n", bam1_qname(b));
	/*
	 * Allocate a record number now, so they're used in the same order
	 * as the input data regardless of whether the read is mapped or
	 * unmapped.
	 *
	 * Note, for spliced alignments we may need multiple record numbers.
	 */
	if (NULL == (tmp = xmalloc(sizeof(*tmp))))
	    return -1;

	tmp->next = NULL;
	if (bio->rec_head) {
	    bio->rec_tail->next = tmp;
	    bio->rec_tail = tmp;
	} else {
	    bio->rec_head = bio->rec_tail = tmp;
	}
	if (bio->a->data_type == DATA_BLANK) {
	    static int fake_recno = 1;
	    tmp->rec = fake_recno++;
	} else {
	    tmp->rec = sequence_new_from(bio->io, NULL);
	}

	for (k = 0; k < b->core.n_cigar; k++) {
	    int op = bam1_cigar(b)[k] & BAM_CIGAR_MASK;
	    if (op == BAM_CREF_SKIP) {
		if (NULL == (tmp = xmalloc(sizeof(*tmp))))
		    return -1;
		tmp->next = NULL;
		bio->rec_tail->next = tmp;
		bio->rec_tail = tmp;

		if (bio->a->data_type == DATA_BLANK) {
		    static int fake_recno = 1;
		    tmp->rec = fake_recno++;
		} else {
		    tmp->rec = sequence_new_from(bio->io, NULL);
		}
	    }
	}

	bam_plbuf_push(b, plbuf);
	if ((++count & 0xffff) == 0) {
	    putchar('.');
	    fflush(stdout);
	    cache_flush(io);
	}
    }
    putchar('\n');
    bam_plbuf_push(0, plbuf);
    bam_plbuf_destroy(plbuf);
    if (bio->rg2pl_hash)
	sam_tbl_destroy(bio->rg2pl_hash);

    cache_flush(io);
    vmessage("Loaded %d of %d sequences\n", bio->count, count);

    if (bio->pair && !a->fast_mode) {    
	sort_pair_file();
	
	complete_pairs(io);
	
	close_tmp_file();
    }
 
    /* Tidy up */
    if (b) {
	if (b->data)
	    free(b->data);
	free(b);
    }

    if (fp)
	samclose(fp);

    if (bio) {
	if (bio->pair)
	    HacheTableDestroy(bio->pair, 1);

	if (bio->libs) {
	    /* call cache_decr on each lib too */
	    HacheItem *hi;
	    HacheIter *iter;

	    if (!(iter = HacheTableIterCreate()))
		return -1;
	
	    while (hi = HacheTableIterNext(bio->libs, iter)) {
		library_t *lib = hi->data.p;
		cache_decr(io, lib);
	    }

	    HacheTableDestroy(bio->libs, 0);
	}

	if (bio->seqs)
	    free(bio->seqs);

	if (bio->tree)
	    depad_seq_tree_free(bio->tree);

	while (bio->rec_head) {
	    tmp = bio->rec_head->next;
	    free(bio->rec_head);
	    bio->rec_head = tmp;
	}

	free(bio);

	if (bio->c)
	    cache_decr(io, bio->c);
    }

    return 0;
}

int parse_bam(GapIO *io, const char *fn, tg_args *a) {
    return parse_sam_or_bam(io, fn, a, "rb");
}

int parse_sam(GapIO *io, const char *fn, tg_args *a) {
    return parse_sam_or_bam(io, fn, a, "r");
}

#endif /* HAVE_SAMTOOLS */

