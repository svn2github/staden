#ifndef _ZTR_H
#define _ZTR_H

#include "Read.h"

/* The header */
typedef struct {
    unsigned char  magic[8];	  /* 0xae5a54520d0a1a0a (be) */
    unsigned char  version_major; /* ZTR_VERSION_MAJOR */
    unsigned char  version_minor; /* ZTR_VERSION_MINOR */
} ztr_header_t;

/* The ZTR magic numbers */
#define ZTR_MAGIC   		"\256ZTR\r\n\032\n"
#define ZTR_VERSION_MAJOR	1
#define ZTR_VERSION_MINOR	2

/*
 * CHUNKS
 *
 * Chunks consist of a block length followed by the type, format and data.
 */

typedef struct {
    uint4 type;			/* chunk type (be) */
    uint4 mdlength;		/* length of meta data field (be) */
    char *mdata;		/* meta data */
    uint4 dlength;		/* length of data field (be) */
    char *data;			/* a format byte and the data itself */
} ztr_chunk_t;

/* Format types */
#define ZTR_FORM_RAW		0
#define ZTR_FORM_RLE		1
#define ZTR_FORM_ZLIB		2
#define ZTR_FORM_XRLE		3
#define ZTR_FORM_DELTA1		64
#define ZTR_FORM_DELTA2		65
#define ZTR_FORM_DELTA4		66
#define ZTR_FORM_DDELTA1	67
#define ZTR_FORM_DDELTA2	68
#define ZTR_FORM_DDELTA4	69
#define ZTR_FORM_16TO8		70
#define ZTR_FORM_32TO8		71
#define ZTR_FORM_FOLLOW1	72
#define ZTR_FORM_CHEB445	73
#define ZTR_FORM_ICHEB		74

/* Converts a C string to a big-endian 4-byte int */
#define ZTR_STR2BE(str) (((str)[0] << 24) + \
                         ((str)[1] << 16) + \
                         ((str)[2] <<  8) + \
                         ((str)[3] <<  0))

/* Converts a big-endian 4-byte int to a C string */
#define ZTR_BE2STR(i,str) (((str)[0]=((i)>>24)&0xff),\
                           ((str)[1]=((i)>>16)&0xff),\
                           ((str)[2]=((i)>> 8)&0xff),\
                           ((str)[3]=((i)>> 0)&0xff),\
			   (str)[4]='\0',str)\

#define ZTR_TYPE_HEADER	0xae5a5452 /* M-. Z T R */

#define ZTR_TYPE_SAMP	0x53414d50
#define ZTR_TYPE_SMP4	0x534d5034
#define ZTR_TYPE_BASE	0x42415345
#define ZTR_TYPE_BPOS	0x42504f53
#define ZTR_TYPE_CNF4	0x434e4634
#define ZTR_TYPE_CNF1	0x434e4631
#define ZTR_TYPE_CSID	0x43534944
#define ZTR_TYPE_TEXT	0x54455854
#define ZTR_TYPE_CLIP	0x434c4950
#define ZTR_TYPE_COMM	0x434f4d4d
#define ZTR_TYPE_CR32	0x43523332
#define ZTR_TYPE_FLWO	0x464c574f
#define ZTR_TYPE_FLWC	0x464c5743

/* A text segment consists of identifier and value */
typedef struct {
    char *ident; /* Pointer to identifier */
    char *value; /* Pointer to value */
} ztr_text_t;

/* The main ZTR structure, which holds the entire file contents */
typedef struct {
    /* General bits to do with the ZTR file format */
    ztr_header_t header;	/* File Header */
    ztr_chunk_t *chunk;		/* Array of chunks */
    int nchunks;		/* Number of chunks */

    /* Specifics to do with the standard chunk types */
    ztr_text_t *text_segments;
    int ntext_segments;

    /* 'Hint' for delta of SAMP and SMP4 */
    int delta_level;
} ztr_t;

int fwrite_ztr(FILE *fp, ztr_t *ztr);
int mfwrite_ztr(mFILE *fp, ztr_t *ztr);
ztr_t *fread_ztr(FILE *fp);
ztr_t *mfread_ztr(mFILE *fp);
Read *ztr2read(ztr_t *ztr);
ztr_t *read2ztr(Read *r);
int compress_ztr(ztr_t *ztr, int level);
int uncompress_ztr(ztr_t *ztr);
ztr_t *new_ztr(void);
void delete_ztr(ztr_t *ztr);
ztr_chunk_t **ztr_find_chunks(ztr_t *ztr, uint4 type, int *nchunks_p);
void ztr_process_text(ztr_t *ztr);

#endif /* _ZTR_H */
