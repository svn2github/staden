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

typedef struct {
    GapIO *io;
    const char *fn;
    //bam_seq_t *seqs;
    int nseq;
    int max_seq;
    //rec_list_t *rec_head;
    //rec_list_t *rec_tail;
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

typedef union {
    char  *s;
    int    i;
    float  f;
    double d;
} bam_aux_t;

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
static int bam_aux_find(bam1_t *b, char *key, char *type, bam_aux_t *val) {
    char *h = NULL;
    char k[2];

    while (0 == bam_aux_iter(b, &h, k, type, val)) {
	if (k[0] == key[0] && k[1] == key[1])
	    return 0;
    }

    return -1;
}

/*
 * Filters out specific types from the sam aux records.
 * Returns a new aux record, also in the compact binary form.
 * 'len' holds the size of the returned string.
 *
 * Value returned is statically allocated. Do not free.
 */
static char *bam_aux_filter(bam1_t *b, char **types, int ntypes, int *len) {
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
 * Attempts to guess a sequencing technology based on sequence name.
 * This is far from ideal and is also a bit Sanger Institute centric.
 */
static int stech_guess_by_name(char *str) {
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

/*
 * Creates a new contig and updates the bam_io_t struct.
 */
static void bio_new_contig(bam_io_t *bio, int tid) {
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

static int bio_add_seq(bam_io_t *bio, bam1_t *b) {
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
    int start;
    int end;
    static unsigned char *alignment = NULL;
    static size_t alignment_alloc = 0;
    size_t alignment_len = 0;
    char *handle, aux_key[2];

    bio->count++;

    //    printf("Processing seq %d\n", bio->count);

    /* Create a new contig if needed */
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


    /* Parse CIGAR string to get clips and other bits */
    end = start = b->core.pos + 1;
    //printf("%s\t", bam1_qname(b));
    s.left  = 1;
    s.right = b->core.l_qseq;
    for (i = 0; i < b->core.n_cigar; i++) {
	int cig = bam1_cigar(b)[i];
	int op  = cig &  BAM_CIGAR_MASK;
	int len = cig >> BAM_CIGAR_SHIFT;

	char c;
	switch (op) {
	case BAM_CMATCH:     c = 'M'; end += len; break;
	case BAM_CINS:       c = 'I'; break;
	case BAM_CDEL:       c = 'D'; end += len; break;
	case BAM_CREF_SKIP:  c = 'N'; end += len; break;
	case BAM_CPAD:       c = 'P'; puts("CIGAR P not supported yet"); break;
	case BAM_CSOFT_CLIP:
	    c = 'S';
	    if (i == 0) {
		start -= len;
		s.left += len;
	    } else if (i == b->core.n_cigar - 1) {
		s.right -= len;
		end +=len;
	    }
	    break;
	case BAM_CHARD_CLIP: c = 'H'; break;
	}

	//printf(" %d%c", len, c);

	/* Store 'S' as 'M' - we hold soft-clips in another form */
	if (c == 'S') c = 'M';
	if (alignment_len > 0 && alignment[alignment_len-1] == c) {
	    len += alignment[alignment_len-2];
	    alignment_len -= 2;
	}

	do {
	    if (alignment_len+2 >= alignment_alloc) {
		alignment_alloc = alignment_alloc ? alignment_alloc*2 : 8192;
		alignment = realloc(alignment, alignment_alloc);
		if (!alignment)
		    return -1;
	    }
	    alignment[alignment_len++] = len <= 255 ? len : 255;
	    alignment[alignment_len++] = c;
	    len -= 255;
	} while (len > 0);
    }
    end--;
    //printf(" => %d .. %d, %d vs %d\n", start, end, end-start+1, b->core.l_qseq);

    /* Allocate variable sized data and fill them out */
    /* Construct the fixed portions of a seq_t struct */
    name = bam1_qname(b);
    name_len = strlen(name);

    aux = NULL;
    s.aux_len = 0;

    if (bio->a->sam_aux)
	aux = bam_aux_filter(b, filter, 1, &s.aux_len);

    s.pos = b->core.pos + 1;
    s.len = b->core.l_qseq;
    s.rec = 0;
    s.seq_tech = stech != STECH_UNKNOWN ? stech : stech_guess_by_name(name);
    s.flags = 0;
    s.parent_type = 0;
    s.parent_rec = 0;

    if (bio->a->data_type & DATA_NAME) {
	s.name_len = name_len;
	s.name = s.data = (char *)malloc(s.name_len + 3 + 2*s.len +
					 alignment_len + 1 + s.aux_len);
	strcpy(s.name, name);
    } else {
	char *n = "";
	s.name_len = 0;
	s.name = s.data = (char *)malloc(s.name_len + 3 + 2*s.len +
					 alignment_len + 1 + s.aux_len);
	strcpy(s.name, n);
    }
    s.trace_name = s.name + s.name_len + 1;
    *s.trace_name = 0;
    s.trace_name_len = 0;
    s.alignment = (unsigned char *)s.trace_name + s.trace_name_len + 1;
    memcpy(s.alignment, alignment, alignment_len);
    s.alignment[alignment_len] = 0;
    s.alignment_len = alignment_len;
    s.seq = (char *)s.alignment + s.alignment_len+1;
    s.conf = s.seq+s.len;
    s.mapping_qual = b->core.qual;
    s.format = SEQ_FORMAT_MAQ; /* pack bytes */
    s.anno = NULL;
    s.sam_aux = s.conf + s.len;
    memcpy(s.sam_aux, aux, s.aux_len);
    
    for (i = 0; i < b->core.l_qseq; i++) {
	s.seq[i] = bio->a->data_type & DATA_SEQ
	    ? bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)]
	    : 'N';
	s.conf[i] = bio->a->data_type & DATA_QUAL
	    ? bam1_qual(b)[i]
	    : 0;
    }

    /* Complement as necessary */
    if (bam1_strand(b)) {
	complement_seq_t(&s);
	s.flags |= SEQ_COMPLEMENTED;
    }


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

    recno = save_range_sequence(bio->io, &s, start, end,
				s.mapping_qual, bio->pair,
    				is_pair, tname, bio->c, bio->a, flags, lib);

    /* Add tags */
    handle = NULL;
    while (0 == bam_aux_iter(b, &handle, aux_key, &type, &val)) {
	range_t r;
	anno_ele_t *e;
	bin_index_t *bin;
	char *tokens[4], *cp, *tag_text, tag_type[5];
	int ntok;
	int tag_st, tag_st_nth, tag_en, tag_en_nth, tag_len;
	char sep = 0;

	if (!(aux_key[0] == 'Z' && (aux_key[1] == 's' || aux_key[1] == 'c')))
	    continue;

	tokens[0] = val.s;
	for (ntok = 1, cp = val.s; *cp && ntok < 4; cp++) {
	    if ((sep && *cp == sep) || *cp == '|' || *cp == '/') {
		sep = *cp;
		*cp = 0;
		tokens[ntok++] = cp+1;
	    }
	}

	/* Parse it */
	tag_type[0] = tag_type[1] = tag_type[2] = tag_type[3] = '-';
	tag_type[4] = 0;
	strncpy(tag_type, tokens[0], 4);
	if (sep == '|') {
	    tag_st     = ntok >= 2 ? atoi(tokens[1]) : 0;
	    tag_len    = ntok >= 3 ? atoi(tokens[2]) : 0;
	    tag_en     = tag_st + tag_len-1;
	    tag_st_nth = 0;
	    tag_en_nth = 0;
	} else {
	    char *endp;
	    tag_st     = ntok >= 2 ? strtol(tokens[1], &endp, 10) : 0;
	    tag_st_nth = *endp == '+' ? atoi(endp+1) : 0;
	    tag_en     = ntok >= 3 ? strtol(tokens[2], &endp, 10) : 0;
	    tag_en_nth = *endp == '+' ? atoi(endp+1) : 0;
	}
	tag_text = ntok >= 4 ? (unescape_line(tokens[3]),tokens[3]) : NULL;

	/* Create the tag */
	r.mqual    = str2type(tag_type);
	r.start    = tag_st-1 + start;
	r.start_nth= tag_st_nth;
	r.end      = tag_en-1 + start;
	r.end_nth  = tag_en_nth;
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

    return 0;
}

static int parse_unpadded_sam_or_bam(GapIO *io, const char *fn,
				     tg_args *a, char *mode) {
    bam_io_t *bio = (bam_io_t*)calloc(1, sizeof(*bio));
    bam1_t *b;
    samfile_t *fp;

    /* for pair data */
    open_tmp_file();

    /* Setup bam_io_t object */
    bio->io = io;
    //bio->seqs = NULL;
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
    //bio->rec_head = bio->rec_tail = NULL;

    if (a->pair_reads) {
	bio->pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
	bio->pair->name = "pair";
    } else {
	bio->pair = NULL;
    }

    /* Open the file and deal with headers etc */
    fp = samopen(fn, mode, NULL);
    if (!fp)
	return -1;
    bio->header = fp->header;
    if (!bio->header)
	return -1;

    if (!bio->header->dict) {
	bio->header->dict = sam_header_parse2(bio->header->text);
    }
    bio->rg2pl_hash = sam_header2tbl(bio->header->dict, "RG", "ID", "PL");

    /*
     * Loop through reads in an unpadded form.
     */
    b = (bam1_t*)calloc(1, sizeof(bam1_t));
    while (samread(fp, b) >= 0) {
	if (!a->store_unmapped && b->core.flag & BAM_FUNMAP)
	    continue;

	bio_add_seq(bio, b);

	if ((bio->count & 0xffff) == 0) {
	    putchar('.');
	    fflush(stdout);
	    cache_flush(io);
	}
    }
    putchar('\n');
    if (bio->rg2pl_hash)
	sam_tbl_destroy(bio->rg2pl_hash);

    cache_flush(io);
    vmessage("Loaded %d sequences\n", bio->count);

    if (bio->pair && !a->fast_mode) {    
	sort_pair_file();
	
	complete_pairs(io);
	
	close_tmp_file();
    }
 
    return 0;
}

int parse_unpadded_bam(GapIO *io, const char *fn, tg_args *a) {
    return parse_unpadded_sam_or_bam(io, fn, a, "rb");
}

int parse_unpadded_sam(GapIO *io, const char *fn, tg_args *a) {
    return parse_unpadded_sam_or_bam(io, fn, a, "r");
}

#endif /* HAVE_SAMTOOLS */

