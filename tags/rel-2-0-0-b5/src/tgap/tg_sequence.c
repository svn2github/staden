#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "tg_gio.h"

#ifndef ABS
#    define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

/*
 * Given a seq_t struct this updates the internal pointers to be valid offsets
 * into the s->data field. This is useful if the structure has been copied to
 * a new address.
 */
void sequence_reset_ptr(seq_t *s) {
    if (!s) return;

    s->name = (char *)&s->data;
    s->trace_name = s->name + s->name_len + 1;
    s->alignment = s->trace_name + s->trace_name_len + 1;
    s->seq = s->alignment + s->alignment_len + 1;
    s->conf = s->seq + (s->len >= 0 ? s->len : -s->len);
    if (s->aux_len)
	s->sam_aux = s->conf + 
	    (s->format == SEQ_FORMAT_CNF4 ? 4 : 1) *
	    (s->len >= 0 ? s->len : -s->len);
    else
	s->sam_aux = NULL;
}

/*
 * Returns the size needed to store confidence values in this sequence.
 * Ie 1 or 4 per base.
 */
#define sequence_conf_size(s) ((s)->format == SEQ_FORMAT_CNF4 ? 4 : 1)

size_t sequence_extra_len(seq_t *s) {
    return
	(s->name       ? strlen(s->name)       : 0) + 1 +
	(s->trace_name ? strlen(s->trace_name) : 0) + 1 + 
	(s->alignment  ? strlen(s->alignment)  : 0) + 1 + 
	ABS(s->len)                                 + 
	ABS(s->len) * sequence_conf_size(s)         +
	s->aux_len;
}

/*
 * Copies the 'f' seq_t struct to the 's' seq_t struct.
 * Assumes 's' has already been allocated to be large enough to hold 'f'.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int  sequence_copy(seq_t *s, seq_t *f) {
    int rec, idx;
    seq_block_t *block;

    if (!s || !f)
	return -1;

    /* Copy almost all */
    rec   = s->rec;
    block = s->block;
    idx   = s->idx;
    *s = *f;
    s->rec   = rec;
    s->block = block;
    s->idx   = idx;

    /* Fix internal pointers */
    sequence_reset_ptr(s);

    /* Copy data */
    strcpy(s->name, f->name ? f->name : "");
    s->name_len = strlen(s->name);

    strcpy(s->trace_name, f->trace_name ? f->trace_name : "");
    s->trace_name_len = strlen(s->trace_name);

    strcpy(s->alignment, f->alignment ? f->alignment : "");
    s->alignment_len = strlen(s->alignment);

    memcpy(s->seq, f->seq, ABS(f->len));

    memcpy(s->conf, f->conf, ABS(f->len)*
	   (f->format == SEQ_FORMAT_CNF4 ? 4 : 1));
    
    if (s->aux_len)
	memcpy(s->sam_aux, f->sam_aux, s->aux_len);

    if (s->anno) {
	s->anno = ArrayCreate(sizeof(int), ArrayMax(f->anno));
	memcpy(ArrayBase(int, s->anno),
	       ArrayBase(int, f->anno),
	       ArrayMax(f->anno) * sizeof(int));
    }

    return 0;
}


/*
 * Given a seq_t struct this allocates a new sequence in the database
 * and copies the contents of 's' into it. If 's' is NULL it simply allocates
 * the new sequence and does nothing with it.
 *
 * Note if s->rec is non-zero it assumes that a record number has already
 * been allocated for this sequence.
 *
 * Returns the record number on success
 *        -1 on failure
 */
#if 0
int sequence_new_from(GapIO *io, seq_t *s) {
    return io->iface->seq.create(io->dbh, s);
}
#else
int sequence_new_from(GapIO *io, seq_t *s) {
    int rec;
    seq_t *n;

    if (s && s->rec) {
	cache_item_init(io, GT_Seq, s, s->rec);
	rec = s->rec;
    } else { 
	rec = cache_item_create(io, GT_Seq, s);
    }

    if (s) {
	n = (seq_t *)cache_search(io, GT_Seq, rec);
	n = cache_rw(io, n);
	n = cache_item_resize(n, sizeof(*n) + sequence_extra_len(s));

	if (sequence_copy(n, s) == -1)
	    return -1;
    }

    //printf("%d -> %.*s\n", rec, s->name_len, s->name);

    return rec;
}
#endif

/*
 * Sets the sequence position
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_position(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->pos = value;
    *s = n;

    return 0;
}

/*
 * Sets the sequence length
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_len(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->len = value;
    *s = n;

    return 0;
}

/*
 * Sets the sequence left clip
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_left(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->left = value;
    *s = n;

    sequence_invalidate_consensus(io, n);

    return 0;
}

/*
 * Sets the sequence right clip
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_right(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->right = value;
    *s = n;

    sequence_invalidate_consensus(io, n);

    return 0;
}

/*
 * Sets the sequence mapping quality
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_mapping_qual(GapIO *io, seq_t **s, uint8_t value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->mapping_qual = value;
    *s = n;

    return 0;
}


/*
 * Sets the index into the bin ranges array referring to this seq.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_bin_index(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->bin_index = value;
    *s = n;

    return 0;
}

/*
 * Sets the parent record type (GT_Template, GT_Ligation, etc)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_parent_type(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->parent_type = value;
    *s = n;

    return 0;
}

/*
 * Sets the parent record number.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_parent_rec(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->parent_rec = value;
    *s = n;

    return 0;
}

/*
 * Sets the flags. See SEQ_* defines in tg_struct.h
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_flags(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->flags = value;
    *s = n;

    return 0;
}

/*
 * Sets the sequence technology. See STECH_* defines in tg_struct.h
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_seq_tech(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->seq_tech = value;
    *s = n;

    return 0;
}

/*
 * Sets a sequence name.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_name(GapIO *io, seq_t **s, char *name) {
    size_t extra_len;
    seq_t *n;
    char *tmp,*cp;

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    extra_len = sequence_extra_len(*s);
    extra_len += (name       ? strlen(name)       : 0) -
	         ((*s)->name ? strlen((*s)->name) : 0);
    n = cache_item_resize(n, sizeof(*n) + extra_len);
    if (NULL == n)
	return -1;
    *s = n;

    n->name_len = strlen(name);
    sequence_reset_ptr(n);

    /* Shift and insert name */
    cp = tmp = malloc(extra_len);
    strcpy(cp, name);
    cp += n->name_len+1;
    strcpy(cp, n->trace_name);
    cp += n->trace_name_len;
    strcpy(cp, n->alignment);
    cp += n->alignment_len;
    memcpy(cp, n->seq, ABS(n->len));
    cp += ABS(n->len);
    memcpy(cp, n->conf, ABS(n->len) * sequence_conf_size(n));
    cp += ABS(n->len) * sequence_conf_size(n);
    if (n->aux_len)
	memcpy(cp, n->sam_aux, n->aux_len);
    memcpy(&n->data, tmp, extra_len);
    free(tmp);
    
    return 0;
}

/*
 * Sets a trace name.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_trace_name(GapIO *io, seq_t **s, char *trace_name) {
    size_t extra_len;
    seq_t *n;
    char *tmp, *cp;

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    if (!trace_name || 0 == strcmp(n->name, trace_name))
	trace_name = "";

    extra_len = sequence_extra_len(*s);
    extra_len += (trace_name       ? strlen(trace_name)       : 0) -
	         ((*s)->trace_name ? strlen((*s)->trace_name) : 0);

    n = cache_item_resize(n, extra_len);
    if (NULL == n)
	return -1;
    *s = n;

    n->trace_name_len = strlen(trace_name);
    sequence_reset_ptr(n);

    /* Shift and insert name */
    cp = tmp = malloc(extra_len);
    strcpy(cp, n->name);
    cp += n->name_len+1;
    strcpy(cp, trace_name);
    cp += n->trace_name_len;
    strcpy(cp, n->alignment);
    cp += n->alignment_len;
    memcpy(cp, n->seq, ABS(n->len));
    cp += ABS(n->len);
    memcpy(cp, n->conf, ABS(n->len) * sequence_conf_size(n));
    cp += ABS(n->len) * sequence_conf_size(n);
    if (n->aux_len)
	memcpy(cp, n->sam_aux, n->aux_len);
    memcpy(&n->data, tmp, extra_len);
    free(tmp);

    return 0;
}

int sequence_set_seq (GapIO *io, seq_t **s, char *seq) {return -1;}
int sequence_set_conf(GapIO *io, seq_t **s, char *conf) {return -1;}


/* ------------------------------------------------------------------------ 
 * Trivial one-off sequence query functions
 */
int seq_pos(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_pos(&s);
}

int seq_len(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_len(&s);
}

int seq_left(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_left(&s);
}

int seq_right(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_right(&s);
}

int seq_mapping_qual(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_mapping_qual(&s);
}

char *seq_name(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_name(&s);
}

char *seq_seq(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_seq(&s);
}

char *seq_conf(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_conf(&s);
}

/* ------------------------------------------------------------------------ */
/* Sequence manipulation - ripped out of Staden Package's seq_utils.c */

/*
 * Reverses and complements a piece of DNA
 */
static int complementary_base[256];
void complement_seq_conf(char *seq, char *conf, int seq_len, int nconf) {
    int i, middle, j;
    char temp, t[4];
    static int init = 0;

    if (!init) {
	for (i = 0; i < 256; i++)
	    complementary_base[i] = i;

	complementary_base['a'] = 't';
	complementary_base['c'] = 'g';
	complementary_base['g'] = 'c';
	complementary_base['t'] = 'a';
	complementary_base['u'] = 'a';
	complementary_base['A'] = 'T';
	complementary_base['C'] = 'G';
	complementary_base['G'] = 'C';
	complementary_base['T'] = 'A';
	complementary_base['U'] = 'A';

	complementary_base['n'] = 'n';
	complementary_base['-'] = '-';
	complementary_base['b'] = 'v';
	complementary_base['d'] = 'h';
	complementary_base['h'] = 'd';
	complementary_base['k'] = 'm';
	complementary_base['m'] = 'k';
	complementary_base['r'] = 'y';
	complementary_base['s'] = 's';
	complementary_base['v'] = 'b';
	complementary_base['w'] = 'w';
	complementary_base['y'] = 'r';

	complementary_base['B'] = 'V';
	complementary_base['D'] = 'H';
	complementary_base['H'] = 'D';
	complementary_base['K'] = 'M';
	complementary_base['M'] = 'K';
	complementary_base['R'] = 'Y';
	complementary_base['S'] = 'S';
	complementary_base['V'] = 'B';
	complementary_base['W'] = 'W';
	complementary_base['Y'] = 'R';
	init = 1;
    }

    middle = seq_len/2;
    if (nconf == 1) {
	for ( i = 0, j = seq_len-1; i < middle; i++, j--) {
	    temp = complementary_base [ (unsigned char) seq[i] ];
	    seq[i] = complementary_base [ (unsigned char) seq[j] ];
	    seq[j] = temp;
	    temp = conf[i];
	    conf[i] = conf[j];
	    conf[j] = temp;
	}
    } else if (nconf == 4) {
	for ( i = 0, j = seq_len-1; i < middle; i++, j--) {
	    temp = complementary_base [ (unsigned char) seq[i] ];
	    seq[i] = complementary_base [ (unsigned char) seq[j] ];
	    seq[j] = temp;
	    t[0] = conf[i*4+0];
	    t[1] = conf[i*4+1];
	    t[2] = conf[i*4+2];
	    t[3] = conf[i*4+3];
	    conf[i*4+0] = conf[j*4+3];
	    conf[i*4+1] = conf[j*4+2];
	    conf[i*4+2] = conf[j*4+1];
	    conf[i*4+3] = conf[j*4+0];
	    conf[j*4+0] = t[3];
	    conf[j*4+1] = t[2];
	    conf[j*4+2] = t[1];
	    conf[j*4+3] = t[0];
	}
    } else {
	fprintf(stderr, "Unsupported number of confidence values per base\n");
    }

    if ( seq_len % 2 )
      seq[middle] = complementary_base [ (unsigned char) seq[middle] ];
}

seq_t *dup_seq(seq_t *s) {
    size_t len = sizeof(seq_t) - sizeof(char *) + sequence_extra_len(s);
    seq_t *d = (seq_t *)calloc(1, len);


    memcpy(d, s, len);
    sequence_reset_ptr(d);

    /* Dup the annotation */
    if (s->anno) {
	d->anno = ArrayCreate(sizeof(int), ArrayMax(s->anno));
	memcpy(ArrayBase(int, d->anno),
	       ArrayBase(int, s->anno),
	       ArrayMax(s->anno) * sizeof(int));
    }

    return d;
}

void complement_seq_t(seq_t *s) {
    int tmp, alen;

    alen = ABS(s->len);
    complement_seq_conf(s->seq, s->conf, alen,
			s->format == SEQ_FORMAT_CNF4 ? 4 : 1);
    s->len *= -1;

    tmp = s->left;
    s->left  = alen - (s->right-1);
    s->right = alen - (tmp-1);
}

GRec sequence_index_query(GapIO *io, char *name) {
    return io->iface->seq.index_query(io->dbh, name);
}

int sequence_index_update(GapIO *io, char *name, int name_len, GRec rec) {
    char n2[1024];
    GRec r;
    //sprintf(n2, "%.*s", name_len, name);
    strncpy(n2, name, name_len > 1024 ? 1024 : name_len);
    n2[name_len > 1024 ? 1024 : name_len] = 0;

    r = io->iface->seq.index_add(io->dbh, n2, rec);
    if (r == -1)
	return -1;

    if (r != io->db->seq_name_index) {
	io->db = cache_rw(io, io->db);
	io->db->seq_name_index = r;
    }

    return 0;
}

/*
 * Finds the contig number and position of a sequence record number.
 *
 * If non-NULL r_out is filled with the associated range_t struct.
 *
 * If non-NULL s_out is filled with a pointer to the seq_t struct.
 * This will have had cache_incr() run on it, so the caller should
 * use cache_decr() to permit deallocation.
 */
int sequence_get_position2(GapIO *io, GRec snum, int *contig,
			   int *start, int *end, int *orient,
			   range_t *r_out, seq_t **s_out) {
    return bin_get_item_position(io, GT_Seq, snum,
				 contig, start, end, orient, NULL,
				 r_out, (void **)s_out);
}

int sequence_get_position(GapIO *io, GRec snum, int *contig,
			  int *start, int *end, int *orient) {
    return bin_get_item_position(io, GT_Seq, snum,
				 contig, start, end, orient, NULL,
				 NULL, NULL);
}

/*
 * Invalidates the cached consensus for this sequence.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_invalidate_consensus(GapIO *io, seq_t *s) {
    int start, end, contig;
    int i, nr;
    rangec_t *r;
    contig_t *c;

    if (io->read_only)
	return -1;

    if (-1 == sequence_get_position(io, s->rec, &contig, &start, &end, NULL))
	return -1;

    if (NULL == (c = (contig_t *)cache_search(io, GT_Contig, contig)))
	return -1;
    
    r = contig_bins_in_range(io, &c, start, end,
			     CSIR_LEAVES_ONLY, CONS_BIN_SIZE, &nr);

    for (i = 0; i < nr; i++) {
	bin_index_t *bin = (bin_index_t *)cache_search(io, GT_Bin, r[i].rec);
	if (!bin)
	    return -1;

	bin = cache_rw(io, bin);
	bin->flags |=  BIN_BIN_UPDATED;
	bin->flags &= ~BIN_CONS_VALID;
    }

    if (r)
	free(r);

    return 0;
}


/*
 * Given the record number for a sequence this returns the record
 * number for the contig containing it.
 */
int sequence_get_contig(GapIO *io, GRec snum) {
    bin_index_t *bin;
    int bnum;
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, snum);

    /* Bubble up bins until we hit the root */
    for (bnum = s->bin; bnum; bnum = bin->parent) {
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
	if (bin->parent_type != GT_Bin)
	    break;
    }

    assert(bin->parent_type == GT_Contig);
    return bin->parent;
}

/*
 * As per sequence_get_contig, but returns only the relative orientation of
 * this sequence vs the contig.
 */
int sequence_get_orient(GapIO *io, GRec snum) {
    bin_index_t *bin;
    int bnum;
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, snum);
    int comp = s->len < 0;

    /* Bubble up bins until we hit the root */
    for (bnum = s->bin; bnum; bnum = bin->parent) {
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
	if (bin->flags & BIN_COMPLEMENTED)
	    comp ^= 1;
	if (bin->parent_type != GT_Bin)
	    break;
    }

    assert(bin->parent_type == GT_Contig);
    return comp;
}

/*
 * Given a sequence struct, this returns the record number for other end,
 * if paired, or zero if not.
 * Returns -1 on failure.
 */
int sequence_get_pair(GapIO *io, seq_t *s) {
    bin_index_t *b;
    range_t *r;

    /* Get range struct for this seq */
    if (!s->bin)
	return -1;
    if (NULL == (b = (bin_index_t *)cache_search(io, GT_Bin, s->bin)))
	return -1;
    if (!b->rng) {
	cache_decr(io, b);
	return -1;
    }

    /* Jump over to pair */
    r = arrp(range_t, b->rng, s->bin_index);
    assert(r->rec == s->rec);
    cache_decr(io, b);
    return r->pair_rec;
}

/*
 * ---------------------------------------------------------------------------
 * Base editing functions
 */

int sequence_orient_pos(GapIO *io, seq_t **s, int pos, int *comp) {
    int swapped;
    sequence_get_position(io, (*s)->rec, NULL, NULL, NULL, &swapped);

    if (((*s)->len > 0) ^ swapped) {
	swapped = 0;
    } else {
	pos = ABS((*s)->len)-1 - pos;
	swapped = 1;
    }

    if (comp)
	*comp = swapped;

    return pos;
}

#define MAX4(ip)         \
((ip)[0] > (ip)[1]       \
 ? ((ip)[0] > (ip)[2]    \
    ? ((ip)[0] > (ip)[3] \
       ? (ip)[0]         \
       : (ip)[3])        \
    : ((ip)[2] > (ip)[3] \
       ? (ip)[2]         \
       : (ip)[3]))       \
 : ((ip)[1] > (ip)[2]    \
    ? ((ip)[1] > (ip)[3] \
       ? (ip)[1]         \
       : (ip)[3])        \
    : ((ip)[2] > (ip)[3] \
       ? (ip)[2]         \
       : (ip)[3])))


/*
 * Operates on position 'pos' in the displayed orientation rather than the
 * stored orientation.
 */
int sequence_get_base(GapIO *io, seq_t **s, int pos, char *base, int *conf,
		      int *cutoff, int contig_orient) {
    seq_t *n = *s;
    int comp = 0;

    if (pos < 0 || pos >= ABS(n->len))
	return -1;

    if (contig_orient)
	pos = sequence_orient_pos(io, s, pos, &comp);

    if (base) {
	if (comp)
	    *base = complementary_base[(unsigned char)n->seq[pos]];
	else
	    *base = n->seq[pos];
    }
    if (conf) {
	if (n->format != SEQ_FORMAT_CNF4) {
	    *conf = n->conf[pos];
	} else {
	    *conf = MAX4(&n->conf[pos*4]);
	}
    }
    if (cutoff) {
	if (pos < n->left-1 || pos >= n->right)
	    *cutoff = 1;
	else
	    *cutoff = 0;
    }

    return 0;
}

static double logodds2log[256];
static double *lo2l = logodds2log+128;

static double logodds2remainder[256];
static double *lo2r = logodds2remainder+128;

static unsigned char logodds2phred[256];
static unsigned char *lo2ph = logodds2phred+128;

static int lookup_init = 0;

/*
 * Operates on position 'pos' in the displayed orientation rather than the
 * stored orientation.
 *
 * As per sequence_get_base but conf is an array of 4 confidence values
 * returned as log(probabilities).
 * If only one is present it is taken to be a phred score and we compute
 * the remainder using (1-p)/3. Otherwise we assume we store 4 log-odds
 * scores instead.
 */
int sequence_get_base4(GapIO *io, seq_t **s, int pos, char *base, double *conf,
		       int *cutoff, int contig_orient) {
    seq_t *n = *s;
    int comp = 0;

    if (!lookup_init) {
	int i;
	lookup_init = 1;

	/* Log odds value to log(P) */
	for (i = -128; i < 128; i++) {
	    double p = 1 / (1 + pow(10, -i / 10.0));
	    lo2l[i] = log(p);
	    lo2r[i] = log((1-p)/3);
	    lo2ph[i] = 10*log(1+pow(10, i/10.0))/log(10)+0.4999;
	}
    }

    if (pos < 0 || pos >= ABS(n->len))
	return -1;

    if (contig_orient)
	pos = sequence_orient_pos(io, s, pos, &comp);

    if (base) {
	if (comp)
	    *base = complementary_base[(unsigned char)n->seq[pos]];
	else
	    *base = n->seq[pos];
    }

    if (conf) {
	if (n->format != SEQ_FORMAT_CNF4) {
	    switch(n->seq[pos]) {
	    case 'A': case 'a':
		conf[0] = lo2l[n->conf[pos]];
		conf[1] = lo2r[n->conf[pos]];
		conf[2] = lo2r[n->conf[pos]];
		conf[3] = lo2r[n->conf[pos]];
		break;
	    case 'C': case 'c':
		conf[0] = lo2r[n->conf[pos]];
		conf[1] = lo2l[n->conf[pos]];
		conf[2] = lo2r[n->conf[pos]];
		conf[3] = lo2r[n->conf[pos]];
		break;
	    case 'G': case 'g':	
		conf[0] = lo2r[n->conf[pos]];
		conf[1] = lo2r[n->conf[pos]];
		conf[2] = lo2l[n->conf[pos]];
		conf[3] = lo2r[n->conf[pos]];
		break;
	    case 'T': case 't':
		conf[0] = lo2r[n->conf[pos]];
		conf[1] = lo2r[n->conf[pos]];
		conf[2] = lo2r[n->conf[pos]];
		conf[3] = lo2l[n->conf[pos]];
		break;
	    default:
		conf[0] = 0;
		conf[1] = 0;
		conf[2] = 0;
		conf[3] = 0;
		break;
	    }
	} else {
	    int i;
	    for (i = 0; i < 4; i++) {
		double p = n->conf[pos*4+i];
		p /= 10.0;
		conf[i] = p*log(10) - log(1+pow(10, p));
	    }
	}

	//richard_munge_conf(conf);
	//rob_munge_conf(conf);

	if (comp) {
	    double tmp;
	    tmp = conf[0]; conf[0] = conf[3]; conf[3] = tmp;
	    tmp = conf[1]; conf[1] = conf[2]; conf[2] = tmp;
	}
    }

    if (cutoff) {
	if (pos < n->left || pos > n->right)
	    *cutoff = 1;
	else
	    *cutoff = 0;
    }

    return 0;
}


int sequence_replace_base(GapIO *io, seq_t **s, int pos, char base, int conf,
			  int contig_orient) {
    seq_t *n;
    int comp = 0;

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    if (pos < 0 || pos >= ABS(n->len))
	return -1;

    sequence_invalidate_consensus(io, n);

    if (contig_orient)
	pos = sequence_orient_pos(io, s, pos, &comp);

    if (n->format != SEQ_FORMAT_CNF4) {
	if (comp) {
	    n->seq[pos] = complementary_base[(unsigned char)base];
	    n->conf[pos] = conf;
	} else {
	    n->seq[pos] = base;
	    n->conf[pos] = conf;
	}
    } else {
	double remainder = -4.34294482*log(2+3*pow(10, conf/10.0));

	n->seq[pos] = comp ? complementary_base[(unsigned char)base] : base;

	switch(base) {
	case 'A': case 'a':
	    n->conf[pos*4+0] = comp ? remainder : conf;
	    n->conf[pos*4+1] = remainder;
	    n->conf[pos*4+2] = remainder;
	    n->conf[pos*4+3] = comp ? conf : remainder;
	    break;
	case 'C': case 'c':
	    n->conf[pos*4+0] = remainder;
	    n->conf[pos*4+1] = comp ? remainder : conf;
	    n->conf[pos*4+2] = comp ? conf : remainder;
	    n->conf[pos*4+3] = remainder;
	    break;
	case 'G': case 'g':
	    n->conf[pos*4+0] = remainder;
	    n->conf[pos*4+1] = comp ? conf : remainder;
	    n->conf[pos*4+2] = comp ? remainder : conf;
	    n->conf[pos*4+3] = remainder;
	    break;
	case 'T': case 't':
	    n->conf[pos*4+0] = comp ? conf : remainder;
	    n->conf[pos*4+1] = remainder;
	    n->conf[pos*4+2] = remainder;
	    n->conf[pos*4+3] = comp ? remainder : conf;
	    break;
	default:
	    n->conf[pos*4+0] = -5;
	    n->conf[pos*4+1] = -5;
	    n->conf[pos*4+2] = -5;
	    n->conf[pos*4+3] = -5;
	    break;
	}
    }

    return 0;
}

/*
 * Insert into position 'pos' of a sequence, where pos starts from 0
 * (before first base).
 */
int sequence_insert_base(GapIO *io, seq_t **s, int pos, char base, char conf,
			 int contig_orient) {
    seq_t *n;
    int comp = 0;
    size_t extra_len = sequence_extra_len(*s) + 1 + sequence_conf_size(*s);
    char *c_old;

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    sequence_invalidate_consensus(io, n);

    n = cache_item_resize(n, sizeof(*n) + extra_len);
    if (NULL == n)
	return -1;
    *s = n;

    if (contig_orient) {
	pos = sequence_orient_pos(io, &n, pos, &comp);
	if (comp)
	    pos++;
    } else {
	pos = n->len < 0
	    ? -n->len - pos
	    : pos;
    }

    /* Reset internal pointers assuming new length */
    if (n->len < 0)
	n->len--;
    else
	n->len++;
    c_old = n->conf;
    sequence_reset_ptr(n);

    /* Shift */
    memmove(&n->seq[pos+1], &n->seq[pos],
	    extra_len - ((char *)&n->seq[pos] - (char *)&n->data));

    c_old++;
    memmove(&c_old[sequence_conf_size(n)*pos]+1,
	    &n->conf[sequence_conf_size(n)*pos],
	    extra_len - ((char *)&n->conf[sequence_conf_size(n)*(pos)]+1
			 - (char *)&n->data));
	    
    /* Set */
    n->seq[pos] = comp ? complementary_base[(unsigned char)base] : base;
    if (n->format == SEQ_FORMAT_CNF4) {
	double remainder = -4.34294482*log(2+3*pow(10, conf/10.0));
	switch (base) {
	case 'A': case 'a':
	    n->conf[pos*4+0] = comp ? remainder : conf;
	    n->conf[pos*4+1] = remainder;
	    n->conf[pos*4+2] = remainder;
	    n->conf[pos*4+3] = comp ? conf : remainder;
	    break;
	case 'C': case 'c':
	    n->conf[pos*4+0] = remainder;
	    n->conf[pos*4+1] = comp ? remainder : conf;
	    n->conf[pos*4+2] = comp ? conf : remainder;
	    n->conf[pos*4+3] = remainder;
	    break;
	case 'G': case 'g':
	    n->conf[pos*4+0] = remainder;
	    n->conf[pos*4+1] = comp ? conf : remainder;
	    n->conf[pos*4+2] = comp ? remainder : conf;
	    n->conf[pos*4+3] = remainder;
	    break;
	case 'T': case 't':
	    n->conf[pos*4+0] = comp ? conf : remainder;
	    n->conf[pos*4+1] = remainder;
	    n->conf[pos*4+2] = remainder;
	    n->conf[pos*4+3] = comp ? remainder : conf;
	    break;
	default:
	    n->conf[pos*4+0] = -5;
	    n->conf[pos*4+1] = -5;
	    n->conf[pos*4+2] = -5;
	    n->conf[pos*4+3] = -5;
	    break;
	}
    } else {
	n->conf[pos] = conf;
    }

    if (pos < n->left-1)
	n->left++;

    if (pos <= n->right)
	n->right++;

    return 0;
}

/*
 * Delete position 'pos' of a sequence, where pos starts from 0
 */
int sequence_delete_base(GapIO *io, seq_t **s, int pos, int contig_orient) {
    seq_t *n;
    int alen;
    int comp = 0;
    size_t extra_len = sequence_extra_len(*s);

    if (pos >= ABS((*s)->len) || pos < 0)
	return 0;

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    sequence_invalidate_consensus(io, n);

    if (contig_orient) {
	pos = sequence_orient_pos(io, &n, pos, &comp);
    } else {
	pos = n->len < 0
	    ? -n->len - pos -1
	    : pos;
    }

    if (n->len < 0)
	n->len++;
    else
	n->len--;

    if (pos < n->left-1)
	n->left--;

    if (pos < n->right)
	n->right--;

    if (pos >= ABS(n->len) || pos < 0) {
	sequence_reset_ptr(n);
	return 0;
    }

    /* Shift */
    memmove(&n->conf[sequence_conf_size(n) * pos],
	    &n->conf[sequence_conf_size(n) * (pos+1)],
	    extra_len -
	    ((char *)&n->conf[sequence_conf_size(n)*(pos+1)]
	             - (char *)&n->data));
    extra_len -= sequence_conf_size(n);

    memmove(&n->seq[pos], &n->seq[pos+1],
	    extra_len - ((char *)&n->seq[pos+1] - (char *)&n->data));

    sequence_reset_ptr(n);

    return 0;
}

/*
 * Moves all annotations attached to a sequence left or right by a certain
 * amount 'dist'. If dist is negative the annotations move left, otherwise
 * they move right.
 *
 * This is used when moving sequences within the editor, just prior to the
 * move itself as we only look for annotations covering the current
 * coordinates.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_move_annos(GapIO *io, seq_t **s, int dist) {
    contig_t *c;
    int start, end, orient, contig;
    rangec_t *r;
    int nr, i;

    /* Find the position and contig this sequence currently covers */
    cache_incr(io, *s);
    if (0 != sequence_get_position(io, (*s)->rec,
				   &contig, &start, &end, &orient))
	return -1;


    /*
     * Identify annotations spanning this region.
     *
     * NB: We may want a specialist function for this if we attempt moving
     * a really long sequence or a sequence in a really deep region.
     *
     * The way to implement this would be along the lines of
     * contig_seqs_in_range2(). We could then keep track of whether the
     * moved item would still be within the same bin and if so we just
     * update the coordinate, otherwise we fall back on the extract and
     * reinsert method (taking care to note we don't then move the same
     * record multiple times due to it now being in a bin we're about to
     * visit). Alternatively we could store bin record and bin_index in
     * the generated rangec_t struct so we can call anno_in_range and
     * take short cuts here if we can get away with it.
     */
    c = cache_search(io, GT_Contig, contig);
    r = (rangec_t *)contig_anno_in_range(io, &c, start-1, end+1, 0, &nr);

    /* Figure out what to move */
    for (i = 0; i < nr; i++) {
	range_t R, *R_out;
	anno_ele_t *a;
	bin_index_t *bin;

	assert((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO);

	if (r[i].pair_rec != (*s)->rec)
	    continue;

	/* A tag for this seq, so create a new range_t and del/add it */
	bin_remove_item(io, &c, GT_AnnoEle, r[i].rec);
	R.start    = r[i].start + dist;
	R.end      = r[i].end   + dist;
	R.rec      = r[i].rec;
	R.mqual    = r[i].mqual;
	R.pair_rec = r[i].pair_rec;
	R.flags    = r[i].flags;
	bin = bin_add_range(io, &c, &R, &R_out, NULL);
	cache_incr(io, bin);

	/* With luck the bin & index into bin haven't changed */
	a = cache_search(io, GT_AnnoEle, r[i].rec);
	if (a->bin != bin->rec ||
	    a->idx != R_out - ArrayBase(range_t, bin->rng)) {
	    //printf("New tag bin %d->%d %d->%d\n", a->bin, bin->rec,
	    //	   a->idx, R_out - ArrayBase(range_t, bin->rng));
	    a = cache_rw(io, a);
	    a->bin = bin->rec;
	    a->idx = R_out - ArrayBase(range_t, bin->rng);
	}

	cache_decr(io, bin);
    }

    free(r);

    cache_decr(io, *s);
}
