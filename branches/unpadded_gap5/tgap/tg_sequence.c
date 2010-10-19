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
    s->alignment = (unsigned char *)s->trace_name + s->trace_name_len + 1;
    s->seq = (char *)s->alignment + s->alignment_len + 1;
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
	s->alignment_len                            + 1 +
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
    tg_rec rec;
    int idx;
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

    strcpy((char *)s->alignment, f->alignment ? (char *)f->alignment : "");
    s->alignment_len = strlen((char *)s->alignment);

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
tg_rec sequence_new_from(GapIO *io, seq_t *s) {
    tg_rec rec;
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
    strcpy(cp, (char *)n->alignment);
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
    strcpy(cp, (char *)n->alignment);
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
int seq_pos(GapIO *io, tg_rec rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_pos(&s);
}

int seq_len(GapIO *io, tg_rec rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_len(&s);
}

int seq_left(GapIO *io, tg_rec rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_left(&s);
}

int seq_right(GapIO *io, tg_rec rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_right(&s);
}

int seq_mapping_qual(GapIO *io, tg_rec rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_mapping_qual(&s);
}

char *seq_name(GapIO *io, tg_rec rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_name(&s);
}

char *seq_seq(GapIO *io, tg_rec rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_seq(&s);
}

char *seq_conf(GapIO *io, tg_rec rec) {
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
	for ( i = 0, j = seq_len-1; i < j; i++, j--) {
	    temp = (unsigned char) seq[i];
	    seq[i] = complementary_base [ (unsigned char) seq[j] ];
	    seq[j] = complementary_base [ temp ];
	    temp = conf[i];
	    conf[i] = conf[j];
	    conf[j] = temp;
	}
    } else if (nconf == 4) {
	for ( i = 0, j = seq_len-1; i < j; i++, j--) {
	    temp = (unsigned char) seq[i];
	    seq[i] = complementary_base [ (unsigned char) seq[j] ];
	    seq[j] = complementary_base [ temp ];
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

    /* Also complement the alignment - just reverse it */
    if (s->alignment_len > 2) {
	int i,j;

	for (i = 0, j = s->alignment_len-2; i < j; i+=2, j-=2){
	    int tmp1, tmp2;
	    tmp1 = s->alignment[i  ];
	    tmp2 = s->alignment[i+1];
	    s->alignment[i  ] = s->alignment[j  ];
	    s->alignment[i+1] = s->alignment[j+1];
	    s->alignment[j  ] = tmp1;
	    s->alignment[j+1] = tmp2;
	}
    }
}

tg_rec sequence_index_query(GapIO *io, char *name) {
    return io->iface->seq.index_query(io->dbh, name, 0);
}

tg_rec sequence_index_query_prefix(GapIO *io, char *prefix) {
    return io->iface->seq.index_query(io->dbh, prefix, 1);
}

tg_rec *sequence_index_query_all(GapIO *io, char *name, int prefix,
				 int *nrecs) {
    return io->iface->seq.index_query_all(io->dbh, name, prefix, nrecs);
}

int sequence_index_update(GapIO *io, char *name, int name_len, tg_rec rec) {
    char n2[1024];
    tg_rec r;
    //sprintf(n2, "%.*s", name_len, name);
    strncpy(n2, name, name_len > 1024 ? 1024 : name_len);
    n2[name_len > 1024 ? 1024 : name_len] = 0;

    r = io->iface->seq.index_add(io->dbh, n2, rec);
    if (r == -1)
	return -1;

    if (r != io->db->seq_name_index) {
	io->db = cache_rw(io, io->db);
	io->db->seq_name_index = (GCardinal)r;
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
int sequence_get_position2(GapIO *io, tg_rec snum, tg_rec *contig,
			   int *start, int *end, int *orient,
			   range_t *r_out, seq_t **s_out) {
    return bin_get_item_position(io, GT_Seq, snum,
				 contig, start, end, orient, NULL,
				 r_out, (void **)s_out);
}

int sequence_get_position(GapIO *io, tg_rec snum, tg_rec *contig,
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
    int start, end;
    tg_rec contig;

    if (io->read_only)
	return -1;

    if (-1 == sequence_get_position(io, s->rec, &contig, &start, &end, NULL))
	return -1;

    return bin_invalidate_consensus(io, contig, start, end);
}


/*
 * Given the record number for a sequence this returns the record
 * number for the contig containing it.
 */
tg_rec sequence_get_contig(GapIO *io, tg_rec snum) {
    bin_index_t *bin = NULL;
    tg_rec bnum;
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, snum);

    /* Bubble up bins until we hit the root */
    for (bnum = s->bin; bnum; bnum = bin->parent) {
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
	if (bin->parent_type != GT_Bin)
	    break;
    }

    assert(bin && bin->parent_type == GT_Contig);
    return bin->parent;
}

/*
 * As per sequence_get_contig, but returns only the relative orientation of
 * this sequence vs the contig.
 */
int sequence_get_orient(GapIO *io, tg_rec snum) {
    bin_index_t *bin = NULL;
    tg_rec bnum;
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

    assert(bin && bin->parent_type == GT_Contig);
    return comp;
}

/*
 * Replaces the cigar string with a new one. Checks validity of it too.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_cigar(GapIO *io, seq_t **s, unsigned char *cig, int len) {
    size_t extra_len;
    seq_t *n;
    char *old_seq, *new_seq;
    char *old_conf, *new_conf;
    int grow;

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    /* Realloc if needed, don't care if it shrunk */
    if (len > n->alignment_len) {
	extra_len = sequence_extra_len(*s);
	n = cache_item_resize(n, sizeof(*n) + extra_len);
	if (NULL == n)
	    return -1;
	*s = n;
    }

    /* Shift other data */
    if (n->alignment_len != len) {
	sequence_reset_ptr(n);
	old_seq  = n->seq;
	old_conf = n->conf;
	grow = len > n->alignment_len ? 1 : 0;

	n->alignment_len = len;
	sequence_reset_ptr(n);	
	new_seq  = n->seq;
	new_conf = n->conf;

	if (grow) {
	    memmove(new_conf, old_conf, ABS(n->len) * sequence_conf_size(n));
	    memmove(new_seq,  old_seq,  ABS(n->len));
	} else {
	    memmove(new_seq,  old_seq,  ABS(n->len));
	    memmove(new_conf, old_conf, ABS(n->len) * sequence_conf_size(n));
	}
    }

    /* Copy in the new cig */
    n->alignment_len = len;
    memcpy(n->alignment, cig, len);
    n->alignment[len] = 0;

    return 0;
}

/*
 * Given a sequence struct, this returns the record number for other end,
 * if paired, or zero if not.
 * Returns -1 on failure.
 */
tg_rec sequence_get_pair(GapIO *io, seq_t *s) {
    bin_index_t *b;
    range_t *r;

    /* Get range struct for this seq */
    if (!s->bin)
	return -1;
    if (NULL == (b = (bin_index_t *)cache_search(io, GT_Bin, s->bin)))
	return -1;
    if (!b->rng)
	return -1;

    /* Jump over to pair */
    r = arrp(range_t, b->rng, s->bin_index);
    assert(r->rec == s->rec);
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

/*
 * Operates on position 'pos' in the displayed orientation rather than the
 * stored orientation.
 * As above, but pos is an unpadded coordinate into the sequence so we
 * need to parse the CIGAR string too to figure out what base number this
 * is.
 */
int sequence_get_ubase(GapIO *io, seq_t **s, int pos, int nth,
		       char *base, int *conf, int *cutoff) {
    seq_t *n = *s;
    int comp;

    if (pos < 0)
	return -1;

    sequence_get_position(io, n->rec, NULL, NULL, NULL, &comp);
    comp = (n->len > 0) ^ comp ? 0 : 1;

    /* Convert pos from unpadded to in-seq offset */
    if (n->alignment) {
	int op = 0, op_len = 0, spos = 0;
	int more = 1;
	int a_ind = 0;

	/* Find pos first */
	pos++;
	while (pos > 0 && a_ind < n->alignment_len) {
	    int r_ind = comp ? n->alignment_len - a_ind - 2 : a_ind;
	    op_len = n->alignment[r_ind];
	    op = n->alignment[r_ind+1];
	    a_ind += 2;

	    //printf("%d%c\n", op_len, op);
	    switch (op) {
	    case 'H':
		continue;

	    case 'S':
	    case 'M':
		while (op_len && pos) {
		    spos++;
		    op_len--;
		    pos--;
		}
		break;

	    case 'D':
		while (op_len && pos) {
		    op_len--;
		    pos--;
		}
		break;

	    case 'P':
		op_len = 0;
		break;

	    case 'I':
		spos += op_len;
		op_len = 0;
		break;
	    }
	}
	spos--;
	
	//printf("op_len=%d pos=%d nth=%d spos=%d %.10s\n",
	//       op_len, pos, nth, spos, &n->seq[spos]);

	if (pos != 0 && a_ind == n->alignment_len) {
	    /* Off end */
	    return -1;
	}

	/* And now the nth base at pos */
	while (more && nth > 0 && a_ind < n->alignment_len) {
	    if (!op_len) {
		int r_ind = comp ? n->alignment_len - a_ind - 2 : a_ind;
		op_len = n->alignment[r_ind];
		op = n->alignment[r_ind+1];
		a_ind += 2;
	    }
	    switch (op) {
	    case 'H':
		continue;

	    case 'S':
	    case 'M':
	    case 'D':
		more = 0;
		break;

	    case 'P':
		while (op_len && nth) {
		    op_len--;
		    nth--;
		}
		break;

	    case 'I':
		while (op_len && nth) {
		    op_len--;
		    nth--;
		    spos++;
		}
		break;
	    }
	}

	/* CIGAR counts from right for complemented data */
	if (comp)
	    spos = ABS(n->len)-1 - spos;

	/* If still have 'nth' left, it's inbetween spos and spos+1 */
	if (nth || op == 'D' || op == 'P') {
	    if (base)
		*base = nth ? '-' : '*';
	    if (conf) {
		if (spos + 1 < ABS(n->len))
		    *conf = (n->conf[spos] + n->conf[spos+1])/2; // or MIN()?
		else
		    *conf = n->conf[spos];
	    }
	    return 0;
	}

	pos = spos;
    }

    /*
     * Otherwise it's nul alignment (all match) or in an alignment at
     * reference coordinate base.
     */
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


/*
 * Given the nth base at unpadded reference position pos, convert this
 * coordinate to a raw sequence position.
 *
 * Optionally fills out 'exists' fill a flag to indicate if the base
 * is in the raw sequence or not (pos may refer to a deletion or nth may
 * refer to an insertion not within this sequence). Either was the value
 * returned is always the base at pos or to the left of it if it doesn't
 * exist.
 *
 * Returns sequence position on success
 *        -1 on failure
 */
int sequence_get_spos(GapIO *io, seq_t **s, int pos, int nth, int *exists) {
    seq_t *n = *s;
    int comp;

    if (pos < 0)
	return -1;

    sequence_get_position(io, n->rec, NULL, NULL, NULL, &comp);
    comp = (n->len > 0) ^ comp ? 0 : 1;

    if (exists)
	*exists = 1;

    /* Same as in sequence_get_ubase... */

    /* Convert pos from unpadded to in-seq offset */
    if (n->alignment) {
	int op = 0, op_len = 0, spos = 0;
	int more = 1;
	int a_ind = 0;

	/* Find pos first */
	pos++;
	while (pos > 0 && a_ind < n->alignment_len) {
	    int r_ind = comp ? n->alignment_len - a_ind - 2 : a_ind;
	    op_len = n->alignment[r_ind];
	    op = n->alignment[r_ind+1];
	    a_ind += 2;

	    //printf("%d%c\n", op_len, op);
	    switch (op) {
	    case 'H':
		continue;

	    case 'S':
	    case 'M':
		while (op_len && pos) {
		    spos++;
		    op_len--;
		    pos--;
		}
		break;

	    case 'D':
		while (op_len && pos) {
		    op_len--;
		    pos--;
		}
		break;

	    case 'P':
		op_len = 0;
		break;

	    case 'I':
		spos += op_len;
		op_len = 0;
		break;
	    }
	}
	spos--;
	
	//printf("op_len=%d pos=%d nth=%d spos=%d %.10s\n",
	//       op_len, pos, nth, spos, &n->seq[spos]);

	if (pos != 0 && a_ind == n->alignment_len) {
	    /* Off end */
	    return -1;
	}

	/* And now the nth base at pos */
	while (more && nth > 0 && a_ind < n->alignment_len) {
	    if (!op_len) {
		int r_ind = comp ? n->alignment_len - a_ind - 2 : a_ind;
		op_len = n->alignment[r_ind];
		op = n->alignment[r_ind+1];
		a_ind += 2;
	    }
	    switch (op) {
	    case 'H':
		continue;

	    case 'S':
	    case 'M':
	    case 'D':
		more = 0;
		break;

	    case 'P':
		while (op_len && nth) {
		    op_len--;
		    nth--;
		}
		break;

	    case 'I':
		while (op_len && nth) {
		    op_len--;
		    nth--;
		    spos++;
		}
		break;
	    }
	}

	/* CIGAR counts from right for complemented data */
	if (comp)
	    spos = ABS(n->len)-1 - spos;

	/* If still have 'nth' left, it's inbetween spos and spos+1 */
	if (nth || op == 'D' || op == 'P') {
	    if (exists)
		*exists = 0;
	}

	pos = spos;
    }

    return pos;
}


/*
 * Converts a raw position within a sequence to the unpadded coordinate and
 * an nth base.
 * Ie a sequence with CIGAR 10M5I10M and raw_pos of 12 would return
 * upos=10, nth=2.  raw_pos = 17 => upos=12, nth=0.
 * For deletions: 10M5D10M and raw_pos 10=>{10,0}, raw_pos 11=>{16,0}.
 * These are regardless of the start position of the read (FIXME: or
 * the orientation?).
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_get_upos(GapIO *io, seq_t **s, int raw_pos, int *upos, int *unth) {
    seq_t *n = *s;
    int comp;

    if (raw_pos < 0)
	return -1;

    sequence_get_position(io, n->rec, NULL, NULL, NULL, &comp);
    comp = (n->len > 0) ^ comp ? 0 : 1;

    //    if (contig_orient)
    //	pos = sequence_orient_pos(io, s, pos, &comp);

    if (n->alignment) {
	int op, op_len, spos = 0, pos = 0, nth = 0, a_ind = 0;
	
	do {
	    int r_ind = comp ? n->alignment_len - a_ind - 2 : a_ind;
	    op_len = n->alignment[r_ind];
	    op = n->alignment[r_ind+1];
	    a_ind += 2;

	    switch (op) {
	    case 'H':
		continue;

	    case 'S':
	    case 'M':
		while (op_len && raw_pos != spos) {
		    op_len--;
		    spos++;
		    pos++;
		    nth = 0;
		}
		break;

	    case 'D':
		pos += op_len;
		nth = 0;
		op_len = 0;
		break;

	    case 'P':
		while(op_len) {
		    op_len--;
		    nth++;
		}
		break;

	    case 'I':
		while(op_len && raw_pos != spos) {
		    op_len--;
		    spos++;
		    nth++;
		}
		break;
	    }
	} while (spos != raw_pos && a_ind < n->alignment_len);

	if (spos == raw_pos) {
	    *upos = pos;
	    *unth = nth;
	    //printf("Get upos %d => %d\n", raw_pos, pos);
	} else {
	    return -1;
	}
    } else {
	*upos = raw_pos;
	*unth = 0;
    }

    return 0;
}

/*
 * Returns the padded length of a sequence, minus the impact of
 * insertions via other sequences.
 */
int sequence_padded_len(seq_t *s) {
    int i, len;

    if (!s->alignment || s->alignment_len == 0)
	return ABS(s->len);

    for (len = i = 0; i < s->alignment_len; i += 2) {
	switch (s->alignment[i+1]) {
	case 'S':
	case 'M':
	case 'I':
	    len += s->alignment[i];
	}
    }

    return len;
}

/*
 * Returns the unpadded length of a sequence.
 */
int sequence_unpadded_len(seq_t *s, int *nth) {
    int i, len;

    if (nth)
	*nth = 0;

    if (!s->alignment || s->alignment_len == 0)
	return ABS(s->len);

    for (len = i = 0; i < s->alignment_len; i += 2) {
	switch (s->alignment[i+1]) {
	case 'S':
	case 'M':
	case 'D':
	    len += s->alignment[i];
	}
    }

    if (nth && s->alignment[s->alignment_len-1] == 'I')
	*nth = s->alignment[s->alignment_len-2];

    return len;
}

/*
 * Given a CIGAR string and an index & partially used counter into it, 
 * turn this into ref and sequence coords.
 *
 *    0 2 4 6 8
 *    | | | | |
 * Eg 3M2I5M2D6M
 *
 * index 4, count 3 => apply 3M2I3M and see what position this gives
 * us in the original sequence (3+2+3) and in the reference (3+3).
 *
 * Returns 0 for success and fills out rpos / spos (may be NULL ptrs)
 *        -1 for failure
 */
int sequence_cigar2pos(GapIO *io, seq_t *s, int cigar_ind, int cigar_len,
		       int *spos_p, int *snth_p, int *rpos_p, int *rnth_p) {
    int i;
    int comp;
    int spos = 0, rpos = 0;
    int snth = 0, rnth = 0;

    sequence_get_position(io, s->rec, NULL, NULL, NULL, &comp);
    comp = (s->len > 0) ^ comp ? 0 : 1;

    if (cigar_ind > s->alignment_len)
	return -1;

    for (i = 0; i < cigar_ind-2; i += 2) {
	int r_ind = comp ? s->alignment_len - i - 2 : i;
	int op_len = s->alignment[r_ind];
	int op = s->alignment[r_ind+1];

	switch (op) {
	case 'S':
	case 'M':
	    snth = rnth = 0;
	    spos += op_len;
	    rpos += op_len;
	    break;

	case 'P':
	    snth = 0;
	    rnth += op_len;
	    break;

	case 'I':
	    snth = 0;
	    rnth += op_len;
	    spos += op_len;
	    break;

	case 'D':
	    rnth = 0;
	    snth += op_len;
	    rpos += op_len;
	}
    }
    
    {
	int r_ind = comp ? s->alignment_len - i - 2 : i;
	int op = s->alignment[r_ind+1];
	int op_len = s->alignment[r_ind];

	cigar_len = op_len - cigar_len;

	switch (op) {
	case 'S':
	case 'M':
	    snth = rnth = 0;
	    spos += cigar_len;
	    rpos += cigar_len;
	    break;

	case 'P':
	    snth = 0;
	    rnth += cigar_len;
	    break;

	case 'I':
	    snth = 0;
	    rnth += cigar_len;
	    spos += cigar_len;
	    break;

	case 'D':
	    rnth = 0;
	    snth += cigar_len;
	    rpos += cigar_len;
	}
    }

    if (spos_p) *spos_p = spos;
    if (snth_p) *snth_p = snth;
    if (rpos_p) *rpos_p = rpos;
    if (rnth_p) *rnth_p = rnth;

    return 0;
}

/*
 * The reverse of sequence_ciga2pos.
 * This takes unpadded coord (pos,nth) and identifies which cigar operation
 * and offset into op this position corresponds to. It also returns 'spos',
 * the coordinate into the raw s->seq string (or we're between that and
 * spos+1).
 *
 * The cigar_len returned is the *remaining* length of this cigar operation
 * and not how far into it we currently are.
 *
 * nth_p is the remaining number of bases in the original nth specified.
 * If *nth_p is non-zero it implies we're at a position in this sequence
 * that has been inserted due to insertions from other aligned sequences
 * (implicitly known via the input "nth" value).
 *
 * Returns 0 on success
 *        -1 on failure (off end of seq)
 */
int sequence_pos2cigar(GapIO *io, seq_t *n, int pos, int nth,
		       int *cigar_ind, int *cigar_op, int *cigar_len,
		       int *spos_p, int *pos_p, int *nth_p, int *comp_p) {
    int op = 0, op_len = 0, spos = 0;
    int more = 1;
    int a_ind = 0, r_ind = 0;

    /* Based on sequence_get_ubase */
    /* FIXME: make sequence_get_ubase and sequence_get_spos call this
     * instead.
     */

    int comp;

    if (pos < 0)
	return -1;

    /* Or simulate <len>M here? */
    if (!n->alignment || n->alignment_len == 0)
	return -1;

    sequence_get_position(io, n->rec, NULL, NULL, NULL, &comp);
    comp = (n->len >= 0) ^ comp ? 0 : 1;


    /* Convert pos from unpadded to in-seq offset */
    /* Find pos first */
    pos++;
    while (pos > 0 && a_ind < n->alignment_len) {
	r_ind = comp ? n->alignment_len - a_ind - 2 : a_ind;
	op_len = n->alignment[r_ind];
	op = n->alignment[r_ind+1];
	a_ind += 2;

	//printf("%d%c\n", op_len, op);
	switch (op) {
	case 'H':
	    continue;

	case 'S':
	case 'M':
	    while (op_len && pos) {
		spos++;
		op_len--;
		pos--;
	    }
	    break;

	case 'D':
	    while (op_len && pos) {
		op_len--;
		pos--;
	    }
	    break;

	case 'P':
	    op_len = 0;
	    break;

	case 'I':
	    spos += op_len;
	    op_len = 0;
	    break;
	}
    }
    spos--;
	
    //printf("op_len=%d pos=%d nth=%d spos=%d %.10s\n",
    //       op_len, pos, nth, spos, &n->seq[spos]);

    if (pos != 0 && a_ind == n->alignment_len) {
	/* Off end */
	return -1;
    }

    /* And now the nth base at pos */
    while (more && nth > 0 && a_ind < n->alignment_len) {
	if (!op_len) {
	    r_ind = comp ? n->alignment_len - a_ind - 2 : a_ind;
	    op_len = n->alignment[r_ind];
	    op = n->alignment[r_ind+1];
	    a_ind += 2;
	}
	switch (op) {
	case 'H':
	    continue;

	case 'S':
	case 'M':
	case 'D':
	    more = 0;
	    break;

	case 'P':
	    while (op_len && nth) {
		op_len--;
		nth--;
	    }
	    break;

	case 'I':
	    while (op_len && nth) {
		op_len--;
		nth--;
		spos++;
	    }
	    break;
	}
    }

    /* CIGAR counts from right for complemented data */
    if (comp)
	spos = ABS(n->len)-1 - spos;

    if (cigar_ind) *cigar_ind = r_ind;
    if (cigar_op)  *cigar_op  = op;
    if (cigar_len) *cigar_len = op_len;
    if (pos_p)     *pos_p     = pos;
    if (nth_p)     *nth_p     = nth;
    if (spos_p)    *spos_p    = spos;
    if (comp_p)    *comp_p    = comp;

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

	/* Special case for manually edited bases */
	lo2l[100] = 0;    /* log(prob = 1.0) */
	lo2r[100] = -100; /* lof(prob ~ 0.0) */
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
 * Computes a padded form of the sequence in seq_p, with pos_p/nth_p pointing
 * to unpadded pos and nth base at pos for each base in seq.
 * On input, plen is the extra length to use when allocating our buffers.
 * On exit, plen is the padded length of the actual alignment.
 * The caller is expected to call free on pos, nth and seq.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int seq_padded(GapIO *io, seq_t *s, int **pos_p, int **nth_p,
	       char **seq_p, char **conf_p, int *plen, int *comp_p) {
    int l;
    char *seq, *cnf;
    int *nth;
    int *pos;
    int i, j, spos, ppos, rpos, rnth = 0;
    int comp;

    if (!pos_p || !nth_p || !seq_p || !conf_p)
	return -1;

    sequence_get_position(io, s->rec, NULL, NULL, NULL, &comp);
    comp = (s->len >= 0) ^ comp ? 0 : 1;

    if (comp_p)
	*comp_p = comp;

    for (l = i = 0; i < s->alignment_len; i+=2) {
	switch (s->alignment[i+1]) {
	case 'S':
	case 'M':
	case 'I':
	case 'P':
	case 'D':
	    l += s->alignment[i];
	}
    }
    seq = malloc(l+*plen);
    cnf = malloc(l+*plen);
    nth = malloc((l+*plen)*sizeof(int));
    pos = malloc((l+*plen)*sizeof(int));
    

    if (!seq || !nth || !pos || !cnf)
	return -1;

    *pos_p  = pos;
    *nth_p  = nth;
    *seq_p  = seq;
    *conf_p = cnf;
    *plen   = l;

    if (comp) {
	spos = ABS(s->len);
	for (rpos=-1, ppos = 0, i = s->alignment_len-2; i >= 0; i -= 2) {
	    int cigar_op  = s->alignment[i+1];
	    int cigar_len = s->alignment[i];

	    switch (cigar_op) {
	    case 'S':
	    case 'M':
		for (j = 0; j < cigar_len; j++) {
		    seq[ppos] = complementary_base[s->seq[--spos]];
		    cnf[ppos] = s->conf[spos]; /* FIXME: 1/base only */
		    pos[ppos] = ++rpos;
		    nth[ppos] = 0;
		    ppos++;
		}
		rnth = 0;
		break;

	    case 'I':
		for (j = 0; j < cigar_len; j++) {
		    seq[ppos] = complementary_base[s->seq[--spos]];
		    cnf[ppos] = s->conf[spos]; /* FIXME: 1/base only */
		    pos[ppos] = rpos;
		    nth[ppos] = ++rnth;
		    ppos++;
		}
		break;

	    case 'P':
		for (j = 0; j < cigar_len; j++) {
		    seq[ppos] = '+';
		    cnf[ppos] = 0;
		    pos[ppos] = rpos;
		    nth[ppos] = ++rnth;
		    ppos++;
		}
		break;

	    case 'D':
		for (j = 0; j < cigar_len; j++) {
		    seq[ppos] = '*';
		    cnf[ppos] = 0;
		    pos[ppos] = ++rpos;
		    nth[ppos] = 0;
		    ppos++;
		}
		rnth = 0;
		break;

	    default:
		fprintf(stderr, "Unhandled cigar opcode in seq_padded()\n");
	    }
	}
    } else {
	for (rpos=-1, spos = ppos = i = 0; i < s->alignment_len; i += 2) {
	    int cigar_op  = s->alignment[i+1];
	    int cigar_len = s->alignment[i];

	    switch (cigar_op) {
	    case 'S':
	    case 'M':
		for (j = 0; j < cigar_len; j++) {
		    seq[ppos] = s->seq[spos];
		    cnf[ppos] = s->conf[spos++];
		    pos[ppos] = ++rpos;
		    nth[ppos] = 0;
		    ppos++;
		}
		rnth = 0;
		break;

	    case 'I':
		for (j = 0; j < cigar_len; j++) {
		    seq[ppos] = s->seq[spos];
		    cnf[ppos] = s->conf[spos++];
		    pos[ppos] = rpos;
		    nth[ppos] = ++rnth;
		    ppos++;
		}
		break;

	    case 'P':
		for (j = 0; j < cigar_len; j++) {
		    seq[ppos] = '+';
		    cnf[ppos] = 0;
		    pos[ppos] = rpos;
		    nth[ppos] = ++rnth;
		    ppos++;
		}
		break;

	    case 'D':
		for (j = 0; j < cigar_len; j++) {
		    seq[ppos] = '*';
		    cnf[ppos] = 0;
		    pos[ppos] = ++rpos;
		    nth[ppos] = 0;
		    ppos++;
		}
		rnth = 0;
		break;

	    default:
		fprintf(stderr, "Unhandled cigar opcode in seq_padded()\n");
	    }
	}
    }
    
    return 0;
}

/*
 * Replaces the sequence in s with a new sequence, using 'seq' as the basis
 * to work from.
 */
int seq_depad(GapIO *io, seq_t **s, char *seq, char *conf, int *pos, int *nth,
	      int len, int comp) {
    int i, j;
    int dlen = 0;
    unsigned char *new_al, *al;
    int al_size;
    seq_t *n = *s;

    /* For use later in this function */
#define EXTEND_ALIGNMENT(op, cur_len) \
    do {                                                        \
	int l = cur_len < 255 ? cur_len : 255;			\
	if (al+2 - new_al > al_size) {				\
	    size_t diff = al-new_al;				\
	    al_size *= 2;					\
	    if (NULL == (new_al = realloc(new_al, al_size)))	\
		return -1;					\
	    al = new_al + diff;					\
	}							\
	*al++ = l; cur_len -= l;				\
	*al++ = op;                                             \
    } while (cur_len > 0);


    /* Find depadded length */
    for (i = 0; i < len; i++) {
	if (seq[i] == '*' || seq[i] == '+' || seq[i] == '-')
	    continue;
	dlen++;
    }

    /* Generate new CIGAR from seq/pos/nth arrays */
    al_size = (*s)->alignment_len + 4;
    if (NULL == (al = new_al = malloc(al_size)))
	return -1;

    for (i = 0; i < len;) {
	int cur_len;

	switch(seq[i]) {
	case '*':
	    for (cur_len = 0; i < len && seq[i] == '*'; i++)
		cur_len++;
	    printf("%dD ", cur_len);
	    EXTEND_ALIGNMENT('D', cur_len);
	    break;

	case '+':
	    for (cur_len = 0; i < len && seq[i] == '+'; i++)
		cur_len++;
	    printf("%dP ", cur_len);
	    EXTEND_ALIGNMENT('P', cur_len);
	    break;

	case '-':
	    /* Also remove any padding to the left? */
	    if (al-new_al > 2 && al[-1] == 'P') {
		al -= 2;
	    }
	    i++;
	    break;

	default:
	    if (nth[i] == 0) {
		for (cur_len = 0; i < len && nth[i] == 0 && seq[i] != '*' && seq[i] != '-'; i++)
		    cur_len++;
		printf("%dM ", cur_len);
		EXTEND_ALIGNMENT('M', cur_len);
	    } else {
		for (cur_len = 0; i < len && nth[i] != 0 && seq[i] != '+' && seq[i] != '-'; i++)
		    cur_len++;
		printf("%dI ", cur_len);
		EXTEND_ALIGNMENT('I', cur_len);
	    }
	}
    }
    putchar('\n');
    al_size = al - new_al;

    if (comp) {
	while (al > new_al) {
	    al -= 2;
	    printf("%d%c ", al[0], al[1]);
	}
	putchar('\n');
    } else {
	al = new_al;
	for (i = 0; i < al_size; i+= 2) {
	    printf("%d%c ", al[i], al[i+1]);
	}
	putchar('\n');
    }


    /* Resize struct if needed */
    if (al_size != n->alignment_len || dlen != ABS(n->len)) {
	n->alignment_len = al_size;
	n->len = n->len > 0 ? dlen : -dlen;

	n = cache_item_resize(n, sizeof(*n) + sequence_extra_len(n));
	sequence_reset_ptr(n);
    }

    /* Write seq/conf/cigar back again */
    if (comp) {
	for (i = al_size-2, j = 0; j < al_size; i-=2, j += 2) {
	    n->alignment[j  ] = new_al[i  ];
	    n->alignment[j+1] = new_al[i+1];
	}

	for (i = len-1, j = 0; i >= 0; i--) {
	    if (seq[i] == '*' || seq[i] == '+' || seq[i] == '-')
		continue;
	    n->seq[j]  = complementary_base[seq[i]];
	    n->conf[j] = conf[i];
	    j++;
	}
    } else {
	for (j = 0; j < al_size; j += 2) {
	    n->alignment[j  ] = new_al[j  ];
	    n->alignment[j+1] = new_al[j+1];
	}

	for (i = j = 0; i < len; i++) {
	    if (seq[i] == '*' || seq[i] == '+' || seq[i] == '-')
		continue;
	    n->seq[j]  = seq[i];
	    n->conf[j] = conf[i];
	    j++;
	}
    }

    free(new_al);

    return 0;
}


/* Cigar editing ops */
int seq_insert_cigar(seq_t **s, int idx, int len, char op, int merge) {
    seq_t *n = *s;
    int i, extra_len;
    char *old_seq, *new_seq;

    if (merge && idx >= 0) {
	if (n->alignment[idx+1] == op && n->alignment[idx]+len <= 255) {
	    n->alignment[idx] += len;

	    if (idx >= 2 && n->alignment[idx-2+1] == op &&
		    n->alignment[idx-2] + n->alignment[idx] <= 255) {
		/* eg 1M (+1M) 1M to 1M 2M to 3M */
		puts("Can merge CIGAR further");
	    }
	    return 0;
	} else if (idx >= 2 && n->alignment[idx-2+1] == op &&
		   n->alignment[idx-2] + len <= 255) {
	    n->alignment[idx-2] += len;
	    return 0;
	}
    }

    if (idx < 0) idx = 0;

    /* Realloc */
    extra_len = sequence_extra_len(*s) + 2;
    n = cache_item_resize(n, sizeof(*n) + extra_len);
    if (NULL == n)
	return -1;
    *s = n;

    /* Shifting */
    old_seq = n->seq;
    n->alignment_len += 2;
    sequence_reset_ptr(n);
    new_seq = n->seq;

    extra_len  = ABS(n->len);
    extra_len += (n->format == SEQ_FORMAT_CNF4 ? 4 : 1) * ABS(n->len);
    extra_len += n->aux_len;
    memmove(new_seq, old_seq, extra_len);

    /* Insert cigar op */
    for (i = n->alignment_len-4; i >= idx; i -= 2) {
	n->alignment[i+2] = n->alignment[i];
	n->alignment[i+3] = n->alignment[i+1];
    }
    n->alignment[idx] = len;
    n->alignment[idx+1] = op;

    return 0;
}

int seq_delete_cigar(seq_t **s, int idx, int merge) {
    seq_t *n = *s;
    int i = idx, j = idx+2;
    char *old_seq, *new_seq;
    int len;

    if (merge) {
	if (idx > 2 && idx < n->alignment_len-2) {
	    if (n->alignment[idx+3] == n->alignment[idx-1] &&
		n->alignment[idx+2] + n->alignment[idx-2] <= 255) {
		n->alignment[idx-2] += n->alignment[idx+2];
		j = idx+4;
	    }
	}
    }

    /* Remove the item */
    while (j < n->alignment_len) {
	n->alignment[i]   = n->alignment[j];
	n->alignment[i+1] = n->alignment[j+1];
	i += 2;
	j += 2;
    }
    n->alignment_len = i;

    /* Shuffle down memory */
    old_seq = n->seq;
    sequence_reset_ptr(n);
    new_seq = n->seq;

    len  = ABS(n->len); /* seq */
    len += (n->format == SEQ_FORMAT_CNF4 ? 4 : 1) * ABS(n->len); /* conf */
    len += n->aux_len;
    memmove(new_seq, old_seq, len);

    return 0;
}

/* Debugging aid */
static void sequence_dump_cigar(seq_t *s) {
    int i;

    printf("CIGAR:");
    for (i = 0; i < s->alignment_len; i+=2) {
	printf(" %d%c", s->alignment[i], s->alignment[i+1]);
    }
    putchar('\n');
}

//#define CIGAR_IND(s,comp,ind,off) (((comp) == 0) ? (ind) + (off) : (s)->alignment_len - ((ind) + (off)*!(comp)) - 2)


#define CIGAR_IND(s,comp,ind,off) ((ind)+(off))

/*
 * Replaces a cigar op with another. Eg changing M into I (breaking the
 * alignment of course). This allows us to do things like turn 10M into
 * 5M 1I 5M, essentially marking a base as an insertion adjusting the
 * alignment of all other bases accordingly.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_edit_op(GapIO *io, seq_t **s, int pos, int nth, int op,
		     int indel) {
    seq_t *n = *s;
    int cigar_ind, cigar_op, cigar_len, comp;
    int ret = 0;

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    sequence_dump_cigar(n);

    if (-1 == sequence_pos2cigar(io, n, pos, nth,
				 &cigar_ind, &cigar_op, &cigar_len,
				 NULL, NULL, NULL, &comp))
	return -1;

    if (comp) {
	complement_seq_t(n);
	if (-1 == sequence_pos2cigar(io, n, pos, nth,
				     &cigar_ind, &cigar_op, &cigar_len,
				     NULL, NULL, NULL, NULL)) {
	    complement_seq_t(n);
	    return -1;
	}
    }

    /*
     * Not all combinations are valid. The ones implemented are:
     *
     * We can replace 'M' with 'I' and 'I' with 'M'.
     * We can replace 'D' with 'P' and vice versa (indel=0)
     * We can insert 'D' or 'P' (indel=1)
     * We can remove 'D' or 'P' (op = '\0' && indel=1).
     *
     * Combinations like indel=1 and current op being M/I are disallowed.
     */
    if (indel == 0 &&
	((cigar_op == 'M' && op == 'I') ||
	 (cigar_op == 'I' && op == 'M') ||
	 (cigar_op == 'D' && op == 'P') ||
	 (cigar_op == 'P' && op == 'D'))) {
	/* Toggle I/M or D/P */
	if (cigar_len+1 == n->alignment[cigar_ind]) {
	    /* At start */
	    n->alignment[cigar_ind]--;
	    seq_insert_cigar(&n, cigar_ind, 1, op, 1);
	} else if (cigar_len == 0 && n->alignment[cigar_ind] > 1) {
	    /* At end */
	    n->alignment[cigar_ind]--;
	    seq_insert_cigar(&n, cigar_ind+2, 1, op, 1);
	} else if (cigar_len == 0) {
	    /* Total replacement as current op is len 1 anyway */
	    seq_delete_cigar(&n, cigar_ind, 0);
	    seq_insert_cigar(&n, cigar_ind, 1, op, 1);
	} else {
	    /* In the middle of a run, so split */
	    seq_insert_cigar(&n, cigar_ind+2, cigar_len, cigar_op, 0);
	    seq_insert_cigar(&n, cigar_ind+2, 1, op, 0);
	    n->alignment[cigar_ind] -= (cigar_len+1);
	}

    } else if (indel == 1 && (op == 'D' || op == 'P')) {
	/* Insertion of D or P operation */
	if (cigar_len+1 == n->alignment[cigar_ind]) {
	    /* At start */
	    seq_insert_cigar(&n, cigar_ind, 1, op, 1);
	} else if (cigar_len == 0 && n->alignment[cigar_ind] > 1) {
	    /* At end */
	    seq_insert_cigar(&n, cigar_ind+2, 1, op, 1);
	} else {
	    /* In the middle of a run, so split */
	    seq_insert_cigar(&n, cigar_ind+2, cigar_len+1, cigar_op, 0);
	    seq_insert_cigar(&n, cigar_ind+2, 1, op, 0);
	    n->alignment[cigar_ind] -= (cigar_len+1);
	}

    } else if (indel == 1 && op == '\0' &&
	       (cigar_op == 'D' || cigar_op == 'P')) {
	/* Remove of D or P operation */
	if (n->alignment[cigar_ind] > 1)
	    n->alignment[cigar_ind]--;
	else
	    seq_delete_cigar(&n, cigar_ind, 1);

    } else {
	fprintf(stderr, "Unhandled cigar edit in sequence_edit_op()\n");
	ret = -1;;
    }

    if (comp) {
	complement_seq_t(n);
    }

    return ret;
}


/* 
 * As sequence_replace_base but using unpadded reference coords (relative
 * to start of sequence).
 */
int sequence_replace_ubase(GapIO *io, seq_t **s, int pos, int nth,
			   char base, int conf, int contig_orient) {
    seq_t *n;
    int cigar_ind, cigar_op, cigar_len, cigar_rep;
    int spos, nth_left, comp;
    int is_pad = 0;
    int i;

    printf("Edit pos %d nth %d\n", pos, nth);

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    sequence_dump_cigar(n);

    /*
     * Alternative code method, produced by padding out the sequence,
     * applying the edit, and converting back again.
     * Only in debug form at present.
     */
    {
	char *seq;
	char *qual;
	int *rpos;
	int *rnth;
	int plen = nth;
	int i, j, gap;

	if (-1 == seq_padded(io, n, &rpos, &rnth, &seq, &qual, &plen, &comp)) {
	    fprintf(stderr, "seq_padded failed\n");
	    return -1;
	}
	
	/* Find pos, brute force */
	for (i = 0; i < plen; i++) {
	    if (rpos[i] == pos)
		break;
	}
	if (i == plen) {
	    fprintf(stderr, "ERROR: off end of sequence?\n");
	    return -1;
	}

	/* Find nth base */
	for (; i < plen; i++) {
	    if (rpos[i] != pos) {
		i--;
		break;
	    }
	    if (rnth[i] == nth)
		break;
	}

	/* Expand sequence as needed */
	gap = nth - rnth[i];
	if (gap) {
	    printf("Need to insert %d bases\n", gap);
	    memmove(&seq [i+1+gap], &seq [i+1], (plen-(i+1)) * sizeof(*seq));
	    memmove(&qual[i+1+gap], &qual[i+1], (plen-(i+1)) * sizeof(*qual));
	    memmove(&rpos[i+1+gap], &rpos[i+1], (plen-(i+1)) * sizeof(*rpos));
	    memmove(&rnth[i+1+gap], &rnth[i+1], (plen-(i+1)) * sizeof(*rnth));

	    if (i+1 < n->right)
		n->right++;
	    if (i+1 < n->left)
		n->left++;

	    for (j = 1; j <= gap; j++) {
		seq [i+j] = '+';
		qual[i+j] = 0;
		rpos[i+j] = rpos[i];
		rnth[i+j] = rnth[i]+j;
	    }
	    i += gap;
	    plen += gap;
	}

	//printf("i=%d => %d/%d %c\n", i, rpos[i], rnth[i], seq[i]);

	/* Adjust clip points */
	if (nth && (base == '*' /*|| base == '-'*/))
	    base = '+';
	if (!nth && (base == '+' || base == '-'))
	    base = '*';

	if ((base == '*' || base == '+') && seq[i] != '*' && seq[i] != '+') {
	    if (i-1 < n->right)
		n->right--;
	    if (i+1 < n->left)
		n->left--;
	}

	if (base != '*' && base != '+' && (seq[i] == '*' || seq[i] == '+')) {
	    if (i-1 < n->right)
		n->right++;
	    if (i+1 < n->left)
		n->left++;
	}

	/* Make the edit */
	seq [i] = base;
	qual[i] = conf;

//	for (j = 0; j < plen; j++) {
//	    printf("%d/%d %c%s\n",
//		   rpos[j], rnth[j], seq[j], i==j ? " <==" : "");
//	}

	/* Convert back to unpadded form */
	if (-1 == seq_depad(io, &n, seq, qual, rpos, rnth, plen, comp))
	    return -1;

	free(seq);
	free(qual);
	free(rpos);
	free(rnth);

	return 0;
    }

    /*
     * Old method follows.
     *
     * This is more efficient than the above, but vastly more complicated
     * and very tricky to cover every single corner case. For safety I'll
     * stick with the belt and braces method above and decide whether to
     * use the scarier version later if and when we find efficiency is a
     * problem.
     */
#if 0
    if (-1 == sequence_pos2cigar(io, n, pos, nth, &cigar_ind, &cigar_op,
				 &cigar_len, &spos, NULL, &nth_left, &comp))
	return -1;

    if (base == '-' || base == '*' || base == '+')
	is_pad = 1;

    if (nth_left && is_pad) /* null op */
	return 0; 

    sequence_invalidate_consensus(io, n);

    if (comp) {
	/*
	 * This is nasty, but the code to edit the cigar is already complex
	 * and making it work in both orientations is horrendous. So instead
	 * we cheat by complementing it, recomputing the positions again,
	 * make the edit as if uncomplemented data, and then complement back 
	 * again.
	 */
	complement_seq_t(n);

	if (-1 == sequence_pos2cigar(io, n, pos, nth, &cigar_ind, &cigar_op,
				     &cigar_len, &spos, NULL, &nth_left, NULL))
	    return -1;
    }


    printf("Seq=%.*s\n", ABS(n->len), n->seq);
    printf("cigar_ind=%d\n", cigar_ind);
    printf("cigar_len=%d/%d\n", cigar_len, n->alignment[cigar_ind]);
    printf("cigar_op=%d\n", cigar_op);

    if (!is_pad && nth_left == 0 && cigar_op != 'D' && cigar_op != 'P') {
	/* Just a replacement */
	puts("sub");

	n->seq[spos] = base;
	n->conf[spos] = conf;

    } else if (is_pad && (cigar_op == 'D' || cigar_op == 'P')) {
	puts("nop");
	/* Null op */

    } else if (is_pad && nth_left == 0) {
	/* Replacing existing base with a pad => deletion */
	puts("del");
	sequence_delete_base(io, &n, spos, -1);

	switch (cigar_op) {
	case 'I':
	case 'M':
	    /* Replace M or I by D or P */
	    cigar_rep = cigar_op == 'M' ? 'D' : 'P';
	    if (cigar_len == 0) {
		/* At end of [MI] */
		if (n->alignment[cigar_ind] > 1) {
		    n->alignment[cigar_ind]--;
		    seq_insert_cigar(&n, cigar_ind+2, 1, cigar_rep, 1);
		} else {
		    /* Swap 1[MI] for 1P */
		    seq_delete_cigar(&n, cigar_ind, 0);
		    seq_insert_cigar(&n, cigar_ind, 1, cigar_rep, 1);
		}
	    } else if (cigar_len+1 == n->alignment[cigar_ind]) {
		/* At start of [MI] */
		n->alignment[cigar_ind]--;
		seq_insert_cigar(&n, cigar_ind, 1, cigar_rep, 1);
	    } else {
		/* Split [MI] in middle */
		seq_insert_cigar(&n, cigar_ind+2, cigar_len, cigar_op, 0);
		seq_insert_cigar(&n, cigar_ind+2, 1, cigar_rep, 0);
		assert(n->alignment[cigar_ind] > 1);
		n->alignment[cigar_ind] -= (cigar_len+1);
	    }
	    break;

	case 'S':
	case 'N':
	    puts("HELP - disallowed?!");
	    break;
	}

    } else {
	/* nth_left != 0 => Replacing a gap with a base => insertion */
	puts("ins");
	spos++;
	sequence_insert_base(io, &n, spos, base, conf, -1);

	if (n->alignment[cigar_ind] == cigar_len) {
	    /*
	     * Edit is between this location and previous base - so in an
	     * insertion caused by another sequence.
	     */
	    cigar_ind -= 2;
	    cigar_len = 0;
	    cigar_op  = n->alignment[cigar_ind+1];
	}

	switch (cigar_op) {
	case 'M':
	    cigar_rep = 'I';
	    if (nth) {
		/* Split */
		if (cigar_len) {
		    n->alignment[cigar_ind] -= cigar_len;
		    seq_insert_cigar(&n, cigar_ind+2, cigar_len, cigar_op, 0);
		}
		seq_insert_cigar(&n, cigar_ind+2, 1, cigar_rep, 0);
		if (nth_left > 1)
		    seq_insert_cigar(&n, cigar_ind+2, nth_left-1, 'P', 0);
	    } else {
		puts("nth == 0: How is this possible?");
	    }
	    break;

	case 'I':
	    cigar_rep = 'I';
	    /* Assume appending? */
	    if (cigar_len == 0) {
		seq_insert_cigar(&n, cigar_ind+2, 1, cigar_rep, 0);
		if (nth_left > 1)
		    seq_insert_cigar(&n, cigar_ind+2, nth_left-1, 'P', 0);
	    } else {
		puts("cigar_len != 0: How is this possible?");
	    }
	    break;

	case 'D':
	case 'P':
	    /* Both D and P can be inserted into the middle of a run */
	    cigar_rep = cigar_op == 'D' ? 'M' : 'I';
	    if (cigar_len == 0) {
		/* End of a run */
		if (nth_left >= 1) {
		    seq_insert_cigar(&n, cigar_ind+2, 1,
				     nth ? 'I' : 'M', 0);
		    if (nth_left > 1)
			seq_insert_cigar(&n, cigar_ind+2, nth_left-1, 'P', 1);
		} else {
		    if (n->alignment[cigar_ind] > 1) {
			n->alignment[cigar_ind]--;
			seq_insert_cigar(&n, cigar_ind+2, 1, cigar_rep, 1);
		    } else {
			/* The entire D/P run */
			seq_delete_cigar(&n, cigar_ind, 0);
			seq_insert_cigar(&n, cigar_ind, 1, cigar_rep, 1);
		    }
		}
	    } else if (nth_left == 0 && cigar_len+1 == n->alignment[cigar_ind]) {
		/* Start of a run */
		n->alignment[cigar_ind]--;
		seq_insert_cigar(&n, cigar_ind, 1, cigar_rep, 1);
	    } else {
		/* Split a run in the middle */
		seq_insert_cigar(&n, cigar_ind+2, cigar_len, cigar_op, 0);
		if (nth) {
		    seq_insert_cigar(&n, cigar_ind+2, 1, 'I', 0);
		    if (nth_left > 1)
			seq_insert_cigar(&n, cigar_ind+2, nth_left-1, 'P', 0);
		} else {
		    seq_insert_cigar(&n, cigar_ind+2, 1, 'M', 0);
		}
		assert(n->alignment[cigar_ind] > 1);
		n->alignment[cigar_ind] -= cigar_len + (nth_left == 0);
	    }
	    break;

	default:
	    printf("Unsupport cigar op '%d' in insertion\n",
		   cigar_op);
	}
    }

    if (comp) {
	complement_seq_t(n);
    }

    sequence_dump_cigar(n);
    *s = n;

    return 0;
#endif
}


/*
 * Insert into position 'pos' of a sequence, where pos starts from 0
 * (before first base). We insert just after the position as this allows
 * us to determine the number of pads to incorporate position pos as the
 * new inserted base.
 */
int sequence_insert_ubase(GapIO *io, seq_t **s, int pos, int nth,
			  char base, char conf, int contig_orient) {
    
    seq_t *n;
    int cigar_ind, cigar_op, cigar_len, cigar_rep;
    int spos, nth_left, comp;

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    sequence_invalidate_consensus(io, n);

    sequence_dump_cigar(n);

    /* Pad & depad method */
    {
	char *seq;
	char *qual;
	int *rpos;
	int *rnth;
	int plen = nth+1;
	int i, j, gap;

	if (-1 == seq_padded(io, n, &rpos, &rnth, &seq, &qual, &plen, &comp)) {
	    fprintf(stderr, "seq_padded failed\n");
	    return -1;
	}
	
	/* Find pos, brute force */
	for (i = 0; i < plen; i++) {
	    if (rpos[i] == pos)
		break;
	}
	if (i == plen) {
	    fprintf(stderr, "ERROR: off end of sequence?\n");
	    return -1;
	}

	/* Find nth base */
	for (; i < plen; i++) {
	    if (rpos[i] != pos) {
		i--;
		break;
	    }
	    if (rnth[i] == nth)
		break;
	}

	/* Expand sequence as needed */
	gap = nth - rnth[i];
	if (!gap) gap=1;
	if (gap) {
	    printf("Need to insert %d bases\n", gap);
	    memmove(&seq [i+1+gap], &seq [i+1], (plen-(i+1)) * sizeof(*seq));
	    memmove(&qual[i+1+gap], &qual[i+1], (plen-(i+1)) * sizeof(*qual));
	    memmove(&rpos[i+1+gap], &rpos[i+1], (plen-(i+1)) * sizeof(*rpos));
	    memmove(&rnth[i+1+gap], &rnth[i+1], (plen-(i+1)) * sizeof(*rnth));

	    if (i < n->right)
		n->right++;
	    if (i < n->left)
		n->left++;

	    for (j = 1; j <= gap; j++) {
		seq [i+j] = '+';
		qual[i+j] = 0;
		rpos[i+j] = rpos[i];
		rnth[i+j] = rnth[i]+j;
	    }
	    i += gap;
	    plen += gap;
	}
	
	/* Make the edit */
	if (base == '*' || base == '-')
	    base = '+';
	seq [i] = base;
	qual[i] = conf;

//	for (j = 0; j < plen; j++) {
//	    printf("%d/%d %c%s\n",
//		   rpos[j], rnth[j], seq[j], i==j ? " <==" : "");
//	}

	/* Convert back to unpadded form */
	if (-1 == seq_depad(io, &n, seq, qual, rpos, rnth, plen, comp))
	    return -1;

	free(seq);
	free(qual);
	free(rpos);
	free(rnth);

	return 0;
    }

    /* Old method - see sequence_replace_ubase for discussion */
#if 0
    if (-1 == sequence_pos2cigar(io, n, pos, nth, &cigar_ind, &cigar_op,
				 &cigar_len, &spos, NULL, &nth_left, &comp))
	return -1;


    if (n->alignment[cigar_ind] == cigar_len) {
	/*
	 * Edit is between this location and previous base - so in an
	 * insertion caused by another sequence.
	 */
	cigar_ind -= 2;
	cigar_len = 0;
	cigar_op  = n->alignment[cigar_ind+1];
    }

    if (base == '-' || base == '+' || base == '*') {
	cigar_rep = 'P';
    } else {
	sequence_insert_base(io, &n, spos + (nth > 0), base, conf, 0);
	cigar_rep = 'I';
    }

    /*
     * All edits are insertions, which makes this easier than replace.
     * However the occur to the left of the base we've identified.
     */
    switch (cigar_op) {
    case 'M':
    case 'I':
    case 'P':
    case 'D':
	if (cigar_len == 0 && nth > 0) {
	    /* End of a run */
	    seq_insert_cigar(&n, cigar_ind+2, 1, cigar_rep, 1);
	    if (nth_left > 1)
		seq_insert_cigar(&n, cigar_ind+2, nth_left-1, 'P', 0);
	} else if (cigar_len+1 == n->alignment[cigar_ind] && nth_left == 0) {
	    /* Start of a run */
	    /* Can't tell how many nth are left at present */
	    seq_insert_cigar(&n, cigar_ind, 1, cigar_rep, 1);
	} else {
	    /* In middle */
	    if (nth) {
		seq_insert_cigar(&n, cigar_ind+2, cigar_len, cigar_op, 0);
		seq_insert_cigar(&n, cigar_ind+2, 1, cigar_rep, 0);
		if (nth_left > 1)
		    seq_insert_cigar(&n, cigar_ind+2, nth_left-1, 'P', 0);
		n->alignment[cigar_ind] -= cigar_len;
	    } else {
		seq_insert_cigar(&n, cigar_ind+2, cigar_len+1, cigar_op, 0);
		seq_insert_cigar(&n, cigar_ind+2, 1, cigar_rep, 0);
		n->alignment[cigar_ind] -= cigar_len+1;
	    }
	}
	break;
    }

    sequence_dump_cigar(n);
    *s = n;

    return 0;
#endif
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

    if (contig_orient > 0) {
	pos = sequence_orient_pos(io, &n, pos, &comp);
	if (comp)
	    pos++;
    } else if (contig_orient == 0) {
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
	    ABS(n->len)-1 - pos +
	    sequence_conf_size(n) * (ABS(n->len)-1));

    c_old++;
    memmove(&c_old[sequence_conf_size(n)*pos]+1,
	    &n->conf[sequence_conf_size(n)*pos],
	    (ABS(n->len)-1 - pos) * sequence_conf_size(n));
	    
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
int sequence_delete_ubase(GapIO *io, seq_t **s, int pos, int nth,
			  int contig_orient) {
    seq_t *n;
    int cigar_ind, cigar_op, cigar_len, cigar_rep;
    int spos, nth_left, comp;

    if (pos < 0)
	return 0;

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    sequence_invalidate_consensus(io, n);

    sequence_dump_cigar(n);

    /* Pad & depad method */
    {
	char *seq;
	char *qual;
	int *rpos;
	int *rnth;
	int plen = nth;
	int i, j;

	if (-1 == seq_padded(io, n, &rpos, &rnth, &seq, &qual, &plen, &comp)) {
	    fprintf(stderr, "seq_padded failed\n");
	    return -1;
	}
	
	/* Find pos, brute force */
	for (i = 0; i < plen; i++) {
	    if (rpos[i] == pos)
		break;
	}
	if (i == plen) {
	    fprintf(stderr, "ERROR: off end of sequence?\n");
	    return -1;
	}

	/* Find nth base */
	for (; i < plen; i++) {
	    if (rpos[i] != pos) {
		i--;
		break;
	    }
	    if (rnth[i] == nth)
		break;
	}

//	for (j = 0; j < plen; j++) {
//	    printf("%d/%d %c%s\n",
//		   rpos[j], rnth[j], seq[j], i==j ? " <==" : "");
//	}


	/* Expand sequence as needed */
	if (nth - rnth[i])
	    goto nowt_to_do;

	if (seq[i] == '+' || nth) {
	    memmove(&seq [i], &seq [i+1], (plen-(i+1)) * sizeof(*seq));
	    memmove(&qual[i], &qual[i+1], (plen-(i+1)) * sizeof(*qual));
	    memmove(&rpos[i], &rpos[i+1], (plen-(i+1)) * sizeof(*rpos));
	    memmove(&rnth[i], &rnth[i+1], (plen-(i+1)) * sizeof(*rnth));
	    plen--;
	} else {
	    seq[i] = '*';
	}

//	puts("");
//	for (j = 0; j < plen; j++) {
//	    printf("%d/%d %c%s\n",
//		   rpos[j], rnth[j], seq[j], i==j ? " <==" : "");
//	}

	/* Convert back to unpadded form */
	if (-1 == seq_depad(io, &n, seq, qual, rpos, rnth, plen, comp))
	    return -1;

    nowt_to_do:
	free(seq);
	free(qual);
	free(rpos);
	free(rnth);

	return 0;
    }

    /* Old method - see sequence_replace_ubase for discussion */
#if 0
    if (-1 == sequence_pos2cigar(io, n, pos, nth, &cigar_ind, &cigar_op,
				 &cigar_len, &spos, NULL, &nth_left, &comp))
	return -1;


    printf("cigar_op  = %c\n", cigar_op);
    printf("cigar_len = %d/%d\n", cigar_len, n->alignment[cigar_ind]);
    printf("cigar_ind = %d\n", cigar_ind);
    printf("nth = %d, left = %d\n", nth, nth_left);
    printf("spos = %d\n", spos);

    if (nth_left)
	/* Nothing at this position to delete */
	return -1;

    switch(cigar_op) {
    case 'I':
    case 'M':
	/* Removing I or D also implies removing the base associated with it */
	sequence_delete_base(io, &n, spos, 0);
	
	/* NOTE: fall through intended */

    case 'D':
    case 'P':
	/* Deleting a D or P is easy enough - just remove it */
	if (n->alignment[cigar_ind] > 1)
	    n->alignment[cigar_ind]--;
	else
	    seq_delete_cigar(&n, cigar_ind, 1);
	break;
    }

    sequence_dump_cigar(n);
    *s = n;

    return 0;
#endif
}

/*
 * Delete position 'pos' of a sequence, where pos starts from 0
 */
int sequence_delete_base(GapIO *io, seq_t **s, int pos, int contig_orient) {
    seq_t *n;
    int comp = 0;
    size_t extra_len = sequence_extra_len(*s);

    if (pos >= ABS((*s)->len) || pos < 0)
	return 0;

    if (!(n = cache_rw(io, *s)))
	return -1;
    *s = n;

    sequence_invalidate_consensus(io, n);

    if (contig_orient > 0) {
	pos = sequence_orient_pos(io, &n, pos, &comp);
    } else if (contig_orient == 0) {
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
	    (ABS(n->len) - pos) * sequence_conf_size(n));

    memmove(&n->seq[pos], &n->seq[pos+1],
	    (ABS(n->len) - pos) + pos * sequence_conf_size(n));

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
    int start, end, orient;
    tg_rec contig;
    rangec_t *r;
    int nr, i;

    /* Find the position and contig this sequence currently covers */
    cache_incr(io, *s);
    if (0 != sequence_get_position(io, (*s)->rec,
				   &contig, &start, &end, &orient)) {
	cache_decr(io, *s);
	return -1;
    }


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

    return 0;
}
