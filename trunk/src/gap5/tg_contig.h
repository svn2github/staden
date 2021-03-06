#ifndef _TG_CONTIG_H_
#define _TG_CONTIG_H_

#include <limits.h>
#include "tree.h"

/*
 * Define this to make the cached consensus visible within the views,
 * by allocating it a Y coordinate.  This aids debugging of the
 * caching algorithm, but should be disabled in production versions to
 * avoid confusing the user.
 */
/* #define CACHED_CONS_VISIBLE */

/*
 * 'get' functions - simply returns the structure member.
 *
 * <type> contig_get_XXX(contig_t **c)
 */
#define contig_get_start(c)  ((*(c))->start)
#define contig_get_end(c)    ((*(c))->end)
#define contig_get_bin(c)    ((*(c))->bin)
#define contig_get_name(c)   ((*(c))->name)
#define contig_get_length(c) (contig_get_end(c) - contig_get_start(c) + 1)


/*
 * 'set' functions all have prototype:
 *
 * int contig_set_XXX(GapIO *io, contig_t **c, <type> new_value)
 *
 * Returns 0 for success, possibly also modifying *c pointer
 *        -1 for failure
 */
int contig_set_start(GapIO *io, contig_t **c, int value);
int contig_set_end  (GapIO *io, contig_t **c, int value);
int contig_set_bin  (GapIO *io, contig_t **c, tg_rec value);
int contig_set_name (GapIO *io, contig_t **c, char *name);

/*
 * Returns the offset to apply to the root bin of the contig.
 * This takes into account whether it have been complemented or not.
 */
int contig_offset(GapIO *io, contig_t **c);

int contig_insert_base(GapIO *io, contig_t **c, int pos, char base, int conf);

int contig_insert_bases(GapIO *io, contig_t **c, int pos, char base, int conf,
			int nbases);

typedef struct {
    tg_rec rec;
    int    pos;
    char   base;
    int8_t conf;
} col_inserted_base;

int contig_insert_column(GapIO *io, contig_t **c, int pos,
			 size_t count, col_inserted_base *bases);

int contig_delete_base(GapIO *io, contig_t **c, int pos);

int contig_delete_pad(GapIO *io, contig_t **c, int pos);

int contig_shift_base(GapIO *io, contig_t **c, int pos, int dir);


contig_t *find_contig_by_name(GapIO *io, char *name);
tg_rec contig_index_query(GapIO *io, char *name);
tg_rec contig_index_query_prefix(GapIO *io, char *prefix);
int contig_index_update(GapIO *io, char *name, int name_len, tg_rec rec);

contig_t *contig_new(GapIO *io, char *name);

/* FIXME: move this elsewhere */
#define get_bin(io, bnum) ((bin_index_t *)cache_search((io), GT_Bin, (bnum)))

rangec_t *contig_items_in_range(GapIO *io, contig_t **c, seq_sort_t *settings, int start, int end,
				int first_order, int second_order, int *count);
rangec_t *contig_seqs_in_range(GapIO *io, contig_t **c, int start, int end, int job, int *count);
rangec_t *contig_bins_in_range(GapIO *io, contig_t **c, int start, int end,
			       int job, int min_size, int *count);
rangec_t *contig_anno_in_range(GapIO *io, contig_t **c, int start, int end,
			       int job, int *count);
rangec_t *contig_cons_in_range(GapIO *io, contig_t **c, int start, int end,
			       int job, int *count);
rangec_t *contig_refpos_in_range(GapIO *io, contig_t **c, int start, int end,
				 int job, int *count);
			       
void contig_set_default_sort(seq_sort_t *set, int primary, int secondary);
void contig_set_base_sort_point(int pos);

#define CSIR_PAIR                  (1<<0)
#define CSIR_ALLOCATE_Y_SINGLE     (1<<1)
#define CSIR_ALLOCATE_Y_MULTIPLE   (1<<2)
#define CSIR_ALLOCATE_Y \
    (CSIR_ALLOCATE_Y_SINGLE | CSIR_ALLOCATE_Y_MULTIPLE)

/* Sorry this is a mess - bit fields for things that aren't compatible */
#define CSIR_SORT_BY_X             (1<<3)
#define CSIR_SORT_BY_XEND          (1<<8)
#define CSIR_SORT_BY_Y             (1<<4)
#define CSIR_SORT_BY_SEQ_TECH      (1<<7)
#define CSIR_SORT_BY_CLIPPED       (1<<9) /* bit to combine with X vs XEND */

#define CSIR_COUNT_ONLY            (1<<5)
#define CSIR_LEAVES_ONLY           (1<<6)
#define CSIR_DEFAULT               (1<<10)
#define CSIR_SORT_BY_TEMPLATE      (1<<11)
#define CSIR_SORT_BY_STRAND        (1<<12)
#define CSIR_SORT_BY_BASE          (1<<13)
#define CSIR_SORT_BY_SEQUENCE      (1<<14)
#define CSIR_SORT_BY_TEMPLATE_STATUS (1<<15)
#define CSIR_SORT_BY_LIBRARY       (1<<16)

/* ---------------------------------------------------------------------- */

/*
 * The iterator is basically a way to easily loop through all sequences
 * in a given range. It can optionally fetches the next range of data too
 * if you iterate off the edge.
 */
typedef struct {
    rangec_t *r;     /* r[] array itself and size of it */
    int nitems;
    int index;       /* current index into r[] array */
    tg_rec cnum;     /* contig number */
    int start;       /* current region fetch coords */
    int end;
    int cstart;      /* requested limits of region to fetch */
    int cend;
    int auto_extend; /* whether to extend past cstart..cend */
    int first_r;     /* True if r[] is our first search */
    int type;
    int sort_mode;   /* either CSIR_SORT_BY_X or CSIR_SORT_BY_XEND */
} contig_iterator;

#define CITER_FIRST   (0<<0)
#define CITER_LAST    (1<<0)
#define CITER_FL_MASK (1<<0)

#define CITER_ISTART        (0<<1)
#define CITER_IEND          (1<<1)
#define CITER_ICLIPPEDSTART (2<<1)
#define CITER_ICLIPPEDEND   (3<<1)
#define CITER_SE_MASK       (3<<1)

#define CITER_SMALL_BS (1<<3)
#define CITER_PAIR     (1<<4)

#define CITER_CSTART INT_MIN
#define CITER_CEND   INT_MAX



/*
 * Y position allocation functions/data structs
 */
typedef struct xy_pair {
    SPLAY_ENTRY(xy_pair) x_link; /* Sorted on X */
    SPLAY_ENTRY(xy_pair) y_link; /* Sorted on Y */
    int x;
    int y;
} xy_pair;
int x_cmp(struct xy_pair *y1, struct xy_pair *y2);
int y_cmp(struct xy_pair *y1, struct xy_pair *y2);


/*
 * Allocates and initialises a new contig_iterator struct.
 *
 * 'whence' may be either CITER_FIRST or CITER_LAST and it controls whether
 * we start the iteration point from the beginning or end of the list.
 *
 * The start and end parameters dictate the initial region to query. We
 * may specify them as either coordinates or use CITER_CSTART and CITER_CEND
 * as synonyms for the first and last coordinate in the contig.
 *
 * Whence with CITER_PAIR bit set will use the CSIR_PAIR option when
 * requesting blocks of sequences. This won't necessarily set all pairings
 * correctly (only those within a block), so be sure to call
 * sequence_get_range_pair_position() if you need updated values. This will
 * be a nop for most cases the CSIR_PAIR will have updated the timestamp for
 * reads fetched within the same block.
 *
 * 'type' can be either GRANGE_FLAG_ISSEQ or GRANGE_FLAG_ISANNO to only
 * iterate around data of that specific type, or GRANGE_FLAG_ISANY to
 * iterate around all data.
 *
 * Finally auto_extend controls whether the start..end range is just a
 * location to start iterating from (auto_extend == 1) or a hard limit
 * with no desire to iterate past that range (auto_extend == 0).
 *
 * Returns contig_iterator pointer on success
 *         NULL on failure
 */
contig_iterator *contig_iter_new_by_type(GapIO *io, tg_rec cnum,
					 int auto_extend, int whence,
					 int start, int end, int type);
contig_iterator *contig_iter_new(GapIO *io, tg_rec cnum, int auto_extend,
				 int whence, int start, int end);

/*
 * Dealloctes memory used by a contig_iterator including the cached ranges.
 */
void contig_iter_del(contig_iterator *ci);

/*
 * Returns seq range_t struct pointer on success
 *        NULL on failure (typically fallen off the end of the contig)
 */
rangec_t *contig_iter_next(GapIO *io, contig_iterator *ci);

/*
 * Returns seq range_t struct pointer on success
 *        NULL on failure (typically fallen off the end of the contig)
 */
rangec_t *contig_iter_prev(GapIO *io, contig_iterator *ci);

/*
 * Given an iterator from start..end we'll find sequences that may cover
 * e_start..e_end where e_start and e_end may possibly be larger than start
 * to end. Pictorially:
 *
 *                start      end
 *A-------------- |          |
 *B      -----    |          |
 *C   ---------------------  |
 *D   |         -----------------------------
 *E   |                        --------     |
 *F   |                            ---------------------
 *    |                                     |
 *    e_start                               e_end
 *
 * The original query can fetch back seqs C & D, but annotations entirely
 * outside this could be missed if we're doing GRANGE_FLAG_ISANY queries.
 *
 * So we expand the range to e_start to e_end to pick up annotations.
 * NOTE: the caller then need to manually filter the start/end range
 * itself to avoid then picking up seqs B and E which aren't in the
 * original range.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int iterator_expand_range(GapIO *io, tg_rec crec, int start, int end,
			  int *e_start, int *e_end);

/*
 * Queries and/or creates a track of a given type from a region of the
 * contig at a require resolution (bpv = bases per value).
 *
 * Returns track_t to free with free_track on success
 *         NULL on failure
 */
track_t *contig_get_track(GapIO *io, contig_t **c, int start, int end,
			  int type, double bpv);

/*
 * Produces a postscript file containing a plot of the contig bin structure.
 */
int contig_dump_ps(GapIO *io, contig_t **c, char *fn,
		    int depth_first, int draw_seqs);

/*
 * Produces a graphviz file containing the contig bin structure.
 */
int contig_dump_graph(GapIO *io, contig_t **c, char *fn);

/*
 * Destroys contig 'rec'.
 * Returns 0 for success
 *        -1 for failure
 */
int contig_destroy(GapIO *io, tg_rec rec);

/* temp */
int plot_seqs_in_range(GapIO *io, contig_t **c, int start, int end, void *template, void *drawing, 
			void (*plot_func)(void *, void *, void *));

/* -------------------------------------------------------------------------
 * Padded / reference coordinate mappings.
 */
int padded_to_reference_pos(GapIO *io, tg_rec cnum, int ppos, int *dir_p,
			    int *ref_id);

/*
 * Looks for a refpos marker at position ppos. If it finds one it returns
 * 0 and sets the bin and bin_idx fields.
 *
 * If it doesn't, it returns -1;
 */
int find_refpos_marker(GapIO *io, tg_rec cnum, int ppos,
		       tg_rec *bin_r, int *bin_idx_r, rangec_t *rp);

/*
 * Remove refpos marker if present.
 * 
 * io   is the GapIO struct for the database.
 * crec is the contig record number
 * pos  is the padded position in the contig
 *
 * Returns  0 if marker removed or no marker found
 *         -1 on failure 
 */
int delete_refpos_marker(GapIO *io, tg_rec crec, int pos);

/*
 * Set a refpos marker.  Will alter an existing one, or create a new one
 * as necessary.
 *
 * io   is the GapIO struct for the database.
 * c    is the contig_t struct ** for the contig
 * pos  is the padded position on the contig
 * type is the type of refpos (GRANGE_FLAG_REFPOS_INS or GRANGE_FLAG_REFPOS_DEL)
 * dir  is the direction (GRANGE_FLAG_REFPOS_FWD or GRANGE_FLAG_REFPOS_REV)
 * id   is the reference ID
 * rpos is the reference position
 * len  is the number of deleted bases (GRANGE_FLAG_REFPOS_DEL only)
 *
 * Returns  0 on success
 *         -1 on failure
 */

int set_refpos_marker(GapIO *io, contig_t **c, int pos,
		      int type, int dir, int id, int rpos, int len);


/*
 * Given a contig record and a reference position, attempt to return
 * the padded coordinate. Note this may not exist, it may in extreme cases
 * exist multiple times (after breaking and rejoining), or it may exist
 * only once but be too hard to discover. If we don't care which specific
 * reference ID to search, pass in ref_id == -1.
 *
 * We only attempt to tackle easy cases of a single reference in this contig
 * with monotonically increasing or decreasing values. (We need more data
 * stored to allow arbitrary queries to be fast.)
 *
 * Returns 0 on success, position stored in *padded_pos.
 *        -1 on failure, *padded_pos undefined.
 */
int reference_to_padded_pos(GapIO *io, tg_rec cnum, int ref_id, int ref_pos,
			    int *padded_pos);

/*
 * As above, but starting from a single known point on that reference.
 * This allows for reference positions to occur more than once.
 */
int reference_to_padded_pos2(GapIO *io, tg_rec cnum, int ref_id, int ref_pos,
			     int cur_padded_pos, int *padded_pos);

/*
 * Converts a range of padded coordinates to reference coordinates.
 * Optionally also fills out a reference seq ID array too, if non-NULL.
 *
 * ref_pos and ref_id should be allocated by the caller to be of
 * appropriate size (paddeed_end - padded_start + 1).
 * Insertions get ref_id of -1 (if non NULL) and ref_pos[] element of INT_MIN.
 *
 * If non-NULL start_pos is the first reference coordinate used. Note that the
 * read may start in an insertion, in which case nP (if non NULL) is the
 * number of preceeding padding characters before the first base in order
 * to keep the multiple sequence alignment.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int padded_to_reference_array(GapIO *io, tg_rec cnum,
			      int padded_start, int padded_end,
			      int *ref_pos, int *ref_id,
			      int *start_pos, int *nP);


/* 
 * Moves an entire contig by a relative amount left (-ve) or right (+ve).
 * Returns 0 on success
 *        -1 on failure
 */
int move_contig(GapIO *io, tg_rec crec, int distance);

/*
 * Sets the visible start of a contig, both in the contig structure and
 * also the root bin to keep everything internally consistent.
 *
 * Ie it's a more end-user version of contig_set_start.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_set_visible_start(GapIO *io, tg_rec contig, int pos);

/*
 * Copies the nseq, nanno and nrefpos from the root bin to the contig
 * struct.
 *
 * This can be necessary during algorithms that move data around, for example
 * break_contig().
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_fix_nseq(GapIO *io, contig_t *c);

/*
 * Adds a bi-directional contig link.
 * The data to link from/to is in the passed in abs_link. Coordinates here
 * are all contig absolute coords.
 *
 * Internally this link gets created at both ends (rec1, rec2) and the
 * positions are converted into relative coordinates. (The end1 and end2
 * fields in abs_link are ignored and computed as required.)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_add_link(GapIO *io, contig_link_t *abs_link);

/*
 * Converts a specific link number from relative coordinates to absolute
 * contig coordinates.
 *
 * Input is rel_link, output is abs_link.
 * Returns 0 on success
 *        -1 on failure
 */
int contig_get_link_positions(GapIO *io,
			      contig_link_t *rel_link,
			      contig_link_t *abs_link);

#endif /* _TG_CONTIG_H_ */
