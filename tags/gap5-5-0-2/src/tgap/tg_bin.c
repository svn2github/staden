#include <assert.h>
#include <string.h>
#include <math.h>

#include "xalloc.h"
#include "tg_gio.h"
#include "tg_tracks.h"

#define get_bin(io, bnum) ((bin_index_t *)cache_search((io), GT_Bin, (bnum)))

/*
 * Allocates a new bin record.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
int bin_new(GapIO *io, int pos, int sz, int parent, int parent_type) {
    int rec;
    bin_index_t bin;

    /* Initialise disk struct */
    bin.pos         = pos;
    bin.size        = sz;
    bin.start_used  = 0;
    bin.end_used    = 0;
    bin.parent      = parent;
    bin.parent_type = parent_type;
    bin.child[0]    = 0;
    bin.child[1]    = 0;
    bin.bin_id      = 0;
    bin.rng	    = NULL;
    bin.rng_rec     = 0;
    bin.flags       = BIN_BIN_UPDATED;
    bin.track       = NULL;
    bin.track_rec   = 0;

    if (-1 == (rec = io->iface->bin.create(io->dbh, &bin)))
	return -1;

    return rec;
}

/*
 * Doubles up the number of bins by adding a new root node and duplicating
 * the tree.
 *
 * It takes the old root_id as an argument and returns the new one.
 * Returns -1 on failure.
 */
static bin_index_t *contig_extend_bins_right(GapIO *io, contig_t **c) {
    int old_root_id = contig_get_bin(c);
    bin_index_t *oroot = get_bin(io, old_root_id), *nroot;
    int root_id;
    int sz = oroot->size;

    cache_incr(io, oroot);
    if (!(oroot = cache_rw(io, oroot)))
	return NULL;

    /* Create a new root */
    root_id = bin_new(io, oroot->pos, sz*2, oroot->parent, oroot->parent_type);
    nroot = get_bin(io, root_id);
    cache_incr(io, nroot);
    if (!(nroot = cache_rw(io, nroot)))
	return NULL;

    nroot->child[0] = old_root_id;
    nroot->child[1] = 0;

    nroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, nroot);

    /* Move old left bin */
    oroot->parent = root_id;
    oroot->parent_type = GT_Bin;
    oroot->pos = 0;

    oroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, oroot);

    contig_set_bin(io, c, root_id);

    return nroot;
}

static bin_index_t *contig_extend_bins_left(GapIO *io, contig_t **c) {
    int old_root_id = contig_get_bin(c);
    bin_index_t *oroot = get_bin(io, old_root_id), *nroot;
    int root_id;
    int sz = oroot->size;

    cache_incr(io, oroot);
    if (!(oroot = cache_rw(io, oroot)))
	return NULL;

    /* Create a new root */
    root_id = bin_new(io, oroot->pos-sz, sz*2, oroot->parent, oroot->parent_type);
    nroot = get_bin(io, root_id);
    cache_incr(io, nroot);
    if (!(nroot = cache_rw(io, nroot)))
	return NULL;

    nroot->child[0] = 0;
    nroot->child[1] = old_root_id;

    nroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, nroot);

    /* Move old right bin */
    oroot->parent = root_id;
    oroot->parent_type = GT_Bin;
    oroot->pos = sz;

    oroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, oroot);

    contig_set_bin(io, c, root_id);

    return nroot;
}

/*
 * Allocates and returns next range from a range array
 * Returns -1 for error.
 */
static int next_range(Array ra) {
    if (NULL == ArrayRef(ra, ArrayMax(ra)))
	return -1;
    
    return ArrayMax(ra)-1;
}

/*
 * This finds a bin suitable for a given range. If such a bin doesn't
 * exist it can optionally create one if the 'extend' flag is true.
 *
 * When a bin is found the absolute offset of that bin is returned
 * in 'offset_r' (may be NULL).
 *
 * Returns the bin pointer on success
 *         NULL on failure
 */
bin_index_t *bin_for_range(GapIO *io, contig_t **c,
			   int start, int end, int extend,
			   int *offset_r) {
    int offset;
    bin_index_t *bin = get_bin(io, contig_get_bin(c));

    static int last_c = 0;
    static GapIO *last_io = NULL;
    static bin_index_t *last_bin = NULL;
    static int last_start = 0, last_end = 0;
    static int last_offset;

    //cache_incr(io, bin);

    /*
     * If we're trying to insert a new node beyond the total bounds of the
     * root node then we need to extend the bin structure first to either
     * the left and/or the right.
     */
    while (end >= bin->pos + bin->size) {
	//cache_decr(io, bin);
	if (extend)
	    bin = contig_extend_bins_right(io, c);
	else
	    return NULL;
	//cache_incr(io, bin);
    }

    while (start < bin->pos) {
	//cache_decr(io, bin);
	if (extend)
	    bin = contig_extend_bins_left(io, c);
	else
	    return NULL;
	//cache_incr(io, bin);
    }

    /*
     * In theory we can jump straight to a candidate starting bin, possibly
     * even returning it right here if it's the min bin size, saving about
     * 10% of our CPU time in this function
     */
#if 0
    if (last_bin && last_c == (*c)->rec && last_io == io) {
	if (start >= last_start && end <= last_end) {
	    if (last_bin && last_bin->size == MIN_BIN_SIZE) {
		/* leaf node, so we can return right now */
		if (offset_r)
		    *offset_r = last_offset;
		return last_bin;
	    }

	    /* Maybe a smaller bin, but start the search from here on */
	    bin = last_bin;
	    offset = last_offset;
	    cache_incr(io, bin);
	    goto jump;
	}
    }
#endif

    /* Now recurse down the bin hierachy searching for the smallest bin */
    offset = bin->pos;
 jump:
    for (;;) {
	int i;
	bin_index_t *ch;

	/* Find which child bin is most suitable */
	cache_incr(io, bin);
	for (i = 0; i < 2;) {
	    if (bin->child[i] <= 0) {
		i++;
		continue;
	    }

	    ch = get_bin(io, bin->child[i]);
	    if (start >= offset + ch->pos &&
		end   <= offset + ch->pos + ch->size-1) {
		cache_decr(io, bin);
		bin = ch;
		cache_incr(io, ch);
		offset += bin->pos;
		i = 0; /* restart loop */
	    } else {
		i++;
	    }
	}

	if (!extend) {
	    if (offset_r)
		*offset_r = offset;
	    return bin;
	}

	/*
	 * We now have the smallest bin available holding this range, but 
	 * as we delay creation of the sub-bins until required then we
	 * should perhaps create a child bin.
	 */
	if (bin->size > MIN_BIN_SIZE && (!bin->child[0] || !bin->child[1])) {
	    /* Construct left child if needed */
	    if (!bin->child[0]) {
		int pos = 0;
		int sz = bin->size/2;

		if (bin->child[1]) {
		    ch = get_bin(io, bin->child[1]);
		    sz = ch->pos;
		}

		if (start >= offset + pos &&
		    end   <= offset + pos + sz-1) {
		    /* It would fit - create it and continue recursing */

		    if (!(bin = cache_rw(io, bin)))
			return NULL;

		    bin->child[0] = bin_new(io, pos, sz, bin->rec, GT_Bin);
		    bin->flags |= BIN_BIN_UPDATED;
		    cache_decr(io, bin);

		    bin = get_bin(io, bin->child[0]);
		    offset += bin->pos;
		    continue;
		}
	    }

	    /* Construct right child if needed */
	    if (!bin->child[1]) {
		int pos = bin->size/2;
		int sz;

		if (bin->child[0]) {
		    ch = get_bin(io, bin->child[0]);
		    pos = ch->size;
		}
		sz = bin->size - pos;

		if (start >= offset + pos &&
		    end   <= offset + pos + sz-1) {
		    /* It would fit - create it and continue recursing */

		    if (!(bin = cache_rw(io, bin)))
			return NULL;

		    bin->child[1] = bin_new(io, pos, sz, bin->rec, GT_Bin);
		    bin->flags |= BIN_BIN_UPDATED;
		    cache_decr(io, bin);

		    bin = get_bin(io, bin->child[1]);
		    offset += bin->pos;
		    continue;
		}
	    }
	}
	
	/* At the smallest bin already */
	break;
    }

    if (offset_r)
	*offset_r = offset;

    last_io = io;
    last_c = (*c)->rec;
    last_bin = bin;
    last_start = offset;
    last_end = offset + bin->size-1;
    last_offset = offset;

    cache_decr(io, bin);
    return bin;
}

/*
 * Adds a range to the contig.
 *
 * Returns the bin we added the range to on success
 *         NULL on failure
 */
bin_index_t *bin_add_range(GapIO *io, contig_t **c, range_t *r,
			   range_t **r_out) {
    bin_index_t *bin;
    range_t *r2;
    int nr, offset;

    if (contig_get_start(c) > r->start)
	contig_set_start(io, c, r->start);

    if (contig_get_end(c) < r->end)
	contig_set_end(io, c, r->end);

    if (!(bin = bin_for_range(io, c, r->start, r->end, 1, &offset)))
	return NULL;

    if (!(bin = cache_rw(io, bin)))
	return NULL;

    /* Adjust start/end used in bin */
    if (bin->start_used != bin->end_used) {
	if (bin->start_used > r->start - offset)
	    bin->start_used = r->start - offset;
	if (bin->end_used < r->end - offset)
	    bin->end_used = r->end - offset;
    } else {
	bin->start_used = r->start - offset;
	bin->end_used = r->end - offset;
    }

    /* Update Range array */
    bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;
    if (!bin->rng)
	bin->rng = ArrayCreate(sizeof(range_t), 0);

    nr = next_range(bin->rng);
    r2 = arrp(range_t, bin->rng, nr);
    *r2 = *r; /* struct copy */
    r2->start -= offset;
    r2->end -= offset;

    if (r_out)
	*r_out = r2;

    return bin;
}

/*
 * Removes a record referenced from a bin.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_item(GapIO *io, contig_t **c, int rec) {
    bin_index_t *bin;
    int i, start, end;
    int cnum, pos;
    seq_t *s;

    s = (seq_t *)cache_search(io, GT_Seq, rec);
    if (-1 == sequence_get_position(io, rec, &cnum, &pos, NULL))
	return -1;

    start = pos;
    end   = pos + (s->len > 0 ? s->len : -s->len) - 1;

    if (!(bin = bin_for_range(io, c, start, end, 0, NULL)))
	return -1;

    if (!(bin = cache_rw(io, bin)))
	return -1;

    /* FIXME: we should check and update bin->start_used and bin->end_used */

    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);
	if (r->rec != rec)
	    continue;

	/*
	 * Found it, move it down or should we just label it as a hole?
	 * The latter is more efficient, but we need to implement the other
	 * end for this support.
	 */
	memmove(arrp(range_t, bin->rng, i), arrp(range_t, bin->rng, i+1), 
		(ArrayMax(bin->rng) - (i+1)) * sizeof(range_t));
	ArrayMax(bin->rng)--;
	break;
    }

    return 0;
}

/*
 * Finds the contig number and position of the start of a bin.
 * The right position is obviously this + bin->size (regardless of
 * whether it has been complemented).
 *
 * Returns 0 on success (plus *contig & *pos)
 *        -1 on failure
 */
int bin_get_position(GapIO *io, bin_index_t *bin, int *contig, int *pos) {
    int bnum;
    int offset = 0;

    /* FIXME: Complemented coordinates need more fixes here */
    if (bin->flags & BIN_COMPLEMENTED) {
	offset = bin->size-1 - offset;
    }

    offset += bin->pos;

    /* Find the position of this bin relative to the contig itself */
    while (bin->parent_type == GT_Bin) {
	bnum = bin->parent;
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
	offset += bin->pos;
    }
    
    assert(bin->parent_type == GT_Contig);
    *contig = bin->parent;

    *pos = offset;

    return 0;
}

/*
 * Some tracks are "official" cached objects and are deallocated as part of
 * the tg_cache system. Others are temporary on-the-fly structs generated
 * for the purpose of temporary display. These we need to free.
 *
 * This function knows which is which and does the appropriate thing.
 */
void track_free(track_t *t) {
    if (t->flag & TRACK_FLAG_FREEME) {
	if (t->data)
	    ArrayDestroy(t->data);
	free(t);
    }
}

/*
 * Creates a fake track struct, to be freed with track_free
 *
 * Returns track_t on success
 *         NULL on failure
 */
track_t *track_create_fake(int type, int size) {
    track_t *t = (track_t *)calloc(1, sizeof(*t));
    if (!t)
	return 0;
    t->type = type;
    t->nitems = size;
    t->item_size = sizeof(int);
    t->data = ArrayCreate(sizeof(int), size);
    t->flag |= TRACK_FLAG_FREEME;

    return t;
}

/*
 * Creates a track of a given type for this bin.
 * Note, this does not actually add it to the bin (but probably should
 * otherwise it's nothing more than the non-bin track_create).
 *
 * Returns track_t pointer on success
 *         NULL on failure
 */
track_t *bin_create_track(GapIO *io, bin_index_t *bin, int type) {
    int rec;
    track_t *t;

    if (-1 == (rec = io->iface->track.create(io->dbh, NULL)))
	return NULL;

    t = (track_t *)cache_search(io, GT_Track, rec);
    track_set_type(io, &t, type);

    return t;
}

/*
 * Adds a track to this bin.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_add_track(GapIO *io, bin_index_t **bin, track_t *track) {
    bin_index_t *n;
    int i;
    GBinTrack *bt;

    if (!(n = cache_rw(io, *bin)))
	return -1;
    *bin = n;

    /* Create new bin-track, or error if already found */
    if (!n->track) {
	n->track = ArrayCreate(sizeof(GBinTrack), 0);
	n->flags |= BIN_TRACK_UPDATED;
    }

    for (i = 0; i < ArrayMax(n->track); i++) {
	bt = arrp(GBinTrack, n->track, i);
	if (bt->type == track->type)
	    return -1;
    }

    /* Add the track pointer */
    bt = (GBinTrack *)ArrayRef(n->track, ArrayMax(n->track));
    bt->type = track->type;
    bt->flags = 0;
    bt->rec = track->rec;

    return 0;
}

/*
 * Finds the track of a given type for this bin.
 *
 * Returns track_t pointer on success (do not free)
 *         NULL on failure (eg bin too small)
 */
track_t *bin_get_track(GapIO *io, bin_index_t *bin, int type) {
    int i;

    /* If it exists and is up to date, return it */
    if (bin->track) {
	for (i = 0; i < ArrayMax(bin->track); i++) {
	    GBinTrack *bt = arrp(GBinTrack, bin->track, i);
	    if (bt->type == type) {
		return (track_t *)cache_search(io, GT_Track, bt->rec);
	    }
	}
    }

    return NULL;
}

/*
 * A bit like bin_get_track, but this is designed to auto-generate and
 * update the track as desired. The expectation is that this will always
 * succeed and anything else is a fatal error.
 */
track_t *bin_query_track(GapIO *io, bin_index_t *bin, int type) {
    track_t *bin_recalculate_track(GapIO *io, bin_index_t *bin, int type);
    int i;

    /* If it exists and is up to date, return it */
    if (bin->track) {
	for (i = 0; i < ArrayMax(bin->track); i++) {
	    GBinTrack *bt = arrp(GBinTrack, bin->track, i);
	    if (bt->type == type /*&& (bt->flags&1) != 0*/)
		return (track_t *)cache_search(io, GT_Track, bt->rec);
	}
    }

    /* Otherwise generate and maybe cache */
    return bin_recalculate_track(io, bin, type);
}

/*
 * ---------------------------------------------------------------------------
 * Track handling
 */

#define RD_ELEMENTS 1024
track_t *bin_recalculate_track(GapIO *io, bin_index_t *bin, int type) {
    int pos, cnum;
    track_t *track, *child;
    int nele;
    double bpv;
    contig_t *c;

    /*
     * So we have a bin of a given size in which we wish to have at least
     * RD_ELEMENTS of track samples, but it's out of date. We query the
     * contig for track data at a resolution of double what we need and
     * then average/downsample it to generate our new bin stats.
     */
    bpv = (double)bin->size / RD_ELEMENTS;
    if (bpv < 1) bpv = 1;
    nele = bin->size / bpv;
    if (nele & 1) nele++;
    bpv = (double)bin->size / nele;

    /*
     * Bottom layer, so no point querying a child - we just create it
     * ourselves now.
     */
    if (bpv <= 1) {
	track_t *fake;
	int rec, *depth;
	fake = track_create_fake(type, bin->size);
	fake->flag = TRACK_FLAG_VALID | TRACK_FLAG_FREEME;

	/* FIXME: type specific code here - or the size at least (int) */
	switch (type) {
	case TRACK_READ_DEPTH:
	    depth = track_read_depth_r1(io, bin);
	    break;
	default:
	    fprintf(stderr, "Unknown track type %d\n", type);
	    return NULL;
	}
	memcpy(ArrayBase(int, fake->data), depth, bin->size * sizeof(int));
	free(depth);

	rec = io->iface->track.create(io->dbh, fake);
	track = (track_t *)cache_search(io, GT_Track, rec);

	printf("Initialising track %d in bin %d at bpv of 1\n",
	       rec, bin->rec);

	bin_add_track(io, &bin, track);
	track_free(fake);

	return track;
    }

    /* Else use child bin tracks, in ever-decreasing circles */
    if (-1 == bin_get_position(io, bin, &cnum, &pos))
	return NULL;
    c = (contig_t *)cache_search(io, GT_Contig, cnum);
    child = contig_get_track(io, &c, pos, pos + bin->size-1, type, bpv);
    if (NULL == child)
	return NULL;

    track = bin_get_track(io, bin, type);
    if (!track) {
	track = bin_create_track(io, bin, type);
	bin_add_track(io, &bin, track);
    }
    cache_incr(io, track);

    /* Copy child 'fake track' to our real track */
    track_set_data(io, &track, ArrayCreate(sizeof(int), nele));
    track_set_nitems(io, &track, nele);
    track_set_item_size(io, &track, sizeof(int));
    memcpy(ArrayBase(int, track->data), ArrayBase(int, child->data),
	   nele * sizeof(int));

    track_free(child);

    track->flag |= TRACK_FLAG_VALID;
    cache_decr(io, track);

    return track;
}
