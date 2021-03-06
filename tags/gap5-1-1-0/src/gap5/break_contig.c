#include <stdio.h>
#include <string.h>

#include "break_contig.h"
#include "misc.h"

#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

static int break_contig_move_bin(GapIO *io, bin_index_t *bin,
				 contig_t *cfrom, int pfrom,
				 contig_t *cto,   int pto,
				 int child_no) {
    /* Add to */
    if (pto == cto->rec) {
	/* Parent is a contig */
	if (bin->rec != cto->bin)
	    printf("Destroy old bin for contig %d (bin %d). New root=%d\n",
		   cto->rec, cto->bin, bin->rec);
	cto->bin = bin->rec;
	cto->start = 1;
	cto->end = bin->size;

	bin->parent = cto->rec;
	bin->parent_type = GT_Contig;
	bin->flags |= BIN_BIN_UPDATED;

    } else {
	/* Parent is a bin */
	bin_index_t *pb;

	if (!(pb = get_bin(io, pto)))
	    return -1;
	if (!(pb = cache_rw(io, pb)))
	    return -1;

	pb->child[child_no] = bin->rec;
	pb->flags |= BIN_BIN_UPDATED;

	bin->parent = pto;
	bin->parent_type = GT_Bin;
	bin->flags |= BIN_BIN_UPDATED;
    }

    /* Remove from: NB it may not exist? */
    if (pfrom == cfrom->rec) {
	/* Parent is a contig */
	if (cfrom->bin != bin->rec) {
	    fprintf(stderr, "pfrom incorrect\n");
	    return -1;
	}

	cfrom->bin = 0;
    } else if (pfrom > 0) {
	/* Parent is a bin */
	bin_index_t *pb;

	if (!(pb = get_bin(io, pfrom)))
	    return -1;
	if (!(pb = cache_rw(io, pb)))
	    return -1;

	if (pb->child[0] != bin->rec && pb->child[1] != bin->rec) {
	    fprintf(stderr, "pfrom incorrect\n");
	    return -1;
	}

	if (!(pb = cache_rw(io, pb)))
	    return -1;
	
	if (pb->child[0] == bin->rec)
	    pb->child[0] = 0;
	else
	    pb->child[1] = 0;
	pb->flags |= BIN_BIN_UPDATED;
    }

    return 0;
}

/*
 * Given ranges contained within a bin this makes sure that all sequences
 * referred to in these ranges have their parent listed as the new bin.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int break_contig_reparent_seqs(GapIO *io, bin_index_t *bin) {
    int i, nr = bin->rng ? ArrayMax(bin->rng) : 0;

    for (i = 0; i < nr; i++) {
	range_t *r = arrp(range_t, bin->rng, i);
	seq_t *seq = (seq_t *)cache_search(io, GT_Seq, r->rec);
	if (seq->bin != bin->rec) {
	    seq = cache_rw(io, seq);
	    seq->bin = bin->rec;
	}
    }

    return 0;
}

/*
 * A recursive break contig function.
 * bin_num	The current bin being moved or split.
 * pos		The contig break point.
 * offset	The absolute positional offset of this bin in original contig
 * pleft	The parent bin/contig record num in the left new contig
 * pright	The parent bin/contig record num in the right new contig
 * child_no     0 or 1 - whether this bin is the left/right child of its parent
 */
static int break_contig_recurse(GapIO *io,
				contig_t *cl, contig_t *cr,
				int bin_num, int pos, int offset,
				int level, int pleft, int pright,
				int child_no, int complement,
				int *left_end, int *right_start) {
    int i, f_a, f_b, rbin;
    bin_index_t *bin = get_bin(io, bin_num), *bin_dup ;
    int bin_min, bin_max;

    cache_incr(io, bin);

    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }

    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    printf("%*sBreak offset %d pos %d => test bin %d: %d..%d\n",
	   level*4, "",
	   offset, pos, bin->rec,
	   NMIN(bin->start_used, bin->end_used),
	   NMAX(bin->start_used, bin->end_used));

    bin = cache_rw(io, bin);
    bin_invalidate_track(io, bin, TRACK_ALL);

    bin_min = bin->rng ? NMIN(bin->start_used, bin->end_used) : offset;
    bin_max = bin->rng ? NMAX(bin->start_used, bin->end_used) : offset;

    /*
     * Add to right parent if this bin is to the right of pos,
     * or if the used portion is to the right and we have no left child.
     *
     * FIXME: Not a valid assumption!
     * The used portion of a bin is not a placeholder for the used portion
     * of all the the children beneath it. Therefore if the used portion of
     * this bin is > pos (and we have no left child) it still doesn't mean
     * that the absolute positions of the used portion of the right child
     * won't be < pos.
     */
    if (offset >= pos /*|| (bin_min >= pos && !bin->child[0])*/) {
	printf("%*sADD_TO_RIGHT pl=%d pr=%d\n", level*4, "", pleft, pright);
	if (0 != break_contig_move_bin(io, bin,
				       cl, pleft, cr, pright, 
				       child_no))
	    return -1;

	if (*right_start > bin_min)
	    *right_start = bin_min;

	cache_decr(io, bin);
	return 0;
    }

    /*
     * Add to left parent if this bin is entirely to the left of pos,
     * or if the used portion is to the left and we have no right child.
     */
    if (offset + bin->size < pos /*|| (bin_max < pos && !bin->child[1])*/) {
	printf("%*sADD_TO_LEFT\n", level*4, "");

	//if (0 != break_contig_move_bin(io, bin, cr, pright, cl, pleft, child_no))
	//return -1;
	if (*left_end < bin_max)
	    *left_end = bin_max;

	cache_decr(io, bin);
	return 0;
    }

    /*
     * Nominally the bin overlaps both left and right and so needs duplicating.
     * There are cases though at the roots of our trees where duplicating is
     * unnecessary as it leads to empty bins at the root. In this case
     * we skip creating a duplicate for the right, or alternatively steal
     * the left root bin and use that instead.
     *
     * Similarly the range_t array will either be left where it is, moved to
     * the right contig, or split in half (creating a new one for the right).
     */
    if (pright != cr->rec ||
	(bin->rng && NMAX(bin->start_used, bin->end_used) >= pos)) {
	//printf("NMAX=%d >= %d\n", NMAX(bin->start_used, bin->end_used), pos);

	rbin = 0;

	/* Possibly steal left contig's bin */
	if (pleft == cl->rec && NMIN(bin->start_used, bin->end_used) >= pos) {
#if 0
	    /* Currently this doesn't always work */
	    if (bin->child[1]) {
		bin_index_t *ch = get_bin(io, bin->child[1]);
		if (NMIN(ch->pos, ch->pos + ch->size-1) >= pos) {
		    rbin = cl->bin;
		    cl->bin = bin->child[0];
		}
	    }
#else
	    pleft = bin->rec;
#endif
	} else {
	    pleft = bin->rec;
	}

	/* Create new bin, or use root of contig if it's unused so far */
	if (!rbin && pright == cr->rec) {
	    rbin = cr->bin;
	}

	/* Otherwise we genuingly need a duplicate */
	if (!rbin)
	    rbin = bin_new(io, 0, 0, 0, GT_Bin);

	/* Initialise with duplicate values from left bin */
	bin_dup = get_bin(io, rbin);
	bin_dup = cache_rw(io, bin_dup);
	bin_dup->size = bin->size;
	bin_dup->pos = bin->pos;
	bin_dup->parent = pright;
	bin_dup->parent_type = (pright == cr->rec ? GT_Contig : GT_Bin);
	bin_dup->flags = bin->flags | BIN_BIN_UPDATED;
	bin_dup->start_used = bin->start_used;
	bin_dup->end_used = bin->end_used;

	/*
	 * Shift bin by offset if it's the contig root.
	 * It'll be shifted back by the correct amount later.
	 */
	if (pright == cr->rec) {
	    printf("moving root bin via offset=%d comp=%d\n", offset, complement);
	    bin_dup->pos += offset;
	}

	printf("%*sCreated dup for right, rec %d\n", level*4,"", bin_dup->rec);
	break_contig_move_bin(io, bin_dup, cl, 0, cr, pright, child_no);
	pright = bin_dup->rec;
    } else {
	bin_dup = NULL;
	pleft = bin->rec;
    }

    if (!bin->rng) {
	/* Empty bin */
	printf("%*sEMPTY range\n", level*4, "");
	bin->start_used = bin->end_used = 0;
	bin->flags |= BIN_BIN_UPDATED;
	if (bin_dup) {
	    bin_dup->start_used = bin_dup->end_used = 0;
	    bin_dup->flags |= BIN_BIN_UPDATED;
	}
	    
    } else if (NMIN(bin->start_used, bin->end_used) >= pos) {
	/* Move range to right contig */
	printf("%*sDUP %d, MOVE Array to right\n", level*4, "", bin_dup->rec);

	bin_dup->rng = bin->rng;
	bin_dup->rng_rec = bin->rng_rec;
	if (bin_dup->rng_rec)
	    bin_dup->flags |= BIN_RANGE_UPDATED;

	if (bin->rec != bin_dup->rec) {
	    bin->rng = NULL;
	    bin->rng_rec = 0;
	    bin->flags |= BIN_BIN_UPDATED;
	}

	if (*right_start > NMIN(bin->start_used, bin->end_used))
	    *right_start = NMIN(bin->start_used, bin->end_used);

	bin->start_used = bin->end_used = 0;
	break_contig_reparent_seqs(io, bin_dup);

    } else if (NMAX(bin->start_used, bin->end_used) < pos) {
	/* Range array already in left contig, so do nothing */
	printf("%*sMOVE Array to left\n", level*4, "");

	if (*left_end < NMAX(bin->start_used, bin->end_used))
	    *left_end = NMAX(bin->start_used, bin->end_used);

	if (bin_dup)
	    bin_dup->start_used = bin_dup->end_used = 0;

    } else {
	/* Range array covers pos, so split in two */
	int i, j, n;
	int lmin = bin->size, lmax = 0, rmin = bin->size, rmax = 0;
	printf("%*sDUP %d, SPLIT array\n", level*4, "", bin_dup->rec);

	bin->flags |= BIN_RANGE_UPDATED;
	bin_dup->flags |= BIN_RANGE_UPDATED;

	bin_dup->rng = ArrayCreate(sizeof(range_t), 0);
	n = ArrayMax(bin->rng);
	for (i = j = 0; i < n; i++) {
	    range_t *r = arrp(range_t, bin->rng, i), *r2;
	    seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	    int cstart; /* clipped sequence positions */

	    if ((s->len < 0) ^ complement) {
		cstart = NMAX(r->start, r->end) - (s->right-1);
		/* cend   = NMAX(r->start, r->end) - (s->left-1); */
	    } else {
		cstart = NMIN(r->start, r->end) + s->left-1;
		/* cend   = NMIN(r->start, r->end) + s->right-1; */
	    }
	    
	    if (cstart >= pos)  {
		r2 = (range_t *)ArrayRef(bin_dup->rng, ArrayMax(bin_dup->rng));
		*r2 = *r;
		if (rmin > r->start) rmin = r->start;
		if (rmin > r->end)   rmin = r->end;
		if (rmax < r->start) rmax = r->start;
		if (rmax < r->end)   rmax = r->end;
	    } else {
		if (lmin > r->start) lmin = r->start;
		if (lmin > r->end)   lmin = r->end;
		if (lmax < r->start) lmax = r->start;
		if (lmax < r->end)   lmax = r->end;

		if (j != i) {
		    r2 = arrp(range_t, bin->rng, j);
		    *r2 = *r;
		}
		j++;
	    }
	}

	ArrayMax(bin->rng) = j;

	printf("%d seqs => left=%ld, right=%ld\n",
	       n, ArrayMax(bin->rng), bin_dup ? ArrayMax(bin_dup->rng) : 0);

	if (bin_dup)
	    break_contig_reparent_seqs(io, bin_dup);

	if (lmin < lmax) {
	    bin->start_used     = lmin;
	    bin->end_used       = lmax;
	} else {
	    /* No data left in bin */
	    bin->start_used = 0;
	    bin->end_used = 0;
	}

	printf("%*sLeft=>%d..%d right=>%d..%d\n", level*4, "",
	       lmin, lmax, rmin, rmax);

	if (bin_dup) {
	    if (rmin < rmax) {
		bin_dup->start_used = rmin;
		bin_dup->end_used   = rmax;
	    } else {
		/* No data moved in bin */
		bin_dup->start_used = 0;
		bin_dup->end_used   = 0;
	    }
	}

	if (lmin < lmax) {
	    if (*left_end < NMAX(bin->start_used, bin->end_used))
		*left_end = NMAX(bin->start_used, bin->end_used);
	}
	if (bin_dup && rmin < rmax) {
	    if (*right_start > NMIN(bin_dup->start_used, bin_dup->end_used))
		*right_start = NMIN(bin_dup->start_used, bin_dup->end_used);
	}
    }


    /* Recurse */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;

	ch = get_bin(io, bin->child[i]);
	if (0 != break_contig_recurse(io, cl, cr, bin->child[i], pos,
				      NMIN(ch->pos, ch->pos + ch->size-1),
				      level+1, pleft, pright,
				      i, complement,
				      left_end, right_start))
	    return -1;
    }

    cache_decr(io, bin);
    if (bin_dup)
	cache_decr(io, bin_dup);

    return 0;
}

/*
 * Looks for redundant bins at the root containing no data and just a single
 * child.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int remove_redundant_bins(GapIO *io, contig_t *c) {
    int bnum;

    if (!(c = cache_rw(io, c)))
	return -1;

    for (bnum = c->bin; bnum;) {
	bin_index_t *bin = get_bin(io, bnum);
	if (bin->rng || (bin->child[0] && bin->child[1]))
	    break;

	/* Empty */
	c->bin = bin->child[0] ? bin->child[0] : bin->child[1];
	printf("Remove bin %d\n", bin->rec);
	bnum = c->bin;
    }

    return 0;
}

/*
 * Breaks a contig in two such that snum is the right-most reading of
 * a new contig.
 */
int break_contig(GapIO *io, int crec, int cpos) {
    contig_t *cl = (contig_t *)cache_search(io, GT_Contig, crec);
    contig_t *cr;
    int cid;
    char cname[1024], *cname_end;
    int left_end, right_start;
    bin_index_t *bin;
    int do_comp = 0;

    contig_dump_ps(io, &cl, "/tmp/tree.ps");

    strncpy(cname, contig_get_name(&cl), 1000);
    cname_end = cname + strlen(cname);
    cid = 1;
    do {
	sprintf(cname_end, "#%d", cid++);
    } while (contig_index_query(io, cname) > 0);

    if (!(cr = contig_new(io, cname)))
	return -1;
    cr = cache_rw(io, cr);
    if (0 != contig_index_update(io, cname, strlen(cname), cr->rec))
	return -1;
    printf("Break in contig %d, pos %d\n", crec, cpos);

    printf("Existing left bin = %d, right bin = %d\n",
	   cl->bin, cr->bin);

    bin = get_bin(io, cl->bin);
    do_comp = bin->flags & BIN_COMPLEMENTED;
    cache_decr(io, bin);

    left_end = 0;
    right_start = INT_MAX;
    break_contig_recurse(io, cl, cr,
			 contig_get_bin(&cl), cpos, contig_offset(io, &cl),
			 0, cl->rec, cr->rec, 0, 0, &left_end, &right_start);

    printf("New left end = %d, right start = %d\n", left_end, right_start);

    /* Ensure start/end positions of contigs work out */
    bin = cache_rw(io, get_bin(io, cr->bin));
    //#define KEEP_POSITIONS 1
#ifndef KEEP_POSITIONS
    cr->start = 1;
    cr->end = cl->end - right_start + 1;
    bin->pos -= right_start-1;
#endif
    if ((do_comp && !(bin->flags & BIN_COMPLEMENTED)) ||
	(!do_comp && (bin->flags & BIN_COMPLEMENTED))) {
	bin->flags ^= BIN_COMPLEMENTED;
    }

#ifdef KEEP_POSITIONS
    cr->start = right_start;
    cr->end = cl->end;
#endif

    cl = cache_rw(io, cl);
    cl->end = left_end;

    remove_redundant_bins(io, cl);
    remove_redundant_bins(io, cr);

    cache_flush(io);

    printf("Final left bin = %d, right bin = %d\n",
	   cl->bin, cr->bin);

    contig_dump_ps(io, &cl, "/tmp/tree_l.ps");
    contig_dump_ps(io, &cr, "/tmp/tree_r.ps");

    return 0;
}
