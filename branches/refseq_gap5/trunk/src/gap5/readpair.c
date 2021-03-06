/*
 * File: readpair.c:
 * Version: 2.0 (1.0 == FORTRAN version)
 *
 * Author: Staden Package group
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: "Find read pair" code.
 *
 * Created: 17 March 1994
 * Updated:
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>

#include "xalloc.h"
#include "IO.h"
#include "gap4_compat.h"
#include "cs-object.h"
#include "contig_selector.h"
#include "newgap_cmds.h"
#include "gap_globals.h"
#include "misc.h"
#include "text_output.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "readpair.h"
#include "editor_view.h"
#include "tk-io-reg.h"
#include "io_lib/hash_table.h"

typedef struct {
    tg_rec template; /* maybe 0 if only a direct pair */
    // tg_rec library;  /* ? want this */
    tg_rec rec[2];
    int start[2], end[2];
    int contig[2];
    int mqual[2];
} read_pair_t;

/*
 * Match callback.
 * 'obj' is a match contained within the 'repeat' list.
 */
void *readpair_obj_func(int job, void *jdata, obj_read_pair *obj,
			mobj_template *template) {
    static char buf[80];
    char cmd[1024];
    f_int *handle;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(template->io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(template->io, cs_id);

    switch(job) {
    case OBJ_LIST_OPERATIONS:
	if (io_rdonly(template->io) && ((obj->c1 > 0 && obj->c2 < 0) ||
					(obj->c1 < 0 && obj->c2 > 0))) {
	    return "Information\0Hide\0IGNORE\0"
		"IGNORE\0SEPARATOR\0Remove\0";
	} else {
	    return "Information\0Hide\0Invoke join editor *\0SEPARATOR\0Remove\0";
	}
	break;

    case OBJ_INVOKE_OPERATION:
	switch(*((int *)jdata)) {
	case 0: /* Information */
	    vfuncgroup(1, "2D plot matches");
	case -1: /* Information from results manager */
	    start_message();

	    vmessage("Read pair:\n");
	    vmessage("    From contig %s(#%"PRIrec") at %d reading %s(#%"PRIrec")\n",
		     get_contig_name(template->io, ABS(obj->c1)),
		     io_clnbr(template->io, ABS(obj->c1)), obj->pos1,
		     get_read_name(template->io, obj->read1), obj->read1);
	    vmessage("    With contig %s(#%"PRIrec") at %d reading %s(#%"PRIrec")\n",
		     get_contig_name(template->io, ABS(obj->c2)),
		     io_clnbr(template->io, ABS(obj->c2)), obj->pos2,
		     get_read_name(template->io, obj->read2), obj->read2);
	    {
		seq_t *s;;

		s = cache_search(template->io, GT_Seq, obj->read1);
		vmessage("    Direction of first read is %swards\n",
			 (s->flags & SEQ_END_MASK) == SEQ_END_FWD
			 ? "for" : "back");

		s = cache_search(template->io, GT_Seq, obj->read2);
		vmessage("    Direction of second read is %swards\n",
			 (s->flags & SEQ_END_MASK) == SEQ_END_FWD
			 ? "for" : "back");
	    }
	    vmessage("    Length %d\n\n", obj->length);
	    end_message(cs->window);
	    break;

	case 1: /* Hide */
	    obj_hide(GetInterp(), cs->window, (obj_match *)obj,
		     (mobj_repeat *)template, csplot_hash);
	    break;

	case -2: /* default */
	case 2: /* Join editor */
	    {
	        tg_rec cnum[2], llino[2];
		int pos[2];
		int pt1, pt2;

		cnum[0] = ABS(obj->c1);
		cnum[1] = ABS(obj->c2);

		/* Complement a contig if needed */
		if ((obj->c1 > 0) != (obj->c2 > 0)) {
		    if (cnum[0] == cnum[1]) {
			verror(ERR_WARN, "join_editor",
			       "cannot display the same contig in two "
			       "different orientations");
			break;
		    }

		    if (-1 == complement_contig(template->io, cnum[0]))
			if (-1 == complement_contig(template->io, cnum[1]))
			    return NULL;
		}

		pos[0] = 1;
		pos[1] = 1;
		join_contig(template->io, cnum, llino, pos);

		break;
	    }

	case 3: /* Remove */
	    obj_remove(GetInterp(), cs->window, (obj_match *)obj,
		       (mobj_repeat *)template, csplot_hash);
	    break;

	}
	break;

    case OBJ_GET_BRIEF:
	sprintf(buf,
		"Read pair: %c#%d@%d (mq %d) with %c#%d@%d (mq %d), len %d",
		obj->c1 > 0 ? '+' : '-', obj->read1, obj->pos1, obj->mq1,
		obj->c2 > 0 ? '+' : '-', obj->read2, obj->pos2, obj->mq2,
		obj->length);
	return buf;
    }

    return NULL;
}

void readpair_callback(GapIO *io, tg_rec contig, void *fdata, reg_data *jdata) {
    mobj_template *r = (mobj_template *)fdata;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id);

    switch (jdata->job) {

    case REG_QUERY_NAME:

	sprintf(jdata->name.line, "Find read pairs");
	break;


    case REG_JOIN_TO:

	csmatch_join_to(io, contig, &jdata->join, (mobj_repeat *)r,
			csplot_hash, cs->window);
	break;


    case REG_COMPLEMENT:

	csmatch_complement(io, contig, r, csplot_hash, cs->window);
	break;


    case REG_GET_OPS:

	if (r->all_hidden)
	    jdata->get_ops.ops = "PLACEHOLDER\0PLACEHOLDER\0Information\0"
		"PLACEHOLDER\0Hide all\0Reveal all\0SEPARATOR\0Remove\0";
	else
	    jdata->get_ops.ops = "Use for 'Next'\0Reset 'Next'\0Information\0"
		"Configure\0Hide all\0Reveal all\0SEPARATOR\0Remove\0";
	break;


    case REG_INVOKE_OP:

	switch (jdata->invoke_op.op) {
	case 0: /* Next */
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(r), NULL);
	    break;
	case 1: /* Reset Next */
	    csmatch_reset_next((mobj_repeat *)r);
	    break;
	case 2: /* Information */
	    csmatch_info((mobj_repeat *)r, "Read pair");
	    break;
	case 3: /* Configure */
	    csmatch_configure(io, cs->window, (mobj_repeat *)r);
	    break;
	case 4: /* Hide all */
	    csmatch_hide(GetInterp(), cs->window, (mobj_repeat *)r,
			 csplot_hash);
	    break;
	case 5: /* Reveal all */
	    csmatch_reveal(GetInterp(), cs->window, (mobj_repeat *)r,
			   csplot_hash);
	    break;
	case 6: /* Remove */
	    csmatch_remove(io, cs->window,
			   (mobj_repeat *)r,
			   csplot_hash);
	    break;
	}
	break;


    case REG_PARAMS:

	jdata->params.string = r->params;
	break;


    case REG_NUMBER_CHANGE:

	csmatch_renumber(io, contig, jdata->number.number,
			 (mobj_repeat *)r, csplot_hash, cs->window);
	break;


    case REG_ORDER:

	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;

    case REG_QUIT:

	csmatch_remove(io, cs->window,
		       (mobj_repeat *)r,
		       csplot_hash);
	break;


    case REG_DELETE:

	csmatch_contig_delete(io, (mobj_repeat *)r, contig,
			      cs->window, csplot_hash);
	break;

    case REG_LENGTH:
	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;
    }
}

/*
 * Plot the templates which span more than one contig on contig selector plot.
 * The direction of the line is governed by the "direction" of the template
 * A line going from top right to bottom left (/) indicates the second contig
 * needs completing
 * A line going from top left to bottom right (\) indicates the second contig
 * does not need completing
 */
int PlotTempMatches(GapIO *io, read_pair_t *rp) {
    int i, j, k, id;
    mobj_template *template;
    obj_read_pair *matches;
    int count;
    int n_matches = 0;
    int max_matches = 64; /* dynamically grows */
    item_t *item, *it;
    char *val;

    if (!rp)
	return 0;

    if (NULL == (template = (mobj_template *)xmalloc(sizeof(mobj_template))))
	return -1;
    if (NULL == (matches = (obj_read_pair *)xmalloc(max_matches *
						    sizeof(obj_read_pair))))
	return -1;

    /* Create cs-object plot array */
    while (rp->rec[0]) {
	matches[n_matches].func =
	    (void *(*)(int, void *, struct obj_match_t *,
		       struct mobj_repeat_t *))readpair_obj_func;
	matches[n_matches].data = template;
	/* Fix dir on c1/c2 */
	matches[n_matches].c1 = rp->contig[0];
	matches[n_matches].c2 = rp->contig[1];
	matches[n_matches].pos1 = rp->start[0];
	matches[n_matches].pos2 = rp->start[1];
	matches[n_matches].length = (ABS(rp->end[0] - rp->start[0]) + 
				     ABS(rp->end[1] - rp->start[1])) / 2;
	matches[n_matches].read1 = rp->rec[0];
	matches[n_matches].read2 = rp->rec[1];
	matches[n_matches].flags = 0;
	matches[n_matches].mq1 = rp->mqual[0];
	matches[n_matches].mq2 = rp->mqual[1];
	if (++n_matches >= max_matches) {
	    obj_read_pair *tmp;
	    
	    max_matches *= 2;
	    tmp = (obj_read_pair *)
		xrealloc(matches,
			 max_matches * sizeof(obj_read_pair));
	    if (NULL == tmp) {
		xfree(template);
		xfree(matches);
		return -1;
	    } else {
		matches = tmp;
	    }
	}

	rp++;
    }


    /*
     * Free memory and return if we've got no matches.
     */
    if (0 == n_matches) {
	xfree(template);
	xfree(matches);

	return 0;
    }

    template->num_match = n_matches;
    template->match = (obj_match *)matches;
    template->io = io;
    strcpy(template->tagname, CPtr2Tcl(template));

    val = get_default_string(GetInterp(), gap5_defs, "READPAIR.COLOUR");
    strcpy(template->colour, val);

    template->linewidth = get_default_int(GetInterp(), gap5_defs,
					  "READPAIR.LINEWIDTH");

    template->params = (char *)xmalloc(10);
    if (template->params)
	sprintf(template->params, "none");
    template->all_hidden = 0;
    template->current = -1;
    template->reg_func = readpair_callback;
    template->match_type = REG_TYPE_READPAIR;

    PlotRepeats(io, template);
    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(template), NULL);

    /*
     * Register the repeat search with each of the contigs used.
     * Currently we assume that this is all.
     */
    id = register_id();
    contig_register(io, 0, readpair_callback, (void *)template, id,
		    REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
		    REG_NUMBER_CHANGE | REG_ORDER, REG_TYPE_READPAIR);
    update_results(io);

    return 0;
}

/*
 * Identifies spanning read pairs within end_size of the end of each contig.
 * Mode is 0 for all vs all.
 *         1 for ends vs all.
 *	   2 for ends vs ends.
 */
read_pair_t *spanning_pairs(GapIO *io, int num_contigs,
			    contig_list_t *contig_array,
			    enum readpair_mode mode,
			    int end_size, int min_mq) {
    int i;
    HashTable *h;
    HashIter *iter;
    HashItem *hi;
    read_pair_t *pairs = NULL;
    int npairs = 0, nalloc = 0;
    int no_large_contigs = 1;
    pool_alloc_t *rp_pool;

    h = HashTableCreate(1024, HASH_FUNC_HSIEH |
			      HASH_DYNAMIC_SIZE |
			      HASH_POOL_ITEMS);
    rp_pool = pool_create(sizeof(rangec_t));

    if (!h || !rp_pool)
	return NULL;

    for (i = 0; i < num_contigs; i++) {
	contig_iterator *ci;
	rangec_t *r;
	contig_t *c;
	int cstart, cend, crec = contig_array[i].contig;
	int large_contig;
	
	c = cache_search(io, GT_Contig, crec);
	cstart = c->start;
	cend   = c->end;

	/*
	 * For end vs end, can we optimise this into one linear scan through
	 * the contig instead of having an iterator at each end?
	 * Small contigs allow this (we just discard the data in the middle).
	 *
	 * For end vs all, we have to trade off scanning through data that
	 * isn't at the end vs finding the "all" part of the match?
	 * sequence_get_position() is very slow so we prefer to scan through
	 * bins far more here than brute force identification of the
	 * other end.
	 *
	 * For all vs all, we have to do all the work anyway so just
	 * face up to it.
	 */
	switch (mode) {
	default: /* To keep gcc happy */
	case all_all:
	    large_contig = 0;
	    break;

	case end_end:
	    large_contig = (cend-cstart)/4  >= end_size ? 1 : 0;
	    break;

	case end_all:
	    large_contig = (cend-cstart)/100 >= end_size ? 1 : 0;
	    break;
	}
	if (large_contig)
	    no_large_contigs = 0;

	/*
	 * Iterate through seqs in contig.
	 *
	 * Pair up reads, removing internal pairs as we go.
	 * We stop accumulating into our pair hash after end_size bases,
	 * but keep going to remove pairs until 2*end_size.
	 */
	ci = contig_iter_new(io, crec, 1, CITER_FIRST,
			     CITER_CSTART, CITER_CEND);

	while (r = contig_iter_next(io, ci)) {
	    if (!r->pair_rec)
		continue; /* unpaired */

	    if (large_contig && r->start - cstart > 2*end_size)
		break;

	    if ((hi = HashTableSearch(h, (char *)&r->pair_rec,
				      sizeof(r->pair_rec)))) {
		rangec_t *r2 = hi->data.p;
		if (r2->orig_rec == crec) {
		    /* internal read pair */
		    pool_free(rp_pool, r2);
		    HashTableDel(h, hi, 0);
		    continue;
		}
	    }

	    if (mode == all_all || mode == end_all ||
		r->start - cstart < end_size || cend - r->start < end_size) {
		HashData hd;
		rangec_t *r2 = pool_alloc(rp_pool);
		*r2 = *r;
		r2->orig_rec = crec; /* convenient place to store contig */
		if (r->start - cstart < end_size ||
		    cend - r->start < end_size)
		    r2->pair_start = 1; /* near end */
		else
		    r2->pair_start = 0;
		hd.p = r2;

		HashTableAdd(h, (char *)&r->rec, sizeof(r->rec), hd, 0);
	    }
	}

	contig_iter_del(ci);

	/* Other end, if not yet done (for large contigs) */
	if (large_contig) {
	    ci = contig_iter_new(io, crec, 1, CITER_FIRST,
				 MAX(cstart, cend - 2*end_size), CITER_CEND);
	    
	    while (r = contig_iter_next(io, ci)) {
		if (!r->pair_rec)
		    continue; /* unpaired */

		if ((hi = HashTableSearch(h, (char *)&r->pair_rec,
					  sizeof(r->pair_rec)))){
		    rangec_t *r2 = hi->data.p;
		    if (r2->orig_rec == crec) {
			/* internal read pair */
			pool_free(rp_pool, r2);
			HashTableDel(h, hi, 0);
		    }
		}

		if (mode == end_all || cend - r->start < end_size) {
		    HashData hd;
		    rangec_t *r2 = pool_alloc(rp_pool);
		    *r2 = *r;
		    r2->orig_rec = crec;
		    if (r->start - cstart < end_size ||
			cend - r->start < end_size)
			r2->pair_start = 1; /* near end */
		    else
			r2->pair_start = 0;
		    hd.p = r2;

		    HashTableAdd(h, (char *)&r->rec, sizeof(r->rec), hd, 0);
		}
	    }

	    contig_iter_del(ci);
	}
    }


    /*
     * Our hash table ideally now contains pairs with rec1 -> rec2 and 
     * rec2 -> rec1. This will be the case for mode == end_end or all_all.
     * Any single rec1 -> rec2 with no rec2 present in the hash will be due
     * to a read-pair outside of the region (end_end) or outside of the
     * set of contigs we chose to scan (all_all).
     *
     * For mode end_all we may have had some large contigs where we chose
     * not to scan through all reads, only effectively adding rec1 for the
     * ends of those contigs. In that case it is expected rec2 may be absent
     * and instead we need to run sequence_get_position instead. This is
     * slower, which is why we only resort to this route on extra-large
     * contigs.
     */
    iter = HashTableIterCreate();
    while ((hi = HashTableIterNext(h, iter))) {
	tg_rec rec1, rec2;
	rangec_t *r1 = (rangec_t *)hi->data.p, *r2, r2_tmp;

	rec1 = r1->rec;
	rec2 = r1->pair_rec;

	if (rec1 == 0)
	    continue; // already added other end */

	if (!(hi = HashTableSearch(h, (char *)&rec2, sizeof(rec1)))) {
	    tg_rec contig;
	    int start, end, orient;

	    if (mode != end_all || no_large_contigs)
		continue; // unpaired

	    /* Maybe rec2 was in a large contig that we didn't scan through */
	    sequence_get_position(io, rec2, &contig, &start, &end, &orient);
	    r2 = &r2_tmp;
	    r2->rec = rec2;
	    r2->orig_rec = contig;
	    r2->start = start;
	    r2->end = end;
	} else {
	    r2 = (rangec_t *)hi->data.p;
	}


	if (r1->orig_rec == r2->orig_rec)
	    /* Same contig */
	    continue;

	if (mode == end_all && !(r1->pair_start || r2->pair_start))
	    /* Spanning pair, but not with a match near an end */
	    continue;

	/* Mapped to a repeat? */
	if (r1->mqual < min_mq || r2->mqual < min_mq)
	    continue;

	/* Accepted - add it to pairs[] array */
	if (npairs >= nalloc) {
	    nalloc = nalloc ? nalloc*2 : 1024;
	    pairs = realloc(pairs, nalloc * sizeof(*pairs));
	    if (!pairs)
		return NULL;
	}

	pairs[npairs].rec[0] = rec1;
	pairs[npairs].rec[1] = rec2;
	pairs[npairs].contig[0] = r1->orig_rec;
	pairs[npairs].contig[1] = r2->orig_rec;
	pairs[npairs].start[0] = r1->start;
	pairs[npairs].start[1] = r2->start;
	pairs[npairs].end[0] = r1->end;
	pairs[npairs].end[1] = r2->end;
	pairs[npairs].mqual[0] = r1->mqual;
	pairs[npairs].mqual[1] = r2->mqual;

	r2->rec = 0;

	npairs++;
    }

    /* Indicate end of list */
    if (npairs >= nalloc) {
	nalloc = nalloc ? nalloc*2 : 1024;
	pairs = realloc(pairs, nalloc * sizeof(*pairs));
	if (!pairs)
	    return NULL;
    }
    pairs[npairs].rec[0] = 0;
	
    HashTableIterDestroy(iter);
    HashTableDestroy(h, 0);
    pool_destroy(rp_pool);

    vmessage("Found %d read-pairs\n", npairs);

    return pairs;
}


/*
 * External callable functions - the C interface
 * ---------------------------------------------
 */
int find_read_pairs(GapIO *io, int num_contigs, contig_list_t *contig_array,
		    enum readpair_mode mode, int end_size, int min_mq) {

    read_pair_t *tarr;

    if (NULL == (tarr = spanning_pairs(io, num_contigs, contig_array,
				       mode, end_size, min_mq)))
	return -1;

    /* Find only those templates spanning multiple contigs and plot them. */
    PlotTempMatches(io, tarr);

    /* Tidy up */
    free(tarr);

    return 0;
}
