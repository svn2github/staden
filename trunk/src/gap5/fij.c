#include <tk.h>
#include <assert.h>

#include "io_utils.h"
#include "gap_globals.h"
#include "fij.h"
#include "contig_selector.h"
#include "newgap_cmds.h" /* GetInterp() */
#include "misc.h"
#include "text_output.h"
#include "consen.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"
#include "gap4_compat.h"
#include "editor_view.h"
#include "tk-io-reg.h"
#include "readpair.h"
#include "editor_join.h"
#include "io_lib/hash_table.h"
#include "dstring.h"

static int counter;
static int counter_max;
static mobj_fij *global_match;

static int auto_join(GapIO *io, mobj_fij *r);

void *fij_obj_func(int job, void *jdata, obj_fij *obj,
		      mobj_fij *fij) {
    static char buf[160];
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(fij->io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(fij->io, cs_id);

    switch(job) {
    case OBJ_LIST_OPERATIONS:
	if (io_rdonly(fij->io) && ((obj->c1 > 0 && obj->c2 < 0) ||
				   (obj->c1 < 0 && obj->c2 > 0))) {
	    return "Information\0Hide\0IGNORE\0"
		"IGNORE\0SEPARATOR\0Remove\0";
	} else {
	    return "Information\0Hide\0IGNORE\0Invoke join editor *\0"
		"Invoke contig editors\0SEPARATOR\0Remove\0";
	}

    case OBJ_INVOKE_OPERATION:
	switch(*((int *)jdata)) {
	case 0: /* Information */
	    vfuncgroup(1, "2D plot matches");
	case -1: /* Information from results manager */

	    start_message();
	    vmessage("FIJ match\n");
	    vmessage("    From contig %s(=%"PRIrec") at %d\n",
		     get_contig_name(fij->io, ABS(obj->c1)),
		     ABS(obj->c1), obj->pos1);
	    vmessage("    With contig %s(=%"PRIrec") at %d\n",
		     get_contig_name(fij->io, ABS(obj->c2)),
		     ABS(obj->c2), obj->pos2);
	    vmessage("    Length %d, score %d, mismatch %2.2f%%\n\n",
		     obj->length, obj->score, ((float)obj->percent)/10000);
	    end_message(cs->window);
	    break;

	case 1: /* Hide */
	    obj_hide(GetInterp(), cs->window, (obj_match *)obj,
		     (mobj_repeat *)fij, csplot_hash);
	    break;

	case 2: /* Make join */ {
	    printf("Make join between %"PRIrec" and %"PRIrec"\n",
		   obj->c1, obj->c2);
	    break;
	}

	case -2: /* default */
        case 3: /* Invoke join editor */ {
	    tg_rec cnum[2], llino[2];
	    int pos[2];
	    int shortest;

	    obj->flags |= OBJ_FLAG_VISITED;
	    fij->current = obj - (obj_fij *)fij->match;
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(fij), NULL);

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
		if (io_rdonly(fij->io)) {
		    bell();
		    break;
		}
		shortest = (io_clength(fij->io, cnum[0])
			    < io_clength(fij->io, cnum[1])) ? 0 : 1;
		if (-1 == complement_contig(fij->io, cnum[shortest]))
		    if (-1 == complement_contig(fij->io, cnum[1 - shortest]))
			    return NULL;
	    }

	    /*
	     * NB: obj->pos1 now may not be the same value as when this
	     * function was entered, due to the complementing!
	     */
	    pos[0] = obj->pos1; 
	    pos[1] = obj->pos2;
	    
	    llino[0] = 0;
	    llino[1] = 0;

	    join_contig(fij->io, cnum, llino, pos);
	    break;
	}

	case 4: /* Invoke contig editors */ {
	    tg_rec cnum;
	    int pos, reveal;
	    
	    cnum  = ABS(obj->c1);
	    pos   = obj->pos1;
	    reveal= (obj->pos1 <= 0 ||
		     obj->pos2 <= 0 ||
		     obj->pos1 >= io_clength(fij->io, ABS(obj->c1)) ||
		     obj->pos2 >= io_clength(fij->io, ABS(obj->c2))) ? 1 : 0;

	    edit_contig(fij->io, cnum, 0, pos);

	    cnum  = ABS(obj->c2);
	    pos   = obj->pos2;
	    edit_contig(fij->io, cnum, 0, pos);
	    break;
	}

	case 5: /* Remove */
	    obj_remove(GetInterp(), cs->window, (obj_match *)obj,
		       (mobj_repeat *)fij, csplot_hash);
	    break;

        }
        break;

    case OBJ_GET_BRIEF:
	sprintf(buf,
		"FIJ: %c=%"PRIrec"@%d with %c=%"PRIrec"@%d, "
		"len %d, score %d, mis %2.2f%%",
		obj->c1 > 0 ? '+' : '-', ABS(obj->c1), obj->pos1,
		obj->c2 > 0 ? '+' : '-', ABS(obj->c2), obj->pos2,
		obj->length, obj->score, ((float)obj->percent)/10000);
	return buf;
    }

    return NULL;
}

static int sort_func_score(const void *p1, const void *p2) {
    obj_fij *m1 = (obj_fij *)p1, *m2 = (obj_fij *)p2;
    return m2->score - m1->score;
}

static int sort_func_c1(const void *p1, const void *p2) {
    obj_fij *m1 = (obj_fij *)p1, *m2 = (obj_fij *)p2;
    return ABS(m2->c1) - ABS(m1->c1);
}

static int sort_func_c2(const void *p1, const void *p2) {
    obj_fij *m1 = (obj_fij *)p1, *m2 = (obj_fij *)p2;
    return ABS(m2->c2) - ABS(m1->c2);
}

/*
 * Match callback.
 * 'obj' is a match contained within the 'fij' list.
 */
void fij_callback(GapIO *io, tg_rec contig, void *fdata, reg_data *jdata) {
    mobj_fij *r = (mobj_fij *)fdata;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id);

    switch (jdata->job) {

    case REG_QUERY_NAME:

	sprintf(jdata->name.line, "Find Internal Joins");
	break;


    case REG_JOIN_TO:

	csmatch_join_to(io, contig, &jdata->join, (mobj_repeat *)r,
			csplot_hash, cs ? cs->window : NULL);
	break;


    case REG_COMPLEMENT:

	csmatch_complement(io, contig, (mobj_repeat *)r,
			   csplot_hash, cs ? cs->window : NULL);
	break;


    case REG_GET_OPS:

	if (r->all_hidden)
	    jdata->get_ops.ops = "PLACEHOLDER\0PLACEHOLDER\0Information\0"
		"PLACEHOLDER\0Hide all\0Reveal all\0Sort matches\0"
		    "Save matches\0SEPARATOR\0Remove\0";
	else
	    jdata->get_ops.ops = "Use for 'Next'\0Reset 'Next'\0Information\0"
		"Configure\0Hide all\0Reveal all\0Sort matches\0"
		    "Save matches\0SEPARATOR\0Remove\0";
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
	    csmatch_info((mobj_repeat *)r, "Find Internal Joins");
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
	case 6: /* Sort */
	    qsort(r->match, r->num_match, sizeof(obj_fij), sort_func_score);
	    csmatch_reset_hash(csplot_hash, (mobj_repeat *)r);
	    r->current = -1;
	    break;
	case 7: { /* Save */
	    char *fn;
	    if (Tcl_VarEval(GetInterp(), "tk_getSaveFile ", "-parent ",
			    cs->window, NULL) != TCL_OK)
		break;
	    fn = Tcl_GetStringResult(GetInterp());
	    if (fn && *fn)
		csmatch_save((mobj_generic *) r, fn);
	    break;
	}
	case 8: /* Remove */
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
			 (mobj_repeat *)r, csplot_hash,
			 cs ? cs->window : NULL);
	break;

    case REG_ORDER:
#ifdef DEBUG
	printf("find internal joins REG_ORDER\n");
#endif
	csmatch_replot(io, (mobj_repeat *)r, csplot_hash,
		       cs ? cs->window : NULL);
	break;

    case REG_QUIT:

	csmatch_remove(io, cs ? cs->window : NULL,
		       (mobj_repeat *)r,
		       csplot_hash);
	break;

    case REG_DELETE:

	csmatch_contig_delete(io, (mobj_repeat *)r, contig,
			      cs ? cs->window : NULL, csplot_hash);
	break;

    case REG_LENGTH:
	csmatch_replot(io, (mobj_repeat *)r, csplot_hash,
		       cs ? cs->window : NULL);
	break;

    case REG_GENERIC:
	switch (jdata->generic.task) {
	    int ret;

	case TASK_CS_PLOT:
	    PlotRepeats(io, (mobj_repeat *)r);
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(r), NULL);
	    break;

	case TASK_CS_SAVE:
	    ret = csmatch_save((mobj_generic *) r, (char *)jdata->generic.data);
	    vTcl_SetResult(GetInterp(), "%d", ret);
	    break;

	case TASK_CS_LOAD:
	    break;

	case TASK_AUTO_JOIN:
	    auto_join(io, r);
	    break;
	}
    }
}

static int cl_compare(const void *va, const void *vb) {
    contig_list_t *a = (contig_list_t *) va;
    contig_list_t *b = (contig_list_t *) vb;
    if (a->contig < b->contig) return -1;
    if (a->contig > b->contig) return  1;
    return 0;
}

static contig_list_t * fij_get_contig_list(GapIO *io,
					   int num1, contig_list_t *contigs1,
					   int num2, contig_list_t *contigs2,
					   int *num_contigs_out,
					   int *list1end_out,
					   int *list2start_out) {
    contig_list_t *combined;
    int i1, i2, num_shared, o1, o2, os;
    tg_rec last_c = 0;

    qsort(contigs1, num1, sizeof(contig_list_t), cl_compare);
    qsort(contigs2, num2, sizeof(contig_list_t), cl_compare);

    /* Ensure contigs1 and contigs2 contain no duplicates */
    for (i1 = o1 = 0; i1 < num1; i1++) {
	if (contigs1[i1].contig == last_c) continue;
	if (o1 != i1) contigs1[o1] = contigs1[i1];
	last_c = contigs1[o1++].contig;
    }
    num1 = o1;

    for (i2 = o2 = 0; i2 < num2; i2++) {
	if (contigs2[i2].contig == last_c) continue;
	if (o2 != i2) contigs2[o2] = contigs2[i2];
	last_c = contigs2[o2++].contig;
    }
    num2 = o2;

    /* Count number of contigs in both contigs1 and contigs2 */
    for (i1 = i2 = num_shared = 0; i1 < num1 && i2 < num2; ) {
	if (contigs1[i1].contig == contigs2[i2].contig) {
	    num_shared++;
	    i1++;
	    i2++;
	} else if (contigs1[i1].contig < contigs2[i2].contig) {
	    i1++;
	} else {
	    i2++;
	}
    }

    /* Sort contigs into combined array, in order:
     * contigs only in contigs1
     * contigs in both contigs1 and contigs2
     * contigs only in contigs2
     */
    combined = xmalloc(sizeof(contig_list_t) * (num1 + num2 - num_shared));
    if (NULL == combined) return NULL;

    o1 = 0;
    os = num1 - num_shared;
    o2 = num1;
    for (i1 = i2 = 0; i1 < num1 && i2 < num2; ) {
	if (contigs1[i1].contig == contigs2[i2].contig) {
	    combined[os++] = contigs1[i1++];
	    i2++;
	} else if (contigs1[i1].contig < contigs2[i2].contig) {
	    combined[o1++] = contigs1[i1++];
	} else {
	    combined[o2++] = contigs2[i2++];
	}
    }
    while (i1 < num1) { combined[o1++] = contigs1[i1++]; }
    while (i2 < num2) { combined[o2++] = contigs2[i2++]; }
    
    *num_contigs_out = num1 + num2 - num_shared;
    *list1end_out = num1;
    *list2start_out = num1 - num_shared;

    return combined;
}

static Contig_parms * contig_list2contig_parms(GapIO *io, int num_contigs,
					       contig_list_t *list) {

    /* Make a contig list */
    Contig_parms *contig_list;
    contig_list = get_contig_list(io, num_contigs, list);
    
    return contig_list;
}

static HashTable * fij_prefilter_repeats(fij_arg *fij_args,
					 HashTable *lib_hash, int *num_contigs,
					 contig_list_t *list, int *list1end,
					 int *list2start) {
    read_pair_t *pairs;
    read_pair_t *pair;
    HashTable *links = NULL;
    HashTable *index = NULL;
    uint8_t *used = NULL;
    int l1e = *list1end, l2s = *list2start;
    int i, j, count = 0;

    pairs = spanning_pairs(fij_args->io, *num_contigs, list,
			   fij_args->rp_mode, fij_args->rp_end_size,
			   fij_args->rp_min_mq, fij_args->rp_min_freq,
			   lib_hash);
    if (NULL == pairs) return NULL;
    links = HashTableCreate(1024, HASH_DYNAMIC_SIZE | HASH_OWN_KEYS);
    if (NULL == links) goto fail;
    index = HashTableCreate(*num_contigs, HASH_NONVOLATILE_KEYS);
    if (NULL == index) goto fail;
    used = calloc((*num_contigs >> 3) + 1, sizeof(uint8_t));
    if (NULL == used) goto fail;

    for (i = 0; i < *num_contigs; i++) {
	HashData hd;
	hd.i = i;
	if (!HashTableAdd(index, (char *)&list[i].contig,
			  sizeof(list[i].contig), hd, NULL)) goto fail;
    }

    for (pair = pairs; pair->rec[0]; pair++) {
	contig_pair cp;
	HashData hd;
	HashItem *hi;
	int idx1, idx2;
	tg_rec c1 = ABS(pair->contig[0]), c2 = ABS(pair->contig[1]);
	int new = 0;

	if (c1 == c2) continue;
	hi = HashTableSearch(index, (char *) &c1, sizeof(c1));
	assert(NULL != hi);
	idx1 = hi->data.i;
	hi = HashTableSearch(index, (char *) &c2, sizeof(c2));
	assert(NULL != hi);
	idx2 = hi->data.i;
	if (!((idx1 < *list1end && idx2 >= *list2start)
	      || (idx2 < *list1end && idx1 >= *list2start))) continue;

	cp.c1 = MIN(c1, c2);
	cp.c2 = MAX(c1, c2);
	hd.i = 0;
	if (!HashTableAdd(links, (char *)&cp, sizeof(cp), hd, &new)) goto fail;
	used[idx1 >> 3] |= 1 << (idx1 & 7);
	used[idx2 >> 3] |= 1 << (idx2 & 7);
	if (new) count++;
    }

    vmessage("%d contig pairs are linked by read pairs\n", count);

    for (i = j = 0; i < *num_contigs; i++) {
	if (0 == (used[i >> 3] & (1 << (i & 7)))) {
	    if (i < *list1end)   --l1e;
	    if (i < *list2start) --l2s;
	    continue;
	}
	if (i != j) list[j] = list[i];
	j++;
    }
    vmessage("After first stage read-pair filter:\n"
	     "  List 1 includes %d contigs (was %d)\n"
	     "  List 2 includes %d contigs (was %d)\n",
	     l1e, *list1end, j - l2s, *num_contigs - *list2start);
    *num_contigs = j;
    *list1end = l1e;
    *list2start = l2s;
    HashTableDestroy(index, 0);
    free(used);
    free(pairs);
    return links;

 fail:
    if (NULL != pairs) free(pairs);
    if (NULL != links) HashTableDestroy(links, 0);
    if (NULL != index) HashTableDestroy(index, 0);
    if (NULL != used)  free(used);
    return NULL;
}


/*
 * Removes FIJ matches based on the number of matching contig ends.
 * If end of contig A joins to start of contig B, but start of contig B joins
 * to both contig A and contig C, then remove the links.
 */
static void fij_postfilter_matches(mobj_fij *FIJMatch) {
    int i, j;
    obj_fij *m;
    GapIO *io = FIJMatch->io;

    vmessage("\nFiltering joins based on duplicate matches per contig end.\n");

    /* Sort by 1st contig */
    qsort(FIJMatch->match, FIJMatch->num_match, sizeof(obj_fij), sort_func_c1);

    /* Remove left or right ends going to the same spot */
    m = FIJMatch->match;
    for (i = j = 0; i < FIJMatch->num_match; ) {
	int k = i, start = 0, end = 0;
	tg_rec c1 = ABS(m[i].c1);
	int cstart1, cend1;
	tg_rec c2 = ABS(m[i].c2);
	int cstart2, cend2;

	consensus_valid_range(io, c1, &cstart1, &cend1);
	consensus_valid_range(io, c2, &cstart2, &cend2);

	/* Count */
	while (k < FIJMatch->num_match && ABS(m[i].c1) == ABS(m[k].c1)) {
	    if (m[k].pos1 <= cstart1)
		start++;
	    if (m[k].end1 >= cend1)
		end++;
	    k++;
	}

	//printf("\n1> ctg=%"PRIrec" start=%d end=%d\n",
	//       ABS(m[i].c1), start, end);

	/* Filter */
	k = i;
	while (k < FIJMatch->num_match && ABS(m[i].c1) == ABS(m[k].c1)) {
	    //printf("%+6"PRIrec" %d (%6d..%6d) %6d /"
	    //       "%+6"PRIrec" %d (%6d..%6d) %6d\t",
	    //       m[k].c1, cstart1, m[k].pos1, m[k].end1, cend1,
	    //       m[k].c2, cstart2, m[k].pos2, m[k].end2, cend2);
	    //if (m[k].pos1 <= cstart1)
	    //    putchar('S');
	    //if (m[k].end1 >= cend1)
	    //    putchar('E');
	    if ((m[k].pos1 <= cstart1 && start > 1) ||
		(m[k].end1 >= cend1   && end   > 1)) {
		vmessage("  Rejecting #%+"PRIrec" / #%+"PRIrec"\n",
			 m[k].c1, m[k].c2);
		//putchar('*');
	    } else {
		m[j++] = m[k];
	    }
	    //putchar('\n');
	    k++;
	}

	i = k;
    }

    //fprintf(stderr, "Filtered from %d to %d matches\n", FIJMatch->num_match, j);

    FIJMatch->num_match = j;


    /* Sort by 2nd contig */
    qsort(FIJMatch->match, FIJMatch->num_match, sizeof(obj_fij), sort_func_c2);

    /* Remove left or right ends going to the same spot */
    m = FIJMatch->match;
    for (i = j = 0; i < FIJMatch->num_match; ) {
	int k = i, start = 0, end = 0;
	tg_rec c1 = ABS(m[i].c1);
	int cstart1, cend1;
	tg_rec c2 = ABS(m[i].c2);
	int cstart2, cend2;

	consensus_valid_range(io, c1, &cstart1, &cend1);
	consensus_valid_range(io, c2, &cstart2, &cend2);

	/* Count */
	while (k < FIJMatch->num_match && ABS(m[i].c2) == ABS(m[k].c2)) {
	    if (m[k].pos2 <= cstart2)
		start++;
	    if (m[k].end2 >= cend2)
		end++;
	    k++;
	}

	//printf("\n2> ctg=%"PRIrec" start=%d end=%d\n",
	//       ABS(m[i].c1), start, end);

	/* Filter */
	k = i;
	while (k < FIJMatch->num_match && ABS(m[i].c2) == ABS(m[k].c2)) {
	    //printf("%+6"PRIrec" %d (%6d..%6d) %6d /"
	    //	   "%+6"PRIrec" %d (%6d..%6d) %6d\t",
	    //	   m[k].c1, cstart1, m[k].pos1, m[k].end1, cend1,
	    //	   m[k].c2, cstart2, m[k].pos2, m[k].end2, cend2);
	    //if (m[k].pos2 <= cstart2)
	    //	putchar('S');
	    //if (m[k].end2 >= cend2)
	    //	putchar('E');
	    if ((m[k].pos2 <= cstart2 && start > 1) ||
		(m[k].end2 >= cend2   && end   > 1)) {
		//putchar('*');
		vmessage("  Rejecting #%"PRIrec" / #%"PRIrec"\n",
			 m[k].c1, m[k].c2);
	    } else {
		m[j++] = m[k];
	    }
	    //putchar('\n');
	    k++;
	}

	i = k;
    }

    //fprintf(stderr, "Filtered from %d to %d matches\n", FIJMatch->num_match, j);
    vmessage("Filtered out %d of %d matches, leaving %d\n",
	     FIJMatch->num_match - j, FIJMatch->num_match, j);

    FIJMatch->num_match = j;
}


/*
 * Constructs a parameter string.
 */
static char *fij_params(fij_arg *a) {
    char *s;
    dstring_t *ds = NULL;

    if (!(ds = dstring_create("Join searching parameters:\n")))
	goto err;
    if (-1 == dstring_append(ds, "  Join types allowed:             "))
			      
	goto err;
    if (a->ends && -1 == dstring_append(ds, " Ends"))
	goto err;
    if (a->containments && -1 == dstring_append(ds, " Containments"))
	goto err;
    dstring_append(ds, "\n");

    if (a->min_match == 0) {
	if (-1 == dstring_appendf(ds,
				  "  Alignment mode:                  "
				  "Sensitive\n"))
	    goto err;
    } else {
	if (-1 == dstring_appendf(ds,
				  "  Alignment mode:                  %s\n",
				  a->fast_mode ? "Fastest" : "Quick"))
	    goto err;

	if (-1 == dstring_appendf(ds,
				  "  Word length:                     %d\n",
				  a->word_len))
	    goto err;
	if (-1 == dstring_appendf(ds,
				  "  Minimum match:                   %d\n",
				  a->min_match))
	    goto err;
    }

    if (-1 == dstring_append(ds, "\nJoin filtering parameters:\n"))
	goto err;
    if (-1 == dstring_appendf(ds,
			      "  Minimum overlap:                 %d\n",
			      a->min_overlap))
	goto err;
    if (-1 == dstring_appendf(ds,
			      "  Maximum overlap:                 %d\n",
			      a->max_overlap))
	goto err;
    if (-1 == dstring_appendf(ds,
			      "  Maximum percentage mismatch:     %.1f%%\n",
			      a->max_mis))
	goto err;
    if (a->rp_mode == -1) {
	if (-1 == dstring_appendf(ds, "  Read-pair mode:                  "
				  "off\n"))
	    goto err;
    } else {
	if (-1 == dstring_appendf(ds,
				  "  Read-pair mode:                  %s\n",
	    a->rp_mode == end_end
	    ? "end_end"
	    : (a->rp_mode == end_all
	       ? "end_all"
	       : "all_all")))
	    goto err;
	if (-1 == dstring_appendf(ds,
				  "  Read-pair min. count:            %d\n",
				  a->rp_min_freq))
	    goto err;
	if (-1 == dstring_appendf(ds,
				  "  Read-pair min. percentage:       %d%%\n",
				  a->rp_min_perc))
	    goto err;
	if (-1 == dstring_appendf(ds,
				  "  Read-pair min. mapping score:    %d\n",
				  a->rp_min_mq))
	    goto err;
    }
    if (a->filter_words) {
	if (-1 == dstring_appendf(ds,
				  "  Filter repeat words >x times:    %.1f\n",
				  a->filter_words))
	    goto err;
    } else {
	if (-1 == dstring_appendf(ds,
				  "  Filter repeat words >x times:    off\n"))
	    goto err;
    }
    if (-1 == dstring_appendf(ds, "  Minimum depth:                   %d\n",
			      a->min_depth))
	goto err;
    if (-1 == dstring_appendf(ds, "  Maximum depth:                   %d\n",
			      a->max_depth ? a->max_depth : INT_MAX))
	goto err;
    if (-1 == dstring_appendf(ds, "  Unique contig joins only:        %s\n",
			      a->unique_ends ? "yes" : "no"))
	goto err;

    s = ds->str; ds->str = NULL; // HACK
    dstring_destroy(ds);

    return s;
    
 err:
    if (ds)
	dstring_destroy(ds);
    return NULL;
}


/*
 * Main find internal joins algorithm entry.
 *
 * All contigs in contig_array1 are compared against all contigs in
 * contig_array2.  Pairs of contigs that occur in both arrays will only
 * be compared once, and a contig will never be compared to itself.
 * 
 */
int
fij(fij_arg *fij_args,
    int num_contigs1,
    contig_list_t *contig_array1,
    int num_contigs2,
    contig_list_t *contig_array2)
{
    GapIO *io = fij_args->io;
    char *consensus = NULL;
    mobj_fij *FIJMatch = NULL;
    int id = 0;
    char *val;
    Contig_parms *contig_list = NULL;
    contig_list_t *combined = NULL;
    HashTable *links = NULL;
    HashTable *lib_hash = NULL;
    int number_of_contigs, list1end, list2start;
    int task_mask;
    int consensus_length, max_read_length;
    static char buf[80];
    Hidden_params p;
    int retval = -1;

    /* FIXME: get these from elsewhere */

    int gap_open, gap_extend;

    gap_open = gopenval; /* note we were using 8 and 6, not 8 and 7 */
    gap_extend = gextendval;

    p.min = p.max = p.verbose = 0;
    p.do_it = fij_args->use_hidden;
    p.use_conf = fij_args->use_conf;
    p.test_mode = 0;
    p.start = 0;
    p.lwin1 = 0;
    p.lcnt1 = 0;
    p.rwin1 = fij_args->win_size;
    p.rcnt1 = fij_args->dash;
    p.qual_val = fij_args->min_conf;
    p.window_len = fij_args->win_size;
    p.gap_open = gopenval;
    p.gap_extend = gextendval;
    p.band = 30; /*FIXME: hardwired 30 bases band for aligning hidden data */

    max_read_length = find_max_gel_len(io, 0, 0);

    if ((FIJMatch = (mobj_fij *)xmalloc(sizeof(mobj_fij)))==NULL){
	return -1;
    }

    counter_max = 14;
    if (NULL == (FIJMatch->match = (obj_fij *)xmalloc(counter_max *
						      sizeof(obj_fij)))) {
	goto out;
    }

    combined = fij_get_contig_list(io, num_contigs1, contig_array1,
				   num_contigs2, contig_array2,
				   &number_of_contigs,
				   &list1end, &list2start);
    if (NULL == combined) goto out;

    if (fij_args->rp_mode >= 0) {
	if (fij_args->rp_library) {
	    lib_hash = create_lib_hash(fij_args->rp_library,
				       fij_args->rp_nlibrary);
	    if (NULL == lib_hash) goto out;
	}
	links = fij_prefilter_repeats(fij_args, lib_hash,
				      &number_of_contigs, combined,
				      &list1end, &list2start);
	if (NULL == links) goto out;
    }

    contig_list = contig_list2contig_parms(io, number_of_contigs, combined);
    if (NULL == contig_list) {
	retval = -5;
	goto out;
    }

    FIJMatch->params = fij_params(fij_args);
    vmessage("\n%s\n", FIJMatch->params);

    global_match = FIJMatch;
    counter = 0;

    task_mask = ADDTITLE | NORMALCONSENSUS;
    if ( fij_args->mask == 2 ) task_mask |= MARKING;
    if ( fij_args->mask == 3 ) task_mask |= MASKING;
    if ( p.do_it ) task_mask |= ADDHIDDENDATA;

    consensus_length = 0;
    if ( make_consensus ( task_mask, io, &consensus, NULL,
			  contig_list, number_of_contigs,
			  &consensus_length,
			  max_read_length,
			  p,
			  consensus_cutoff ) ) {
	goto out;
    }

    if (do_it_fij(fij_args, consensus, consensus_length, gap_open, gap_extend,
		  COMPARE_ALL, lib_hash, links, contig_list, list1end,
		  contig_list + list2start, number_of_contigs - list2start,
		  list1end - list2start)) goto out;

    if (0 == counter) {
	vmessage("No joins found \n");
	retval = 0;
	goto out;
    }

    sprintf(buf, " Number of potential joins found   %d", counter);
    vmessage("%s\n", buf);

    FIJMatch->num_match = counter;
    FIJMatch->io = io;
    strcpy(FIJMatch->tagname, CPtr2Tcl(FIJMatch));

    val = get_default_string(GetInterp(), gap5_defs, "FIJ.COLOUR");
    strcpy(FIJMatch->colour, val);

    FIJMatch->linewidth = get_default_int(GetInterp(), gap5_defs,
					  "FIJ.LINEWIDTH");
    FIJMatch->all_hidden = 0;
    FIJMatch->current = -1;
    FIJMatch->reg_func = fij_callback;
    FIJMatch->match_type = REG_TYPE_FIJ;
    FIJMatch->max_mismatch = fij_args->max_mis;
    FIJMatch->min_length = fij_args->min_match;

    /* Filter matches */
    if (fij_args->unique_ends)
	fij_postfilter_matches(FIJMatch);

    /* Sort matches */
    qsort(FIJMatch->match, FIJMatch->num_match, sizeof(obj_fij),
	  sort_func_score);

    /*
     * Register find internal joins with each of the contigs used.
     * Currently we assume that this is all.
     */
    if (counter) {
	id = register_id();
	contig_register(io, 0, fij_callback, (void *)FIJMatch, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_NUMBER_CHANGE | REG_ORDER | REG_GENERIC,
			REG_TYPE_FIJ);
	update_results(io);
    }

    retval = id;
    FIJMatch = NULL; /* So we don't free it below */
 out:
    if (NULL != FIJMatch) {
	if (NULL != FIJMatch->match) xfree(FIJMatch->match);
	xfree(FIJMatch);
    }
    if (NULL != lib_hash) HashTableDestroy(lib_hash, 0);
    if (NULL != links) HashTableDestroy(links, 0);
    if (NULL != combined) xfree(combined);
    if (NULL != contig_list) xfree(contig_list);
    if (NULL != consensus) xfree(consensus);

    return retval;
} /* end fij */

/* store hits of find internal joins */
void
buffij(tg_rec c1, int pos1, int end1,
       tg_rec c2, int pos2, int end2,
       int len, int score, double percent)
{
    global_match->match[counter].func = fij_obj_func;
    global_match->match[counter].data = global_match;

    global_match->match[counter].c1 = c1;
    global_match->match[counter].pos1 = pos1;
    global_match->match[counter].end1 = end1;
    global_match->match[counter].c2 = c2;
    global_match->match[counter].pos2 = pos2;
    global_match->match[counter].end2 = end2;
    global_match->match[counter].length = len;
    global_match->match[counter].score = score;
    global_match->match[counter].percent = 10000 * percent;
    global_match->match[counter].flags = 0;

    if (++counter >= counter_max) {
	counter_max *= 2;
	global_match->match = (obj_fij *)xrealloc(global_match->match,
						    counter_max *
						    sizeof(obj_fij));
    }
}


/*
 * Automatically attempt to join the FIJ results 'result_id'.
 * This is based around the edJoinAlign() function in editor_join.c.
 *
 * Returns -1 on failure
 *          0 on success
 */
static int auto_join(GapIO *io, mobj_fij *r) {
    int i, j, nr, nr_orig;

    if (!r)
	return -1;

    vfuncheader("Auto-Join");

    nr_orig = r->num_match;
    for (i = j = 0; i < r->num_match; j++, i++) {
	obj_fij *m = &r->match[i];
	alignment_t *a;
	int l1, r1, l2, r2; /* contig used extents */
	int left1,right1,oleft1,oright1;
	int left2,right2,oleft2,oright2;
	int offset;
	int overlapLength;
	int len1,len2;
	int extra;

	vmessage("Processing join %d of %d: =%+"PRIrec" vs =%+"PRIrec"\n",
		 j+1, nr_orig, m->c1, m->c2);

	if (m->c1 < 0)
	    if (-1 == complement_contig(io, -m->c1))
		return -1;

	if (m->c2 < 0)
	    if (-1 == complement_contig(io, -m->c2))
		return -1;

	if (m->c1 < 0 || m->c2 < 0) {
	    /* Should never happen */
	    vmessage("Contig still not complemented - aborting\n");
	    continue;
	}

	/* Pick consensus regions to align over - adds a bit more */
	/* See edJoinAlign in editor_join.c */
	consensus_valid_range(io, ABS(m->c1), &l1, &r1);
	consensus_valid_range(io, ABS(m->c2), &l2, &r2);

	offset = (m->pos1 - l1) - (m->pos2 - l2);

	/* Possibly recompute overlap length - start/ends may have changed */
	if (m->flags & OBJ_FLAG_JOINED) {
	    int start, len;

	    start = MIN(m->pos1 - l1, m->pos2 - l2);
	    m->pos1 -= start;
	    m->pos2 -= start;

	    len = MIN(r1 - m->pos1, r2 - m->pos2);
	    oleft1 = left1 = m->pos1; oright1 = right1 = left1 + len;
	    oleft2 = left2 = m->pos2; oright2 = right2 = left2 + len;
	} else {
	    oleft1 = left1 = m->pos1; oright1 = right1 = m->end1;
	    oleft2 = left2 = m->pos2; oright2 = right2 = m->end2;
	}

	overlapLength = MAX(right1 - left1+1, right2 - left2+1);
	if (overlapLength <= 0) {
	    vmessage("Skipping as contigs seem to no longer overlap.\n");
	    continue;
	}

	/* Double check to avoid any bugs above. Must touch edges */
	if (left1 > l1 && left2 > l2) {
	    vmessage("BUG: Found local alignment - forcing to be global.\n");
	    int dist = MIN(left1-l1, left2-l2);
	    left1 -= dist;
	    left2 -= dist;
	}
	if (right1 < r1 && right2 < r2) {
	    vmessage("BUG: Found local alignment - forcing to be global.\n");
	    int dist = MIN(r1 - right1, r2 - right2);
	    right1 += dist;
	    right2 += dist;
	}

	/* Add on extra data either end to allow for padding */
	extra = set_band_blocks(overlapLength, overlapLength)/2;

	if ((left1  -= extra) < l1) left1  = l1;
	if ((left2  -= extra) < l2) left2  = l2;
	if ((right1 += extra) > r1) right1 = r1;
	if ((right2 += extra) > r2) right2 = r2;

	len1 = right1 - left1+1;
	len2 = right2 - left2+1;

	if (len1 <= 0 || len2 <= 0) {
	    vmessage("Skipping due to negative contig overlaps.\n");
	    continue;
	}

	/* Do the actual alignment */
	a = align_contigs(io, ABS(m->c1), left1, len1,
			  io, ABS(m->c2), left2, len2,
			  0, 0);

	if (!a) {
	    verror(ERR_WARN, "auto_join",
		   "Failed to align =%"PRIrec" with =%"PRIrec"\n",
		   ABS(m->c1), ABS(m->c2));
	    continue;
	}

	/* Check the alignment is still acceptable */
	if (a->match_len < r->min_length) {
	    vmessage("Skipping join as match too short.\n");
	    alignment_free(a);
	    continue;
	}

	if (100 - 100.0*a->match_count / a->match_len > r->max_mismatch) {
	    vmessage("Skipping join as percentage mismatch too high.\n");
	    alignment_free(a);
	    continue;
	}

	/* Apply the edits */
	align_apply_edits(io, ABS(m->c1), io, ABS(m->c2), a);

	/* And finally make the join */
	offset = left2 + a->shift - left1;
	vmessage("Join at offset %d\n", offset);

	alignment_free(a);

	nr = r->num_match;
	if (offset > 0) {
	    contig_t *c;
	    tg_rec crec = ABS(m->c2);

	    if (-1 == join_contigs(io, ABS(m->c2), ABS(m->c1), offset)) {
		vmessage("Join_contigs() failed\n");
		continue;
	    }

	    c = cache_search(io, GT_Contig, crec);
	    printf("Add to =%"PRIrec" at %d..%d\n", c->rec, oleft2, oright2);
	    anno_ele_add(io, GT_Contig, c->rec, 0, str2type("JOIN"), 
			 r->params, oleft2, oright2, ANNO_DIR_NUL);
	} else {
	    contig_t *c;
	    tg_rec crec = ABS(m->c1);

	    if (-1 == join_contigs(io, ABS(m->c1), ABS(m->c2), -offset)) {
		vmessage("Join_contigs() failed\n");
		continue;
	    }

	    c = cache_search(io, GT_Contig, crec);
	    printf("Add to =%"PRIrec" at %d..%d\n", c->rec, oleft1, oright1);
	    anno_ele_add(io, GT_Contig, c->rec, 0, str2type("JOIN"), 
			 r->params, oleft1, oright1, ANNO_DIR_NUL);
	}

	if (nr == 1)
	    break; // r will have been deallocated.

	cache_flush(io);

	i--;
    }

    cache_flush(io);

    vmessage("\nMade %d of %d joins found\n", nr_orig-i, nr_orig);

    return 0;
}

/*
 * Loads a file (already opened in fp) of join results and returns the
 * registered ID.
 *
 * Returns -1 on failure.
 */
int csmatch_load_fij(GapIO *io, FILE *fp) {
    tg_rec c1, c2;
    int pos1, pos2, end1, end2, length, score, n;
    float percent;
    int asize = 0;
    obj_fij *r;
    char *val;
    int id;

    mobj_fij *m = calloc(1, sizeof(*m));
    if (!m)
	return -1;

    strcpy(m->tagname, CPtr2Tcl(m));
    m->num_match = 0;
    m->match = NULL;
    m->io = io;
    m->all_hidden = 0;
    m->current = -1;
    val = get_default_string(GetInterp(), gap5_defs, "FIJ.COLOUR");
    strcpy(m->colour, val);
    m->linewidth = get_default_int(GetInterp(), gap5_defs, "FIJ.LINEWIDTH");
    m->match_type = REG_TYPE_FIJ;
    m->reg_func = fij_callback;

    while (9 == (n = fscanf(fp, "%"PRIrec" %d %d %"PRIrec" %d %d %d %d %f\n",
			    &c1, &pos1, &end1, &c2, &pos2, &end2,
			    &length, &score, &percent))) {
	contig_t *c;

	if (m->num_match >= asize) {
	    asize = asize ? asize*2 : 16;
	    m->match = realloc(m->match,
			       asize * sizeof(*m->match));
	    if (!m->match)
		return -1;
	}

	if (!cache_exists(io, GT_Contig, ABS(c1)) ||
	    !(c = cache_search(io, GT_Contig, ABS(c1)))) {
	    verror(ERR_WARN, "csmatch_load_fij",
		   "Contig =%"PRIrec" does not exist", ABS(c1));
	    continue;
	}

	if (pos1 < c->start)
	    pos1 = c->start;
	if (end1 > c->end)
	    end1 = c->end;

	if (!cache_exists(io, GT_Contig, ABS(c2)) ||
	    !(c = cache_search(io, GT_Contig, ABS(c2)))) {
	    verror(ERR_WARN, "csmatch_load_fij",
		   "Contig =%"PRIrec" does not exist", ABS(c2));
	    continue;
	}

	if (pos2 < c->start)
	    pos2 = c->start;
	if (end2 > c->end)
	    end2 = c->end;

	r = &m->match[m->num_match++];
	r->func = fij_obj_func;
	r->data = m;

	r->c1 = c1;
	r->c2 = c2;
	r->pos1 = pos1;
	r->pos2 = pos2;
	r->end1 = end1;
	r->end2 = end2;
	r->score = score;
	r->percent = 10000*percent;
	r->flags = 0; // fixme
    }
    if (n != EOF)
	verror(ERR_WARN, "csmatch_load_fij", "File malformatted or truncated");

    if (m->num_match) {
	id = register_id();
	contig_register(io, 0, fij_callback, (void *)m, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_NUMBER_CHANGE | REG_ORDER | REG_GENERIC,
			REG_TYPE_FIJ);
	update_results(io);

	return id;
    }

 err:
    if (m) {
	if (m->match)
	    free(m->match);

	free(m);
    }

    return -1;
}
