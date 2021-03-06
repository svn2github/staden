#include <math.h>
#include <ctype.h>

#include "io_utils.h"
#include "IO.h"
#include "tagdb.h"
#include "tagUtils.h"
#include "gap_globals.h"
#include "cs-object.h"
#include "contig_selector.h"
#include "contigEditor.h"
#include "misc.h"
#include "dna_utils.h"
#include "tcl_utils.h"
#include "text_output.h"
#include "complement.h"
#include "tclXkeylist.h"
#include "template_display.h"
#include "search_utils.h"
#include "sequence_formats.h"

#define TAG 0
#define SEQUENCE 1

/*
 * functions to act upon matches generated using TAGs
 */
void *find_oligo_obj_func1(int job,
			  void *jdata,
			  obj_match *obj,
			  mobj_find_oligo *find_oligo)
{
    static char buf[80];
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(find_oligo->io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(find_oligo->io, cs_id, 0);

    switch(job) {
    case OBJ_LIST_OPERATIONS:
	return "Information\0Hide\0Invoke join editor *\0"
	    "Invoke contig editors\0SEPARATOR\0Remove\0";

    case OBJ_INVOKE_OPERATION:
	switch(*((int *)jdata)) {
	case 0: /* Information */
	    vfuncgroup(1, "2D plot matches");
	case -1: /* Information from results manager */

	    start_message();
	    vmessage("Sequence search:\n");
	    vmessage("    From contig %s(#%d) at %d\n",
		     get_contig_name(find_oligo->io, ABS(obj->c1)),
		     io_clnbr(find_oligo->io, ABS(obj->c1)), obj->pos1);
	    vmessage("    With contig %s(#%d) at %d\n",
		     get_contig_name(find_oligo->io, ABS(obj->c2)),
		     io_clnbr(find_oligo->io, ABS(obj->c2)), obj->pos2);
	    vmessage("    Length %d, match %2.2f%%\n\n",
		     obj->length,
		     (float)obj->score / obj->length * 100.0 );
	    end_message(cs->window);
	    break;

	case 1: /* Hide */
	    obj_hide(GetInterp(), cs->window, obj,
		     (mobj_find_oligo *)find_oligo, csplot_hash);
	    break;

	case -2: /* default */
        case 2: /* Invoke join editor */ {
	    int cnum[2], llino[2], pos[2];

	    obj->flags |= OBJ_FLAG_VISITED;
	    find_oligo->current = obj - find_oligo->match;
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(find_oligo), NULL);

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

		if (io_rdonly(find_oligo->io)) {
		    bell();
		    break;
		}

		if (io_clength(find_oligo->io, ABS(obj->c1)) <
		    io_clength(find_oligo->io, ABS(obj->c2))) {
		    if (-1 == complement_contig(find_oligo->io, ABS(obj->c1)))
			if (-1 == complement_contig(find_oligo->io,
						    ABS(obj->c2)))
			    return NULL;
		} else {
		    if (-1 == complement_contig(find_oligo->io, ABS(obj->c2)))
			if (-1 == complement_contig(find_oligo->io,
						    ABS(obj->c1)))
			    return NULL;
		}
	    }

	    /*
	     * NB: obj->pos1 now may not be the same value as when this
	     * function was entered, due to the complementing!
	     */
	    pos[0] = obj->pos1;
	    pos[1] = obj->pos2;

	    llino[0] = io_clnbr(find_oligo->io, cnum[0]);
	    llino[1] = io_clnbr(find_oligo->io, cnum[1]);

	    join_contig(GetInterp(), find_oligo->io, cnum, llino, pos,
			consensus_cutoff, quality_cutoff);

	    break;
	}

	case 3: /* Invoke contig editors */ {
	    int cnum, llino, pos;

	    cnum  = ABS(obj->c1);
	    llino = io_clnbr(find_oligo->io, cnum);
	    pos   = obj->pos1;

	    edit_contig(GetInterp(), find_oligo->io, cnum, llino, pos,
			consensus_cutoff, quality_cutoff, 0, NULL);

	    cnum  = ABS(obj->c2);
	    llino = io_clnbr(find_oligo->io, cnum);
	    pos   = obj->pos2;

	    edit_contig(GetInterp(), find_oligo->io, cnum, llino, pos,
			consensus_cutoff, quality_cutoff, 0, NULL);
	    break;
	}

	case 4: /* Remove */
	    obj_remove(GetInterp(), cs->window, obj,
		       (mobj_find_oligo *)find_oligo, csplot_hash);
	    break;

        }
        break;

    case OBJ_GET_BRIEF:
	sprintf(buf,
		"Oligo: %c#%d@%d with %c#%d@%d, len %d, match %2.2f%%",
		obj->c1 > 0 ? '+' : '-',
		io_clnbr(find_oligo->io, ABS(obj->c1)), obj->pos1,
		obj->c2 > 0 ? '+' : '-',
		io_clnbr(find_oligo->io, ABS(obj->c2)), obj->pos2,
		obj->length, (float)obj->score / obj->length * 100.0);
	return buf;
    }

    return NULL;
}

/*
 * functions to act upon matches generated using SEQUENCE
 */
void *find_oligo_obj_func2(int job,
			  void *jdata,
			  obj_match *obj,
			  mobj_find_oligo *find_oligo)
{
    static char buf[80];
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(find_oligo->io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(find_oligo->io, cs_id, 0);

    switch(job) {
    case OBJ_LIST_OPERATIONS:
	return "Information\0Hide\0Invoke contig editor *\0"
	    "SEPARATOR\0Remove\0";

    case OBJ_INVOKE_OPERATION:
	switch(*((int *)jdata)) {
	case 0: /* Information */
	    vfuncgroup(1, "2D plot matches");
	case -1: /* Information from results manager */

	    start_message();
	    vmessage("Sequence search\n");
	    vmessage("    Contig %s(#%d) at %d\n",
		     get_contig_name(find_oligo->io, ABS(obj->c1)),
		     io_clnbr(find_oligo->io, ABS(obj->c1)), obj->pos1);
	    vmessage("    Length %d, match %2.2f%%\n\n",
		   obj->length, (float)obj->score / obj->length * 100.0 );
	    end_message(cs->window);
	    break;

	case 1: /* Hide */
	    obj_hide(GetInterp(), cs->window, obj,
		     (mobj_find_oligo *)find_oligo, csplot_hash);
	    break;

	case -2: /* default */
	case 2: /* Invoke contig editor */ {
	    int cnum, llino, pos, id;

	    obj->flags |= OBJ_FLAG_VISITED;
	    find_oligo->current = obj - find_oligo->match;

	    cnum  = ABS(obj->c1);
	    llino = 0;
	    pos   = obj->pos1;

	    if ((id = editor_available(cnum, 1)) == -1) {
		edit_contig(GetInterp(), find_oligo->io, cnum, llino, pos,
			    consensus_cutoff, quality_cutoff, 0, NULL);
	    }
	    if ((id = editor_available(cnum, 1)) != -1) {
		move_editor(id, llino, pos);
		editor_select_region(id, llino, pos, obj->length);
	    }
	    break;
	}

	case 3: /* Remove */
	    obj_remove(GetInterp(), cs->window, obj,
		       (mobj_find_oligo *)find_oligo, csplot_hash);
	    break;

        }
        break;

    case OBJ_GET_BRIEF:
	sprintf(buf,
		"Oligo: %c#%d@%d with %c#%d@%d, len %d, match %2.2f%%",
		obj->c1 > 0 ? '+' : '-',
		io_clnbr(find_oligo->io, ABS(obj->c1)), obj->pos1,
		obj->c2 > 0 ? '+' : '-',
		io_clnbr(find_oligo->io, ABS(obj->c2)), obj->pos2,
		obj->length, (float)obj->score / obj->length * 100.0);
	return buf;
    }

    return NULL;
}

static int sort_func(const void *p1, const void *p2) {
    obj_match *m1 = (obj_match *)p1, *m2 = (obj_match *)p2;
    return m2->score - m1->score;
}

/*
 * Match callback.
 * 'obj' is a match contained within the 'find_oligo' list.
 */
void find_oligo_callback(GapIO *io, int contig, void *fdata, reg_data *jdata) {
    mobj_find_oligo *r = (mobj_find_oligo *)fdata;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id, 0);

    switch (jdata->job) {

    case REG_QUERY_NAME:

	sprintf(jdata->name.line, "Find oligo");
	break;

    case REG_JOIN_TO:

	csmatch_join_to(io, contig, &jdata->join, (mobj_find_oligo *)r,
			csplot_hash, cs->window);
	break;

    case REG_COMPLEMENT:

	csmatch_complement(io, contig, r, csplot_hash, cs->window);
	break;


    case REG_GET_OPS:

	if (r->all_hidden)
	    jdata->get_ops.ops = "PLACEHOLDER\0PLACEHOLDER\0Information\0"
		"PLACEHOLDER\0Hide all\0Reveal all\0Sort Matches\0"
		    "SEPARATOR\0Remove\0";
	else
	    jdata->get_ops.ops = "Use for 'Next'\0Reset 'Next'\0Information\0"
		"Configure\0Hide all\0Reveal all\0Sort Matches\0"
		    "SEPARATOR\0Remove\0";
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
	    csmatch_info((mobj_find_oligo *)r, "Find oligo");
	    break;
	case 3: /* Configure */
	    csmatch_configure(io, cs->window, (mobj_find_oligo *)r);
	    break;
	case 4: /* Hide all */
	    csmatch_hide(GetInterp(), cs->window, (mobj_find_oligo *)r,
			 csplot_hash);
	    break;
	case 5: /* Reveal all */
	    csmatch_reveal(GetInterp(), cs->window, (mobj_find_oligo *)r,
			   csplot_hash);
	    break;
	case 6: /* Sort */
	    qsort(r->match, r->num_match, sizeof(obj_match), sort_func);
	    csmatch_reset_hash(csplot_hash, (mobj_repeat *)r);
	    r->current = -1;
	    break;
	case 7: /* Remove */
	    csmatch_remove(io, cs->window,
			   (mobj_find_oligo *)r,
			   csplot_hash);
	    break;
	}
	break;


    case REG_PARAMS:

	jdata->params.string = r->params;
	break;


    case REG_NUMBER_CHANGE:

	csmatch_renumber(io, contig, jdata->number.number,
			 (mobj_find_oligo *)r, csplot_hash, cs->window);
	break;


    case REG_ORDER:

	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;

    case REG_QUIT:

	csmatch_remove(io, cs->window,
		       (mobj_find_oligo *)r,
		       csplot_hash);
	break;

    case REG_DELETE:

	csmatch_contig_delete(io, (mobj_find_oligo *)r, contig,
			      cs->window, csplot_hash);
	break;

    case REG_LENGTH:
	csmatch_replot(io, (mobj_find_oligo *)r, csplot_hash, cs->window);
	break;
    }
}

int
RegFindOligo(GapIO *io,
	     int type,
	     int *pos1,
	     int *pos2,
	     int *score,
	     int *length,
	     int *c1,
	     int *c2,
	     int n_matches)
{
    mobj_find_oligo *find_oligo;
    obj_match *matches = NULL;
    char *val;
    int i, id;

    if (0 == n_matches)
	return 0;

    if (NULL == (find_oligo = (mobj_find_oligo *)xmalloc(sizeof(mobj_find_oligo))))
	return -1;

    if (NULL == (matches = (obj_match *)xmalloc(n_matches * sizeof(obj_match))))
	return -1;

    find_oligo->num_match = n_matches;
    find_oligo->match = matches;
    find_oligo->io = io;
    strcpy(find_oligo->tagname, CPtr2Tcl(find_oligo));

    val = get_default_string(GetInterp(), gap_defs, "FINDOLIGO.COLOUR");
    strcpy(find_oligo->colour, val);

    find_oligo->linewidth = get_default_int(GetInterp(), gap_defs,
					    "FINDOLIGO.LINEWIDTH");

    find_oligo->params = (char *)xmalloc(100);
    if (find_oligo->params)
	sprintf(find_oligo->params, "Unknown at present");
    find_oligo->all_hidden = 0;
    find_oligo->current = -1;
    find_oligo->reg_func = find_oligo_callback;
    find_oligo->match_type = REG_TYPE_OLIGO;

    for (i = 0; i < n_matches; i++) {
	if (type == TAG) {
	    matches[i].func =
		(void *(*)(int, void *, struct obj_match_t *,
			   struct mobj_repeat_t *))find_oligo_obj_func1;
	} else if (type == SEQUENCE) {
	    matches[i].func =
		(void *(*)(int, void *, struct obj_match_t *,
			   struct mobj_repeat_t *))find_oligo_obj_func2;
	} else {
	    return -1;
	}
	matches[i].data = find_oligo;
	matches[i].c1 = c1[i];
	matches[i].c2 = c2[i];
	matches[i].pos1 = pos1[i];
	matches[i].pos2 = pos2[i];
	matches[i].length = length[i];
	matches[i].score = score[i];
	matches[i].flags = 0;
    }

    /* Sort matches */
    qsort(find_oligo->match, find_oligo->num_match, sizeof(obj_match),
	  sort_func);

    PlotRepeats(io, find_oligo);
    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(find_oligo), NULL);

    /*
     * Register the find oligo search with each of the contigs used.
     * Currently we assume that this is all.
     */
    id = register_id();
    for (i = 1; i <= NumContigs(io); i++) {
	contig_register(io, i, find_oligo_callback, (void *)find_oligo, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_NUMBER_CHANGE | REG_ORDER, REG_TYPE_OLIGO);
    }
    return 0;
}

/*
 * return the consensus sequence on contig c_num at position 'position'
 * and length 'length'
 *
 */
static char *
GetTagSequence(GapIO *io,                                             /* in */
	       int c_num,                                             /* in */
	       int position,                                          /* in */
	       int length)                                            /* in */
{
    char *sequence;
    static char seq[1024];

    if (length < 1024)
	sequence = seq;
    else
	if (NULL == (sequence = (char *)xmalloc((length + 1) * sizeof(char ))))
	    return NULL;

    calc_consensus(c_num, position, position+length-1, CON_SUM, sequence,
		   NULL, NULL, NULL, consensus_cutoff, quality_cutoff,
		   database_info, io);

    sequence[length] = '\0';
    return sequence;
}

static void DelTagSequence(char *sequence, int length) {
    if (length >= 1024)
	xfree(sequence);
}

int inexact_pad_match(char *seq,
		      int seq_len,
		      char *string,
		      int string_len,
		      int mis_match,
		      int *match,
		      int *score,
		      int max_matches)
{
    char *pos;
    char *uppert;
    int i;
    int n_matches;
    int n_mis;

    /* remove any pads from the pattern search */
    depad_seq(string, &string_len, NULL);

    /* uppercase search string */
    if (NULL == (uppert = (char *)xmalloc(string_len + 1)))
	return -2;
    uppert[string_len] = 0;
    for (i = string_len-1; i >= 0; i--) {
	uppert[i] = toupper(string[i]);
    }
    for (i = 0; i < seq_len; i++) {
	seq[i] = toupper(seq[i]);
    }
    pos = NULL;

    n_matches = 0;
    pos = pstrstr_inexact(seq, uppert, mis_match, &n_mis);
    while (pos) {
	if (n_matches < max_matches) {
	    match[n_matches] = pos - seq;
	    score[n_matches] = string_len - n_mis;
	    n_matches++;
	} else {
	    /* make positions start at 1 */
	    for (i=0; i < max_matches; i++) {
		match[i]++;
	    }
	    return -1; /* out of match storage */
	}
	pos++;
	pos = pstrstr_inexact(pos, uppert, mis_match, &n_mis);
    }
    /* make positions start at 1 */
    for (i=0; i < n_matches; i++) {
	match[i]++;
    }
    xfree(uppert);
    return n_matches;
}

/*
 * find matches between string covered by a tag and contig list with
 * a minimum match of mis_match
 */
int
TagMatch(GapIO *io,                                                    /* in */
	 int max_clen,                                                 /* in */
	 int num_contigs,                                              /* in */
	 contig_list_t *contig_array,                                  /* in */
	 char **cons_array,                                            /* in */
	 float mis_fmatch,                                             /* in */
	 int *pos1,                                                   /* out */
	 int *pos2,                                                   /* out */
	 int *score,                                                  /* out */
	 int *length,                                                 /* out */
	 int *c1,                                                     /* out */
	 int *c2,                                                     /* out */
	 int max_matches)                                              /* in */
{
    extern char **active_tag_types;
    extern int number_of_active_tags;
    char *string;
    char *cons_match;
    int *match;
    int *scores;
    int mis_match;
    GAnnotations *annotation;
    int i, k, l, c;
    int n_matches = 0;
    int seq_len;
    int cnt = 0;
    char title[1024];
    char name1[10];
    char name2[10];
    int res, too_many = 0;
    int max_imatches = max_matches;

    if (NULL == (scores = (int *)xmalloc((max_matches) * sizeof(int ))))
	return -1;
    if (NULL == (match = (int *)xmalloc((max_matches) * sizeof(int ))))
	return -1;
    if (NULL == (cons_match = (char *)xmalloc((max_clen + 1) * sizeof(char ))))
	return -1;

    /* find tags on each contig in contig_array */
    for (i = 0; i < num_contigs; i++) {

	seq_len = strlen(cons_array[i]);

	/* find tags for current contig */
	annotation = vtagget(io, -contig_array[i].contig,
			     number_of_active_tags, active_tag_types);
	while (annotation && annotation != (GAnnotations *)-1){

	    string = GetTagSequence(io, contig_array[i].contig,
				    annotation->position,
				    annotation->length);

	    /* convert percentage mis-matches into number of mis matches */
	    mis_match = strlen(string) - (ceil(strlen(string) * mis_fmatch / 100.));

	    /* complement string */
	    for (c = 0; c < 2; c++) {
		if (c == 1)
		    complement_seq(string, strlen(string));

		/* loop through list of contigs looking for string */
		for (k = 0; k < num_contigs; k++) {
		    seq_len = strlen(cons_array[k]);
		    n_matches = inexact_pad_match(cons_array[k], seq_len,
						  string,
						  strlen(string),
						  mis_match,
						  match,
						  scores,
						  max_imatches);
		    if (n_matches == -1) {
			verror(ERR_WARN, "find_oligos", "Too many matches");
			n_matches = max_imatches;
		    }

		    /* remove self matches */
		    for (l = 0; l < n_matches; l++) {
			if ((contig_array[i].contig !=
			     contig_array[k].contig) ||
			    (annotation->position !=
			     (match[l] + contig_array[i].start - 1))) {

			    length[cnt] = strlen(string);
			    if (c == 0) {
				c1[cnt] = contig_array[i].contig;
			    } else {
				c1[cnt] = -contig_array[i].contig;
			    }
			    c2[cnt] = contig_array[k].contig;
			    pos1[cnt] = annotation->position;
			    pos2[cnt] = match[l] + contig_array[i].start - 1;
			    score[cnt] = scores[l];
			    /*
			       printf("c1 %d c2 %d pos1 %d pos2 %d \n",
			       c1[cnt], c2[cnt], pos1[cnt], pos2[cnt]);
			       */
			    strncpy(cons_match,&cons_array[k][pos2[cnt]-1],
				    length[cnt]);
			    cons_match[length[cnt]] = '\0';

			    sprintf(title, "Match found between tag on contig "
				    "%d in the %c sense and contig %d",
				    io_clnbr(io, ABS(c1[cnt])),
				    c1[cnt] > 0 ? '+' : '-',
				    io_clnbr(io, c2[cnt]));

			    sprintf(name1, "%d", io_clnbr(io, ABS(c1[cnt])));
			    sprintf(name2, "%d", io_clnbr(io, ABS(c2[cnt])));
			    res = list_alignment(string, cons_match,
						 name1, name2, pos1[cnt],
						 pos2[cnt], title);
			    cnt++;
			    max_imatches--;
			}
		    }

		    if (max_imatches <= 0) {
			too_many = 1;
			break;
		    }
		}

		if (too_many)
		    break;
	    }
	    DelTagSequence(string, annotation->length);

	    if (too_many)
		break;

	    annotation = vtagget(io, 0, number_of_active_tags,
				 active_tag_types);
	} /* end while */

	if (too_many)
	    break;

    } /* end for */

    vmessage("Number of matches found %d \n", cnt);
    xfree(cons_match);
    xfree(match);
    xfree(scores);
    return cnt;
}

/*
 * depad seq until it is of seq_len.
 * Returns the original padded length.
 */
int depad_seq_len(char *strto, char *strfrom, int seq_len)
{
    int len = 0;
    char *a = strto;
    char *b = strfrom;

    while (len < seq_len) {
	if (*b != '*') {
	    *a++ = *b++;
	    len++;
	} else {
	    b++;
	}
    }
    *a = '\0';

    return b-strfrom;
}

/*
 * find matches between user entered sequence string and contig list with
 * a minimum match of mis_match
 */
int
StringMatch(GapIO *io,                                                 /* in */
	    int num_contigs,                                           /* in */
	    contig_list_t *contig_array,                               /* in */
	    char **cons_array,                                         /* in */
	    char *string,                                              /* in */
	    float mis_fmatch,                                          /* in */
	    int *pos1,                                                /* out */
	    int *pos2,                                                /* out */
	    int *score,                                               /* out */
	    int *length,                                              /* out */
	    int *c1,                                                  /* out */
	    int *c2,                                                  /* out */
	    int max_matches,                                           /* in */
	    int consensus_only,                                        /* in */
	    int cutoff_data)					       /* in */
{
    int n_matches = 0;
    int i, j, k, c;
    int mis_match;
    int seq_len;
    int orig;
    int res, too_many = 0;
    char *cons_match;
    char title[1024];
    char name1[10];
    int max_imatches = max_matches;
    size_t stringlen = strlen(string);

    if (NULL == (cons_match = (char *)xmalloc(stringlen + 1)))
	return -1;

    /* convert percentage mis-matches into number of mis matches */
    mis_match = strlen(string) - (ceil(strlen(string) * mis_fmatch / 100.));

    /* complement string */
    for (c = 0; c < 2; c++) {
	int rn = 0;
	if (c == 1)
	    complement_seq(string, stringlen);

	for (i = 0; i < num_contigs; i++) {
	    char *seq_alloc;

	    /*
	     * Consensus first time through loop.
	     * Sequences in that contig on subsequent loops.
	     */
	    do {
		char *seq;
		GReadings r;
		
		seq_alloc = NULL;
		if (rn == 0) {
		    seq = cons_array[i];
		    seq_len = strlen(cons_array[i]);
		} else {
		    gel_read(io, rn, r);

		    /* Basic bounds checking to optimise */
		    if (r.position > contig_array[i].end)
			break;
		    if (r.position + r.sequence_length < contig_array[i].start)
			goto next; /* Sorry! */

		    seq_alloc = (char *)TextAllocRead(io, r.sequence);
		    if (cutoff_data) {
			seq = seq_alloc;
			seq_len = r.length;
		    } else {
			seq = seq_alloc + r.start;
			seq_alloc[r.end-1] = 0;
			seq_len = r.sequence_length;
		    }
		}

		orig = n_matches;
		res = inexact_pad_match(seq, seq_len, string,
					stringlen, mis_match,
					&pos1[n_matches], &score[n_matches],
					max_imatches);
		if (res == -2)
		    return -1;

		if (res < 0) {
		    verror(ERR_WARN, "find_oligos", "Too many matches");
		    too_many = 1;
		    res = max_imatches;
		}
		n_matches += res;
		max_imatches -= res;

		for (j = k = orig; j < n_matches; j++) {
		    int padded_len;

		    c1[j] = contig_array[i].contig;
		    if (c == 0) {
			c2[j] = contig_array[i].contig;
		    } else {
			c2[j] = -contig_array[i].contig;
		    }

		    /*
		     * remove pads such that the final length of cons_match is
		     * of length length[j]
		     */
		    padded_len =
			depad_seq_len(cons_match, &seq[pos1[j]-1], stringlen);

		    if (rn) {
			if (cutoff_data) {
			    pos1[j] += r.position-1 - r.start;
			} else {
			    pos1[j] += r.position-1;
			}
		    }
		    length[j] = padded_len;

		    /* Adjust for searching in a sub-range of the contig */
		    if (rn == 0)
			pos1[j] += contig_array[i].start-1;
		    pos2[j] = pos1[j];

		    /*
		     * The searching above may find hits outside of
		     * contig_array[i].start and contig_array[i].end.
		     *
		     * This happens if we search sequences and the
		     * sequence overlaps the desired range, but has a
		     * hit outside of the desired range.
		     *
		     * Rather than complicate the above code, we post
		     * filter these false hits here.
		     */
		    if (pos1[j] >= contig_array[i].start &&
			pos1[j] <= contig_array[i].end) {
			sprintf(name1, "%d", io_clnbr(io, ABS(c1[j])));
			sprintf(title, "Match found with contig %d read #%d "
				"in the %c sense",
				io_clnbr(io, ABS(c2[j])), rn, 
				c2[j] > 0 ? '+' : '-');

			list_alignment(string, cons_match, "oligo", name1, 1,
				       pos1[j], title);
			
			/*
			 * Copy it from *[j] to *[k].
			 * This code REALLY needs to be using structs!
			 * This is foul.
			 */
			pos1  [k] = pos1  [j];
			pos2  [k] = pos2  [j];
			c1    [k] = c1    [j];
			c2    [k] = c2    [j];
			length[k] = length[j];
			score [k] = score [j];
			k++;
		    }
		}

		n_matches -= j-k;
		max_imatches += j-k;

		if (seq_alloc)
		    xfree(seq_alloc);

		if (too_many)
		    break;

	    next:
		if (consensus_only) {
		    rn = 0;
		} else {
		    rn = rn
			? io_rnbr(io, rn)
			: io_clnbr(io, contig_array[i].contig);
		}
	    } while(rn);

	    if (too_many)
		break;
	}

	if (too_many)
	    break;
    }

    xfree(cons_match);
    vmessage("Number of matches found %d \n", n_matches);
    return n_matches;
}

int
find_oligos(GapIO *io,
	    int num_contigs,
	    contig_list_t *contig_array,
	    float mis_match,
	    char *string,
	    int consensus_only,
	    int in_cutoff)
{
    int i;
    int *pos1 = NULL;
    int *pos2 = NULL;
    int *score = NULL;
    int *length = NULL;
    int *c1 = NULL;
    int *c2 = NULL;
    int max_matches, abs_max;
    int seq_len;
    int n_matches;
    int max_clen;
    char **cons_array = NULL;

    /* Calculate maximum contig length and total contig length */
    max_clen = io_clength(io, CalcLongContig(io));
    for (max_matches = 0, max_clen = 0, i=0; i<num_contigs; i++) {
	if (io_clength(io, contig_array[i].contig) > max_clen)
	    max_clen = io_clength(io, contig_array[i].contig);
	max_matches += io_clength(io, contig_array[i].contig);
    }
    max_matches *= 2; /* both strands */

    abs_max = get_default_int(GetInterp(), gap_defs, "FINDOLIGO.MAX_MATCHES");

    if (max_matches > abs_max)
	max_matches = abs_max;

    if (NULL == (pos1 = (int *)xmalloc((max_matches + 1) * sizeof(int ))))
	goto error;
    if (NULL == (pos2 = (int *)xmalloc((max_matches + 1) * sizeof(int ))))
	goto error;
    if (NULL == (score = (int *)xmalloc((max_matches + 1) * sizeof(int ))))
	goto error;
    if (NULL == (length = (int *)xmalloc((max_matches + 1) * sizeof(int ))))
	goto error;
    if (NULL == (c1 = (int *)xmalloc((max_matches + 1) * sizeof(int ))))
	goto error;
    if (NULL == (c2 = (int *)xmalloc((max_matches + 1) * sizeof(int ))))
	goto error;

    /* save consensus for each contig */
    if (NULL == (cons_array = (char **)xmalloc(num_contigs * sizeof(char *))))
	goto error;

    for (i = 0; i < num_contigs; i++) {
	seq_len = contig_array[i].end - contig_array[i].start + 1;
	if (NULL == (cons_array[i] = (char *)xmalloc(seq_len + 1)))
	    goto error;

	calc_consensus(contig_array[i].contig,
		       contig_array[i].start, contig_array[i].end,
		       CON_SUM, cons_array[i],
		       NULL, NULL, NULL, consensus_cutoff, quality_cutoff,
		       database_info, io);
	cons_array[i][seq_len] = '\0';
    }

    /* do match on either tag(s) or string */
    if (string && *string) {
	n_matches = StringMatch(io, num_contigs, contig_array,
				cons_array, string, mis_match, pos1, pos2,
				score, length, c1, c2, max_matches,
				consensus_only, in_cutoff);
	if (-1 == RegFindOligo(io, SEQUENCE, pos1, pos2, score, length, c1,
			       c2, n_matches))
	    goto error;
    } else {
	if (-1 == (n_matches = TagMatch(io, max_clen, num_contigs,
					contig_array, cons_array,
					mis_match, pos1, pos2,
					score, length, c1, c2, max_matches)))
	    goto error;
	if (-1 == RegFindOligo(io, TAG, pos1, pos2, score, length, c1, c2,
			       n_matches))
	    goto error;
    }

    for (i = 0; i < num_contigs; i++) {
	if (cons_array[i])
	    xfree(cons_array[i]);
    }
    xfree(cons_array);
    xfree(c1);
    xfree(c2);
    xfree(pos1);
    xfree(pos2);
    xfree(score);
    xfree(length);
    return 0;

 error:
    if (c1)
	xfree(c1);
    if (c2)
	xfree(c2);
    if (cons_array)
	xfree(cons_array);
    if (pos1)
	xfree(pos1);
    if (pos2)
	xfree(pos2);
    if (score)
	xfree(score);
    if (length)
	xfree(length);

    return -1;
}


/*
 * Wrapper around find_oligos() to look for all oligos listed in a FASTA file.
 */
int
find_oligo_file(GapIO *io,
		int num_contigs,
		contig_list_t *contig_array,
		float mis_match,
		char *file,
		int consensus_only,
		int in_cutoff)
{
    char **ids;
    int nids;
    char line[8192];
    int i;
    int r = 0; /* ret. code */

    /* Use seq_utils to parse the input file */
    if (0 != get_identifiers(file, &ids, &nids))
	return -1;

    for (i = 0; i < nids; i++) {
	char *seq;
	int seq_len;

	seq = NULL;
	seq_len = 0;

	if (0 != get_seq(&seq, maxseq, &seq_len, file, ids[i]))
	    continue;

	if (seq_len == 0 || !seq || *seq == 0) {
	    if (seq)
		xfree(seq);
	    continue;
	}

	vmessage("Sequence search for ID '%s'\n", ids[i]);
	r |= find_oligos(io, num_contigs, contig_array, mis_match,
			 seq, consensus_only, in_cutoff);
	vmessage("\n");

	if (seq)
	    xfree(seq);
    }

    /* Tidy up memory */
    for (i = 0; i < nids; i++) {
	xfree(ids[i]);
    }
    xfree(ids);

    return r;
}
