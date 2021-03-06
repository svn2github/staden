/*
 *
 * gap_range.c - a tcl command to hold the x range for tracks.  For example
 *           template_display and depth_display
 *
 */




#include <string.h>

#include "gap_range.h"
#include "gap_cli_arg.h"


/*** Tcl command functions ***/ 

static int grange_cmd(ClientData clientData, Tcl_Interp *interp,
			int objc, Tcl_Obj *CONST objv[]) {
    // don't really want it to do anything
    return TCL_OK;
}

static void grange_cmd_delete(ClientData clientData) {
    gap_range_destroy((gap_range_t *)clientData);
    ckfree(clientData);
}
			

typedef struct {
     GapIO *io;
     int cnum;
} grange_arg;

/* forward declaration */
static void update_filter(gap_range_t *gr);

static int tcl_range(ClientData clientData, Tcl_Interp *interp,
			int objc, Tcl_Obj *CONST objv[]) {
    gap_range_t *gr;
    char name[30];
    
    grange_arg args;
    cli_args a[] = {
	{"-io",     ARG_IO,  1, NULL, offsetof(grange_arg, io)},
	{"-cnum",   ARG_INT, 1, NULL, offsetof(grange_arg, cnum)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;
	
    if (NULL == (gr = (gap_range_t *)ckalloc(sizeof(gap_range_t)))) {
    	return TCL_ERROR;
    }
    
    gr->io = args.io;
    gr->contig = (contig_t *)cache_search(gr->io, GT_Contig, args.cnum);
    gr->crec = gr->contig->rec;
    
    if (NULL == gr->contig) {
    	return TCL_ERROR;
    }
    
    cache_incr(gr->io, gr->contig);
    
    gr->r             = NULL;
    gr->tl            = NULL;
    gr->depth         = NULL;
    gr->wx0           = -1;
    gr->wx1           = -1;
    gr->nr            = 0;
    gr->template_mode = -1;
    gr->width         = -1;
    gr->ntl           = -1;
    set_filter(gr, -1, -1, -1, -1, -1);
    update_filter(gr);
    
    sprintf(name, "grange=%p", gr);
    
    if (NULL == Tcl_CreateObjCommand(interp, name, grange_cmd,
    					(ClientData)gr,
					(Tcl_CmdDeleteProc *)grange_cmd_delete)) {
	return TCL_ERROR;
    }
    
    Tcl_SetObjResult(interp, Tcl_NewStringObj(name, -1));
    
    return TCL_OK;
}
    
    
int GRange_Init(Tcl_Interp *interp) {
    if (NULL == Tcl_CreateObjCommand(interp, "g5::range", tcl_range,
					(ClientData)NULL,
					(Tcl_CmdDeleteProc *)NULL)) {
	return TCL_ERROR;
    }
    
    return TCL_OK;
}

/*** range functions ***/

void gap_range_destroy(gap_range_t *gr) {
    if (gr->r) {
    	free(gr->r);
    }
    
    if (gr->tl) {
    	free(gr->tl);
    }
    
    if (gr->depth) {
    	free(gr->depth);
    }
}

int gap_range_test(gap_range_t *gr) {
//    printf("GR - TEST\n");
    printf("r %p wx0 %f wx1 %f nr %d\n", gr->r, gr->wx0, gr->wx1, gr->nr);
    
    return 1;
}




static int is_filter_change(gap_range_t *gr) {
    int changed = 0;
    
    if (gr->old_filter.filter   != gr->new_filter.filter ||
    	gr->old_filter.min_qual != gr->new_filter.min_qual ||
	gr->old_filter.max_qual != gr->new_filter.max_qual ||
	gr->old_filter.c_mode   != gr->new_filter.c_mode) {
	    
    	changed = 1;
    }
    
    return changed;
}

void set_filter(gap_range_t *gr, int filter, int min, int max, int mode, int accuracy) {

    gr->new_filter.filter   = filter;
    gr->new_filter.min_qual = min;
    gr->new_filter.max_qual = max;
    gr->new_filter.c_mode   = mode;
    gr->new_filter.accuracy = accuracy;
}


static void update_filter(gap_range_t *gr) {

    gr->old_filter.filter   = gr->new_filter.filter;
    gr->old_filter.min_qual = gr->new_filter.min_qual;
    gr->old_filter.max_qual = gr->new_filter.max_qual;
    gr->old_filter.c_mode   = gr->new_filter.c_mode;
    gr->old_filter.accuracy = gr->new_filter.accuracy;
}
    


/*
   recalculate the sequence range if neccessary
   if out of memory, returns 1 with a null gr->r pointer
   returns 1 on change, 0 otherwise
*/
int gap_range_recalculate(gap_range_t *gr, int width, double new_wx0, double new_wx1, int new_mode, int force) {
    int changed = 0;

    /* only change if needed or if force is true */
    if (force || gr->r == NULL || new_wx0 != gr->wx0 || new_wx1 != gr->wx1 || new_mode != gr->template_mode ||
    	gr->width != width || gr->new_filter.accuracy != gr->old_filter.accuracy) {

    	if (gr->r) {
	    free(gr->r);
	}
	
	gr->r = contig_seqs_in_range(gr->io, &gr->contig, new_wx0, new_wx1, new_mode, &gr->nr);
	
	if (gr->r) {
	    tline *tmp_line = NULL;
	
	    gr->wx0  = new_wx0;
	    gr->wx1  = new_wx1;
	    gr->template_mode = new_mode;
	    update_filter(gr);
	    
	    tmp_line = (tline *)realloc(gr->tl, gr->nr * sizeof(tline));
	    
	    if (!tmp_line) {
	    	printf("GR tline memory problem\n");
		gr->tl = tmp_line; // but don't stop yet
	    } else {
		gr->tl = tmp_line;
	    }
	    
	    if (gr->width != width) {
	    	gr->width = width;
		gr->depth = (gap_depth_t *)realloc(gr->depth, width * sizeof(gap_depth_t));
	    }
	}
	
	changed = 1;
    }
    
    return changed;
}

int gap_range_x(gap_range_t *gr, double ax_conv, double bx_conv, 
    	    	int fwd_col, int rev_col, int single_col, int span_col, int inconsistent_col,
		int force, int reads_only) {
    
    if (is_filter_change(gr) || force) {
    	int i, j;
    
    	update_filter(gr);
    	gr->ntl = 0;
    	gr->max_height = 0;
	memset(gr->depth, 0, gr->width * sizeof(gap_depth_t));
	
	for (i = 0; i < gr->nr; i++) {
	    int sta, end;
	    int r_sta, r_end;
	    int span   = 0;
	    int single = 0;
	    int col;
	    double mq;
	    rangec_t *r  = &gr->r[i];
	    tline    *tl = &gr->tl[gr->ntl];
	    
	    sta = r->start;
	    end = r->end;
	    mq  = r->mqual;
	    
	    /* depth plot setting (sequence) */
	    r_sta = (sta - bx_conv) * ax_conv; // world to raster conversion
	    r_end = (end - bx_conv) * ax_conv;
	    
	    if (r_sta < 0)  	   r_sta = 0;
	    if (r_end >= gr->width) r_end = gr->width - 1;
	    
	    for (j = r_sta; j <= r_end; j++) {
	    	gr->depth[j].s++;
 		if (gr->depth[j].s > gr->max_height) gr->max_height = gr->depth[j].s;
		
	    }
	    
	    /* deal with filters */
	    if (gr->new_filter.filter & (FILTER_PAIRED | FILTER_SPANNING) && !r->pair_rec) continue;
		
	    if (gr->new_filter.filter & FILTER_SINGLE && r->pair_rec) continue;

	    tl->t_strand = ((r->flags & GRANGE_FLAG_END_REV) != 0) == ((r->flags & GRANGE_FLAG_COMP1) != 0);
	    tl->rec = r->rec;
	    
	    if (!reads_only && r->pair_rec) {
	    	/* only draw once */
		if (!r->rec) continue;
		
		if (r->pair_ind != -1) {
		    gr->r[r->pair_ind].rec = 0;
		}
		
		/* accurate drawing, this will slow things down */
		if (gr->new_filter.accuracy && !(r->pair_start || r->pair_end)) {
		    /* Pair was off-screen, so get the results */
		    seq_t *s;
		    int crec;
		    
		    sequence_get_position(gr->io, r->pair_rec, &crec, &r->pair_start, &r->pair_end, NULL);
		    s = (seq_t *)cache_search(gr->io, GT_Seq, r->pair_rec);
		    r->pair_mqual = s->mapping_qual;
		    
		    if (gr->crec != crec) span = 1;
		} else if (!gr->new_filter.accuracy && !(r->pair_start || r->pair_end)) {
		    span = 1;
		}
		
		if (gr->new_filter.filter & FILTER_SPANNING && !span) continue;
		
		if (gr->new_filter.filter & FILTER_NONSPANNING && span) continue;

		if (!span) {
		    sta = r->start < r->pair_start ? r->start : r->pair_start;
		    end = r->end > r->pair_end ? r->end : r->pair_end;
		    r->flags |= GRANGE_FLAG_CONTIG;
		    
		    switch (gr->new_filter.c_mode) {
			case 0:
			case 3:
			    mq = (r->mqual + r->pair_mqual) / 2.0;
			    break;
			case 1:
			    mq = r->mqual < r->pair_mqual ? r->mqual : r->pair_mqual;
			    break;
			case 2:
			    mq = r->mqual < r->pair_mqual ? r->pair_mqual : r->mqual;
			    break;
		    }
	    	}
	    } else {
	    	mq = r->mqual;
		if (!reads_only) single = 1;
	    }
	    
	    /* now for some colour settings */
	    mq *= 3;
	    if (mq < 0) mq = 0;
	    if (mq > 255) mq = 255;
	    
	    tl->mq = mq;
	    col = (int)(mq / 8);
	    
	    if (single) col = single_col;
	    if (span) 	col = span_col;
	    
	    /* Check consistency */
	    if (r->pair_rec && r->flags & GRANGE_FLAG_CONTIG) {
		if (!((r->start < r->pair_start &&
		       (r->flags & GRANGE_FLAG_COMP1) == 0 &&
		       (r->flags & GRANGE_FLAG_COMP2) != 0) ||
		    (r->pair_start < r->start &&
		       (r->flags & GRANGE_FLAG_COMP1) != 0 &&
		       (r->flags & GRANGE_FLAG_COMP2) == 0)))
		       
		    col = inconsistent_col;
	    }
	    
	    if (gr->new_filter.filter & FILTER_CONSISTENT && col == inconsistent_col) continue;
	    
	    if (gr->new_filter.filter & FILTER_INCONSISTENT && col != inconsistent_col) continue;
	    
	    /* depth plot template pairs */
	    if (!single && !span && col != inconsistent_col) {
	    	r_sta = (sta - bx_conv) * ax_conv; // world to raster conversion
		r_end = (end - bx_conv) * ax_conv;
		
		if (r_sta < 0) r_sta = 0;
		if (r_end >= gr->width) r_end = gr->width - 1;
		
		for (j = r_sta; j <= r_end; j++) {
		    gr->depth[j].t++;
		    if (gr->depth[j].t > gr->max_height) gr->max_height = gr->depth[j].t;
		}
	    }
	    
	    /* Filter */
	    if (tl->mq < gr->new_filter.min_qual || tl->mq > gr->new_filter.max_qual) continue;
	    
	    /* Generate line data */
	    tl->col[1] = col;
	    tl->x[0]   = sta;
	    /* range in col 0 */
	    tl->x[1]   = sta;
	    /* range in col 1 */
	    tl->x[2]   = end; 
	    /* range in col 3 */
	    tl->x[3]   = end;
	    
	    /* Readings too? */
	    if (gr->new_filter.c_mode == 3) {
	    
		col = (r->flags & GRANGE_FLAG_END_MASK)
		    == GRANGE_FLAG_END_FWD ? fwd_col : rev_col;

		if (r->start == tl->x[0]) {
		    tl->x[1]   = r->end;
		    tl->col[0] = col;
		} else if (r->end == tl->x[3]) {
		    tl->x[2]   = r->start;
		    tl->col[2] = col;
		} else {
		    printf("error, start/end do not match template pos (single)\n");
		    printf("start %d/%d end %d/%d\n", r->start, tl->x[0], r->end, tl->x[3]); 
		}
		
		if (r->pair_rec && (r->pair_start || r->pair_end)) {
		
		    col = (r->flags & GRANGE_FLAG_PEND_MASK)
			== GRANGE_FLAG_PEND_FWD ? fwd_col : rev_col;

		    if (r->pair_start == tl->x[0]) {
			tl->x[1]   = r->pair_end;
			tl->col[0] = col;
		    } else if (r->pair_end == tl->x[3]) {
			tl->x[2]   = r->pair_start;
			tl->col[2] = col;
		    } else {
			printf("error, start/end do not match template pos (pair)\n");
			printf("start %d/%d end %d/%d\n", r->pair_start, tl->x[0], r->pair_end, tl->x[3]); 
		    }
		}
		
	    }

	    gr->ntl++;
	}
    }
    
    return gr->ntl;
}



/*
    reset to starting values, free any memory
    leave the gap and contig values
*/
void gap_range_reset(gap_range_t *gr) {
    
    if (gr->r) {
    	free(gr->r);
    }
    
    gr->r             = NULL;
    gr->wx0           = -1;
    gr->wx1           = -1;
    gr->nr            = 0;
    gr->template_mode = -1;
}










































