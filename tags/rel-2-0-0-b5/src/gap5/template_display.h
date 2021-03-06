#ifndef _TEMPLATE_DISPLAY_H_
#define _TEMPLATE_DISPLAY_H_

#include <tcl.h>
#include <tkRaster.h>
#include "tg_gio.h"
#include "template_draw.h"

typedef struct {
    double y;          /* y coord */
    int col[3];        /* colours */
    int x[4];          /* coordinates */
    int t_strand;      /* template strand */
    int mq;            /* mapping qual */
    int rec;	       /* used for yspread */
} tline;

typedef struct {
    GapIO *io;
    contig_t *contig;
    int crec;
    Tk_Raster *raster;
    Tk_Window tkwin;
    Tcl_Interp *interp;
    Tk_OptionTable optionTable;
    int map_col[256];
    int single_col;
    int span_col;
    int inconsistent_col;
    int fwd_col;
    int rev_col;
    int fwd_col3;
    int rev_col3;
    int xhair_col;
    int background;
    int logy;
    int cmode;
    int ymode;
    int accuracy;
    int spread;
    int reads_only;
    double yzoom;
    double xzoom;
    int sep_by_strand;
    int yoffset;
    int ymin, ymax; /* visible extents of data in Y */
    int filter; /* bitmask */
    int *tdepth; /* paired template depth (consistent only) */
    int *sdepth; /* sequence depth */
    int depth_width;
    int plot_depth;
    double xhair_pos;
    double yhair_pos;
    int min_qual, max_qual; /* Filter parameters */
    int min_sz;  /* For stacking mode */
    image_t *image;
    // remember the old parameters for reuse
    double wx0; // the old x range
    double wx1;
    rangec_t *r; // range data
    int nr;
    int mode;
    tline *tl;
    int ntl;
    int old_filter;
    int old_min_qual, old_max_qual;
    int old_cmode;
    int old_accuracy;
    
} template_disp_t;

/* If bit set we filter our this data type */
#define FILTER_PAIRED        (1<<0)
#define FILTER_SINGLE        (1<<1)
#define FILTER_CONSISTENT    (1<<2)
#define FILTER_INCONSISTENT  (1<<3)
#define FILTER_SPANNING      (1<<4)
#define FILTER_NONSPANNING   (1<<5)

int TDisp_Init(Tcl_Interp *interp);

template_disp_t *template_new(GapIO *io, int cnum,
			      Tcl_Interp *interp,
			      Tk_Raster *raster);
void template_destroy(template_disp_t *t);
int template_replot(template_disp_t *t);

#endif /* _TEMPLATE_DISPLAY_H_ */
