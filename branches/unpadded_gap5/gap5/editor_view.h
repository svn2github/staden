#ifndef _EDITOR_VIEW_H_
#define _EDITOR_VIEW_H_

#include "tg_gio.h"
#include "tkSheet.h"
#include "tkEditor.h"
#include "tkEdNames.h"
#include "consensus.h"

#define MAX_DISPLAY_WIDTH 1000
#define MAX_NAME_WIDTH    256
#define WIN_NAME_SIZE 100

/*
 * Display update flags - see xx->refresh_flags
 */
#define ED_DISP_NAMES		(1<<0)
#define ED_DISP_READS		(1<<1)
#define ED_DISP_CONS		(1<<2)
#define ED_DISP_STATUS		(1<<4)
#define ED_DISP_RULER		(1<<5)
#define ED_DISP_CURSOR		(1<<6)
#define ED_DISP_SELECTION	(1<<7)
#define ED_DISP_HEIGHT		(1<<8)
#define ED_DISP_XSCROLL		(1<<3)
#define ED_DISP_YSCROLL		(1<<9)
/* Sequence and consensus */
#define ED_DISP_SEQS		(ED_DISP_READS | ED_DISP_CONS)
/* Colouring changes - quality, cutoffs, etc */
#define ED_DISP_COLOUR		(1<<13)
/* Everything */
#define ED_DISP_ALL		(ED_DISP_SEQS | ED_DISP_NAMES | \
				 ED_DISP_XSCROLL | ED_DISP_YSCROLL | \
                                 ED_DISP_STATUS | \
				 ED_DISP_RULER | ED_DISP_CURSOR | \
				 ED_DISP_SELECTION | ED_DISP_HEIGHT | \
				 ED_DISP_COLOUR)
/* Single items only - see xx->refresh_seq */
#define ED_DISP_NAME		(1<<10)
#define ED_DISP_SEQ		(1<<11)
/* For optimising redisplayDisagreements to remove flicker */
#define ED_DISP_NO_DIFFS	(1<<12)

/* Orig positions not implemented yet */
#define origpos(xx, seq, pos) ((pos)+1)

/* Flags not implemented yet */
#define DB_Flags(xx, seq) (0)
#define DB_setFlags(xx, seq, flags) (0)

struct _EdLink;
enum States {StateDown=0,StateUp};

typedef struct {
    struct _edview *xx;
    int ncols;          // count of columns filled so far
    int nrows;          // maximum number of rows
    int *pos;           // X coordinates, one per ncols
    int *nth;           // X coordinates, one per ncols
    char **str;         // base strings, nrows*ncols
    XawSheetInk **ink;  // ink strings, nrows*ncols
    int *rec;           // record numbers (index into xx->r), nrows
    int *rpos;          // +ve => starts that many bases in from left edge.
                        // -ve => -1st visible base, size xx->nr
    int *cigar_ind;	// Index into s->alignment CIGAR
    int *cigar_len;	// Number of bases left for op at cigar_ind
} edscreen_t;

typedef struct _edview {
    /* A derived IO struct */
    GapIO *io;
    tg_rec cnum;
    contig_t *contig;

    /* Necessary Tcl/Tk bits and bobs */
    Tcl_Interp *interp;
    int editor_id;
    char edname[20];
    char seq_win[WIN_NAME_SIZE];
    char name_win[WIN_NAME_SIZE];
    Editor *ed;
    edNames *names;

    /* Display window size and location */
    int displayPos;
    int displayYPos;
    int displayWidth;
    int displayHeight;
    int names_xPos;

    void (*dispFunc)(void *, int, int, int, void *);
    enum States editorState;
    int refresh_flags;
    tg_rec refresh_seq;

    /* Cached data */
    char        displayedConsensus[MAX_DISPLAY_WIDTH];
    consensus_t cachedConsensus[MAX_DISPLAY_WIDTH];

    /* Cursor coordinates, sequence number and position in seq */
    int    cursor_type;
    tg_rec cursor_rec;
    int    cursor_pos;
    int    cursor_nth;
    int    cursor_apos; /* absolute position in contig */
    int    cursor_anth; /* nth base at abs pos */
    cursor_t *cursor;
    int reg_id; /* registration id */

    /* Y coordinates for the elements in the editor */
    int y_seq_start;
    int y_seq_end;
    int y_cons;
    int y_numbers;

    /* Join editor bits */
    struct _EdLink *link;

    int trace_lock;

    /* Cached range query & Y position */
    rangec_t *r;
    int nr;
    int max_height;
    int r_start, r_end;
    /* FIXME: add cached index into r[] foreach row[y], as it's sorted on y */

    /* Cached screen data, computed when X changes */
    edscreen_t *screen;

    /* Maps r[i].anno.obj_rec to i */
    HacheTable *anno_hash;

    /* Maps record numbers to r[i] indices */
    HacheTable *rec_hash;

    /* Selection */
    int    select_made;
    tg_rec select_seq;
    int    select_start, select_start_nth;
    int    select_end,   select_end_nth;
} edview;

typedef struct _EdLink {
    edview *xx[2];		/* A pair of linked editors */
    int locked;			/* Set if these are locked together */
    int lockOffset;		/* And their relative distances */
    tkSheet *diffs;
} EdLink;

/*
 * A C interface to the edit_contig and join_contig Tcl functions.
 */
int edit_contig(GapIO *io, tg_rec cnum, tg_rec rnum, int pos);
int join_contig(GapIO *io, tg_rec cnum[2], tg_rec rnum[2], int pos[2]);

/*
 * Allocates and initialises a new edview
 */
edview *edview_new(GapIO *io, tg_rec contig, tg_rec crec, int cpos,
		   Editor *ed, edNames *names,
		   void (*dispFunc)(void *, int, int, int, void *),
		   Tcl_Interp *interp);
/*
 * Deallocates an edview
 */
void edview_destroy(edview *xx);

/*
 * Deallocates an edview
 */
int edview_redraw(edview *xx);

/* Update X scrollbar of names display */
void ed_set_nslider_pos(edview *xx, int pos);
int set_displayPos(edview *xx, int pos);

/*
 * Formats tag information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 */
char *edGetBriefTag(edview *xx, tg_rec anno_ele, char *format);

/*
 * Formats reading information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 */
char *edGetBriefSeq(edview *xx, tg_rec seq, int pos, int nth, char *format);

/*
 * Formats consensus information for the status line.
 * This is done using a format string where certain % rules are replaced by
 * appropriate components.
 */
char *edGetBriefCon(edview *xx, tg_rec crec, int pos, int nth, char *format);

/*
 * Given an X,Y coordinate return the reading id under this position.
 *
 * Returns record number on success
 *         -1 on failure.
 */
int edGetGelNumber(edview *xx, int x, int y);

/*
 * Identifies the type of object underneath a specific row and column.
 * 'name' is a boolean which when true indicates the row,col are in the
 * names panel instead of the sequence panel.
 * 'seq_only' forces the item to be a sequence or consensus, and not
 * an object on them (eg annotation).
 *
 * Returns the item type GT_* on success and the record/pos in *rec, *pos
 *         -1 on failure (eg numbers, off screen, etc)
 */
int edview_item_at_pos(edview *xx, int row, int col, int name, int exact,
		       int seq_only, tg_rec *rec, int *pos, int *nth);

/* Cursor movement control */
void edSetApos(edview *xx);
int edSetCursorPos(edview *xx, int type, tg_rec rec, int pos, int nth,
		   int visible);
int edCursorUp(edview *xx);
int edCursorDown(edview *xx);
int edCursorLeft(edview *xx);
int edCursorRight(edview *xx);
int edReadStart(edview *xx);
int edReadStart2(edview *xx);
int edReadEnd(edview *xx);
int edReadEnd2(edview *xx);
int edContigStart(edview *xx);
int edContigEnd(edview *xx);

/*
 * Force the cursor to be visible. If x_safe or y_safe are true then
 * we omit some of the searching and assume there is no reason to check
 * that x or y is still visible.
 *
 * Returns 1 if redraw has taken place
 *         0 if not
 */
int showCursor(edview *xx, int x_safe, int y_safe);


int inJoinMode(edview *xx);
void edDisplayTrace(edview *xx);


/*
 * Given a sequence record number this identifies all other sequence
 * records from the same template. The returned array is malloced and should
 * be freed by the caller once finished with.
 *
 * Returns pointer to array of records of size *nrec on success
 *         NULL on failure (or zero found)
 */
tg_rec *edGetTemplateReads(edview *xx, tg_rec seqrec, int *nrec);

/*
 * Handles the align button in the join editor
 * Returns 0 for success
 *        -1 for failure
 */
int edJoinAlign(edview *xx, int fixed_left, int fixed_right);

/*
 * Perform the actual join process
 * Returns 0 for success
 *        -1 for failure
 */
int edJoin(edview *xx);

/*
 * Populates the cache of visible items in xx->r and xx->nr.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int edview_visible_items(edview *xx, int start, int end);

/*
 * X11 Selection handling
 */
int edSelectClear(edview *xx);
void edSelectFrom(edview *xx, int pos, int nth);
void edSelectTo(edview *xx, int pos, int nth);
void edSelectSet(edview *xx, tg_rec rec,
		 int start, int start_nth,
		 int end, int end_nth);

/*
 * Searching - see edview_search.c
 */
int edview_search(edview *xx, int forwards, int strand,
		  char *type, char *value);

/*
 * Convert from a record and position to a window X,Y coordinate in sheet
 * units.
 *
 * Returns 0 on success and stores via x and y pointers.
 *        -1 on failure (rec/pos not visible).
 */
int edGetXY(edview *xx, int rec_type, tg_rec rec, int pos, int *x, int *y);


/* In editor_join.c */
int edJoinMismatch(edview *xx, int *len, int *mismatch);

#endif /* _EDITOR_VIEW_H_ */
