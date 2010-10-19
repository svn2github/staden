#include <xalloc.h>
#include <math.h>
#include <ctype.h>
#include <X11/Xatom.h> /* XA_PRIMARY - included in Tk distribrution */
#include <assert.h>

#include "editor_view.h"
#include "tkSheet.h"
#include "tman_interface.h"
#include "gap_globals.h"
#include "qualIO.h"
#include "qual.h"
#include "tg_gio.h"
#include "misc.h"
#include "consensus.h"
#include "pileup.h"
#include "tagdb.h"
#include "active_tags.h"
#include "io_utils.h"
#include "dna_utils.h"

static void redisplaySelection(edview *xx);

/*
 * A C interface to the edit_contig and join_contig Tcl functions.
 */
int edit_contig(GapIO *io, tg_rec cnum, tg_rec rnum, int pos) {
    char cmd[1024];

    sprintf(cmd, "edit_contig -io %s -contig %"PRIrec
	    " -reading %"PRIrec" -pos %d\n",
	    io_obj_as_string(io), cnum, rnum, pos);
    return Tcl_Eval(GetInterp(), cmd);
}

int join_contig(GapIO *io, tg_rec cnum[2], tg_rec rnum[2], int pos[2]) {
    char cmd[1024];
    int ret;

    sprintf(cmd, "join_contig -io %s -contig %"PRIrec" -reading %"PRIrec
	    " -pos %d -contig2 %"PRIrec" -reading2 %"PRIrec" -pos2 %d",
	    io_obj_as_string(io),
	    cnum[0], rnum[0], pos[0],
	    cnum[1], rnum[1], pos[1]);
    ret = Tcl_Eval(GetInterp(), cmd);
    if (ret != TCL_OK) {
	fprintf(stderr, "%s\n", Tcl_GetStringResult(GetInterp()));
    }
    return ret;
}

/*
 * Allocates and initialises a new edview
 */
edview *edview_new(GapIO *io, tg_rec contig, tg_rec crec, int cpos,
		   Editor *ed, edNames *names,
		   void (*dispFunc)(void *, int, int, int, void *),
		   Tcl_Interp *interp)
{
    edview *xx;
    static int editor_id = 1;
    char *cp;

    xx = (edview *)xcalloc(1, sizeof(*xx));
    if (!xx)
	return NULL;
    
    xx->editor_id = editor_id++;

    xx->interp = interp;
    xx->io = io; /* model */
    xx->cnum = contig;
    xx->contig = (contig_t *)cache_search(io, GT_Contig, xx->cnum);
    cache_incr(xx->io, xx->contig);

    xx->ed = ed;
    xx->displayYPos = 0;
    xx->displayWidth = xx->ed->sw.columns;
    xx->displayHeight = xx->ed->sw.rows;
    xx->dispFunc = dispFunc;
    xx->editorState = StateUp;

    xx->y_cons = 0;
    xx->y_numbers = 1;
    xx->y_seq_start = 2;
    xx->y_seq_end = 0;

    xx->names = names;
    xx->names_xPos = 0;

    xx->cursor_pos  = cpos;
    xx->cursor_nth  = 0;
    xx->cursor_rec  = crec ? crec : contig;
    xx->cursor_type = (xx->cursor_rec == 0 || xx->cursor_rec == contig)
	? GT_Contig : GT_Seq;

    xx->trace_lock = 1;
    
    if (!xx->ed->consensus_at_top) {
	xx->ed->sw.yflip = 1;
	xx->names->sw.yflip = 1;
    }

    xx->r = NULL;
    xx->anno_hash = NULL;
    xx->rec_hash = NULL;

    /* Private cursor */
    cp = Tcl_GetVar2(xx->interp, Tk_PathName(xx->ed->sw.tkwin), "reg",
		     TCL_GLOBAL_ONLY);
    xx->reg_id = cp ? atoi(cp) : 0;
    xx->cursor = create_contig_cursor(io->base, contig, 1, xx->reg_id);
    edSetApos(xx);
    xx->displayPos = xx->cursor_apos;
    
    return xx;
}

/*
 * Deallocates an edview
 */
void edview_destroy(edview *xx) {
    if (xx->cursor)
	delete_contig_cursor(xx->io->base, xx->cnum, xx->cursor->id, 1);

    if (xx->r)
	free(xx->r);

    if (xx->anno_hash)
	HacheTableDestroy(xx->anno_hash, 0);
    
    if (xx->rec_hash)
	HacheTableDestroy(xx->rec_hash, 0);

    cache_decr(xx->io, xx->contig);

    xfree(xx);
}

static seq_t *get_seq(GapIO *io, tg_rec rec) {
    return (seq_t *)cache_search(io, GT_Seq, rec);
}


/* ----- 'brief' line manipulation ----- */

static void add_number(char *buf, int *j, int l1, int l2, int val) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*d", l1, l2, val);
	else
	    *j += sprintf(buf + *j, "%*d", l1, val);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*d", l2, val);
	else
	    *j += sprintf(buf + *j, "%d", val);
}

static void add_number64(char *buf, int *j, int l1, int l2, int64_t val) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*"PRId64, l1, l2, val);
	else
	    *j += sprintf(buf + *j, "%*"PRId64, l1, val);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*"PRId64, l2, val);
	else
	    *j += sprintf(buf + *j, "%"PRId64, val);
}

static void add_double(char *buf, int *j, int l1, int l2, double val) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*f", l1, l2, val);
	else
	    *j += sprintf(buf + *j, "%*f", l1, val);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*f", l2, val);
	else
	    *j += sprintf(buf + *j, "%f", val);
}

static void add_string(char *buf, int *j, int l1, int l2, char *str) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*s", l1, l2, str);
	else
	    *j += sprintf(buf + *j, "%*s", l1, str);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*s", l2, str);
	else
	    *j += sprintf(buf + *j, "%s", str);
}


/*
 * Formats tag information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 *
 * %%	Single % sign
 * %p	Tag position
 * %t	Tag type (always 4 characters)
 * %l	Tag length
 * %#	Tag number (0 if unknown)
 * %c	Tag comment
 *
 * Additionally, some formats (p, l, n and c) can be specified as
 * %<number><format> (eg %80c) to allow AT MOST that many characters.
 */
char *edGetBriefTag(edview *xx, tg_rec anno_ele, char *format) {
    static char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    GapIO *io = xx->io;
    anno_ele_t *e;

    if (!anno_ele)
	return "";

    e = (anno_ele_t *)cache_search(io, GT_AnnoEle, anno_ele);
    
    for (i = j = 0; format[i]; i++) {
	char type[5];

	if (format[i] != '%') {
	    status_buf[j++] = format[i];
	    continue;
	}

	l1 = strtol(&format[++i], &cp, 10);
	i = cp - format;
	if (format[i] == '.') {
	    l2 = strtol(&format[++i], &cp, 10);
	    i = cp - format;
	} else {
	    l2 = 0;
	}
	if (format[i] == 'R') {
	    raw = 1;
	    i++;
	} else {
	    raw = 0;
	}

	switch(format[i]) {
	case '%':
	    status_buf[j++] = '%';
	    break;

	case 't': /* Type */
	    (void)type2str(e->tag_type, type);
	    status_buf[j++] = type[0];
	    status_buf[j++] = type[1];
	    status_buf[j++] = type[2];
	    status_buf[j++] = type[3];
	    break;

	case 'p': { /* Position */
	    range_t *r = anno_get_range(io, anno_ele, NULL, 0);
	    add_number(status_buf, &j, l1, l2, r->start);
	    break;
	}

	case 'l': { /* Length */
	    range_t *r = anno_get_range(io, anno_ele, NULL, 0);
	    add_number(status_buf, &j, l1, l2, r->end - r->start + 1);
	    break;
	}

	case '#': /* Number */
	    add_number64(status_buf, &j, l1, l2, e->rec);
	    break;


	case 'c': /* Comment */
	    add_string(status_buf, &j, l1, l2,
		       e->comment ? e->comment : "(no comment)");
	    break;

	default:
	    status_buf[j++] = format[i];
	}
    }
    status_buf[j] = 0;

    return status_buf;
}


/*
 * Formats reading information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 *
 * %%	Single % sign
 * %n	Reading name (Raw: number)
 * %#	Reading number
 * %t	Trace name
 * %p	Position
 * %l	Clipped length
 * %L	Total length
 * %s	Start of clip
 * %e	End of clip
 * %m   Mapping quality
 * %S   Sense (whether complemented, +/-, Raw: 0/1)
 * %a	Chemistry (primer/terminator, Raw: integer)
 * %d	Strand (+/-, Raw 0/1)
 * %i	Reading freetext 'info' comment
 * %P	Primer (unknown/forward universal/reverse universal/forward custom/
 *              reverse custom,  Raw: 0/1/2/3/4)
 * %t   Trace name
 * %Tn	Template name (Raw: template number)
 * %T#	Template number
 * %Tv	Template vector (Raw: template vector number)
 * %Tc	Template consistency (Raw: as a number)
 * %Ti	Template insert size
 * %Cn	Clone name (Raw: clone number)
 * %C#	Clone number
 * %Cv	Clone vector (Raw: clone vector number)
 * %b   Base call
 * %c   Base confidence
 * %A   A confidence log-odds (raw for probability value)
 * %C   C confidence log-odds (raw for probability value)
 * %G   G confidence log-odds (raw for probability value)
 * %T   T confidence log-odds (raw for probability value)
 * %V   Vendor/platform
 *
 * Additionally specifying %<number><format> forces AT MOST that many
 * characters to be displayed.
 * Specifying %R<format> (or %<number>R<format>) requests the raw data to
 * be displayed. This only works for some fields. Eg %Rp displays 0 to 4, but
 * %p displays, for instance, "forward universal"
 * Specying %*<format> indicates that the above formats should be applied to
 * the other end of a read-pair instead of this read.
 * The special format %** is used to terminate decoding of the format if
 * the sequence is single-ended.
 */
char *edGetBriefSeq(edview *xx, tg_rec seq, int pos, int nth, char *format) {
    static char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    GapIO *io = xx->io;
    seq_t *s1 = get_seq(io, seq), *s2 = NULL, *s;
    tg_rec pair = 0;
    
    for (i = j = 0; format[i]; i++) {
	if (format[i] != '%') {
	    status_buf[j++] = format[i];
	    continue;
	}

	l1 = strtol(&format[++i], &cp, 10);
	i = cp - format;
	if (format[i] == '.') {
	    l2 = strtol(&format[++i], &cp, 10);
	    i = cp - format;
	} else {
	    l2 = 0;
	}
	if (format[i] == 'R') {
	    raw = 1;
	    i++;
	} else {
	    raw = 0;
	}

	if (format[i] == '*') {
	    if (pair == 0)
		pair = sequence_get_pair(io, s1);
	    if (pair > 0 && !s2)
		s2 = get_seq(io, pair);
	    s = s2 ? s2 : s1;
	    i++;
	} else {
	    s = s1;
	}

	switch(format[i]) {
	case '%':
	    status_buf[j++] = '%';
	    break;

	case '*':
	    if (!s2)
		goto bail_out;
	    break;

	case '#':
	    add_number64(status_buf, &j, l1, l2, s->rec);
	    break;

	case 'n':
	    if (raw)
		add_number64(status_buf, &j, l1, l2, s->rec);
	    else
		add_string(status_buf, &j, l1, l2, s->name);
	    break;

	case 'p': {
	    tg_rec cnum;
	    int cpos;
	    if (0 == sequence_get_position(xx->io, s->rec, &cnum, &cpos, NULL,
					   NULL)) {
		
		if (raw || cnum == xx->contig->rec) {
		    add_number(status_buf, &j, l1, l2, cpos);
		} else {
		    char buf[1024];
		    sprintf(buf, "%d@%s", cpos, get_contig_name(io, cnum));
		    add_string(status_buf, &j, l1, l2, buf);
		}
	    }
	    break;
	}

	case 'l':
	    add_number(status_buf, &j, l1, l2, ABS(s->len));
	    break;

	case 'L':
	    add_number(status_buf, &j, l1, l2, s->right - s->left + 1);
	    break;

	case 's':
	    add_number(status_buf, &j, l1, l2, s->left);
	    break;

	case 'e':
	    add_number(status_buf, &j, l1, l2, s->right);
	    break;

	case 'S':
	    if (raw)
		add_number(status_buf, &j, l1, l2, s->len < 0);
	    else
		add_string(status_buf, &j, l1, l2, s->len < 0 ? "<<" : ">>");
	    break;

	case 'd':
	    {
		int strand = sequence_get_len(&s) < 0 ? 1 : 0;

		if (raw)
		    add_number(status_buf, &j, l1, l2, strand);
		else {
		    char *str;
		    if      (strand == 0) str = "+";
		    else if (strand == 1) str = "-";
		    else                  str = "?";
		    add_string(status_buf, &j, l1, l2, str);
		}
	    }
	    break;

	case 'b': {
	    char base[2];
	    int cut;
	    if (0 == sequence_get_ubase(xx->io, &s, pos, nth,
					&base[0], NULL, &cut)) {
		base[1] = 0;
		if (cut)
		    base[0] = tolower(base[0]);
		else
		    base[0] = toupper(base[0]);
		add_string(status_buf, &j, l1, l2, base);
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;
	}
	case 'c': {
	    int q;
	    if (0 == sequence_get_ubase(xx->io, &s, pos, nth,
					NULL, &q, NULL)) {
		if (raw) {
		    add_double(status_buf, &j, l1, l2, 1 - pow(10, q/-10.0));
		} else {
		    add_number(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;
	}
	case 'A':
	    if (pos >= 0 && pos < ABS(s->len)) {
		double q[4];
		sequence_get_base4(xx->io, &s, pos, NULL, q, NULL, 1);
		if (raw)
		    add_double(status_buf, &j, l1, l2, exp(q[0]));
		else {
		    double p = exp(q[0]);
		    add_double(status_buf, &j, l1, l2, 4.342945*log(p/(1-p)));
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'C':
	    if (pos >= 0 && pos < ABS(s->len)) {
		double q[4];
		sequence_get_base4(xx->io, &s, pos, NULL, q, NULL, 1);
		if (raw)
		    add_double(status_buf, &j, l1, l2, exp(q[1]));
		else {
		    double p = exp(q[1]);
		    add_double(status_buf, &j, l1, l2, 4.342945*log(p/(1-p)));
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'G':
	    if (pos >= 0 && pos < ABS(s->len)) {
		double q[4];
		sequence_get_base4(xx->io, &s, pos, NULL, q, NULL, 1);
		if (raw)
		    add_double(status_buf, &j, l1, l2, exp(q[2]));
		else {
		    double p = exp(q[2]);
		    add_double(status_buf, &j, l1, l2, 4.342945*log(p/(1-p)));
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'T':
	    if (pos >= 0 && pos < ABS(s->len)) {
		double q[4];
		sequence_get_base4(xx->io, &s, pos, NULL, q, NULL, 1);
		if (raw)
		    add_double(status_buf, &j, l1, l2, exp(q[3]));
		else {
		    double p = exp(q[3]);
		    add_double(status_buf, &j, l1, l2, 4.342945*log(p/(1-p)));
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'm':
	    if (raw) {
		add_double(status_buf, &j, l1, l2,
			   1 - pow(10, s->mapping_qual/-10.0));
	    } else {
		add_number(status_buf, &j, l1, l2, s->mapping_qual);
	    }
	    break;

	case 'V':
	    if (raw) {
		add_number(status_buf, &j, l1, l2, s->seq_tech);
	    } else {
		switch(s->seq_tech) {
		case STECH_SANGER:
		    add_string(status_buf, &j, l1, l2, "Sanger");
		    break;
		case STECH_SOLEXA:
		    add_string(status_buf, &j, l1, l2, "Illumina");
		    break;
		case STECH_SOLID:
		    add_string(status_buf, &j, l1, l2, "SOLiD");
		    break;
		case STECH_454:
		    add_string(status_buf, &j, l1, l2, "454");
		    break;
		default:
		    add_string(status_buf, &j, l1, l2, "unknown");
		    break;
		}
	    }
	    break;
	    
	default:
	    status_buf[j++] = format[i];
	}
    }
 bail_out:
    status_buf[j] = 0;

    return status_buf;
}

/*
 * Formats consensus information for the status line.
 * This is done using a format string where certain % rules are replaced by
 * appropriate components.
 *
 * %%	Single % sign
 * %n	Contig name
 * %#	Contig number
 * %p	Position
 * %l	Length
 * %s	Start of clip
 * %e	End of clip
 * %b   Base call
 * %c   Base confidence log-odds (raw for probability value)
 * %A   A confidence log-odds (raw for probability value)
 * %C   C confidence log-odds (raw for probability value)
 * %G   G confidence log-odds (raw for probability value)
 * %T   T confidence log-odds (raw for probability value)
 * %*   * (gap) confidence
 *
 * Additionally specifying %<number><format> forces AT MOST that many
 * characters to be displayed.
 * Specifying %R<format> (or %<number>R<format>) requests the raw data to
 * be displayed. This only works for some fields. Eg %Rp displays 0 to 4, but
 * %p displays, for instance, "forward universal"
 */
char *edGetBriefCon(edview *xx, tg_rec crec, int pos, int nth, char *format) {
    static char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    
    for (i = j = 0; format[i]; i++) {
	if (format[i] != '%') {
	    status_buf[j++] = format[i];
	    continue;
	}

	l1 = strtol(&format[++i], &cp, 10);
	i = cp - format;
	if (format[i] == '.') {
	    l2 = strtol(&format[++i], &cp, 10);
	    i = cp - format;
	} else {
	    l2 = 0;
	}
	if (format[i] == 'R') {
	    raw = 1;
	    i++;
	} else {
	    raw = 0;
	}

	switch(format[i]) {
	case '%':
	    status_buf[j++] = '%';
	    break;

	case '#':
	    add_number64(status_buf, &j, l1, l2, crec);
	    break;

	case 'n':
	    if (raw)
		add_number64(status_buf, &j, l1, l2, crec);
	    else
		add_string(status_buf, &j, l1, l2,
			   contig_get_name(&xx->contig));
	    break;

	case 'p': {
	    add_number(status_buf, &j, l1, l2, pos);
	    break;
	}

	case 'l':
	    add_number(status_buf, &j, l1, l2, contig_get_length(&xx->contig));
	    break;

	case 's':
	    add_number(status_buf, &j, l1, l2, contig_get_start(&xx->contig));
	    break;

	case 'e':
	    add_number(status_buf, &j, l1, l2, contig_get_end(&xx->contig));
	    break;

	case 'b':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		char base[2];
		base[0] =  xx->displayedConsensus[pos - xx->displayPos];
		base[1] = 0;
		add_string(status_buf, &j, l1, l2, base);
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'c':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q =
		    xx->cachedConsensus[p].scores[xx->cachedConsensus[p].call];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'A':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[0];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'C':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[1];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'G':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[2];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'T':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[3];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case '*':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[4];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	default:
	    status_buf[j++] = format[i];
	}
    }
    status_buf[j] = 0;

    return status_buf;
}

static int
tk_redisplaySeqSequences_callback(void *cd, pileup_t *p, int pos, int nth) {
    int i;
    edscreen_t *c = (edscreen_t *)cd;
    edview *xx = c->xx;

    /* Abort once we're past the max column */
    if (c->ncols >= xx->displayWidth+9)
	return 1;

    /* clear column */
    for (i = 0; i < c->nrows; i++) {
	c->str[i][c->ncols] = ' ';
	c->ink[i][c->ncols].sh = sh_default;
    }

    /* Empty columns just get the previous coord */
    if (!p && c->ncols) {
	c->pos[c->ncols] = c->pos[c->ncols-1];
	c->nth[c->ncols] = c->nth[c->ncols-1];
    }

    /* store data from pileup */
    //printf("Pos %d/%d %p ", pos, nth, p);
    for (; p; p = p->next) {
	int row = p->r->y;
	char *s;
	XawSheetInk *i;

	if (p->start) {
	    int sp, sn, rp, rn;
	    sequence_cigar2pos(xx->io, p->s, p->cigar_ind, p->cigar_len,
			       &sp, &sn, &rp, &rn);
	    /*
	    printf("spos=%d/%d rpos=%d/%d\n",
		   sp, sn, rp, rn);
	    printf("%d/%d: New seq at row %d, ind %d, rec %"PRIrec
		   " %d/%d %d%c %.10s\n",
		   pos, nth, row, (int)(p->r - xx->r), p->r->rec,
		   c->ncols, p->seq_offset, p->cigar_len, p->cigar_op,
		   &p->s->seq[p->seq_offset]);
	    */
	    p->start = 0;
	    c->rpos[(int)(p->r - xx->r)] = c->ncols ? c->ncols: -p->seq_offset;
	    c->cigar_ind[(int)(p->r - xx->r)] = p->cigar_ind;
	    c->cigar_len[(int)(p->r - xx->r)] = p->cigar_len;

	    if (c->rec[row] == -1) {
		/*
		printf("Set c->rec[%d] = %d\n", 
		       row, (int)(p->r - xx->r));
		*/
		c->rec[row] = (int)(p->r - xx->r);
	    }
	}

	//putchar(p->base);
	//printf("%d ", row);

	s = &c->str[row][c->ncols];
	i = &c->ink[row][c->ncols];
	
	*s = p->base;
	c->pos[c->ncols] = pos;
	c->nth[c->ncols] = nth;

	/* Cutoffs */
	if (p->sclip) {
	    if (xx->ed->display_cutoffs) {
		i->sh |= sh_light;
	    } else {
		*s = ' ';
	    }
	}

	/* Quality values */
	if (xx->ed->display_quality &&
	    (p->sclip == 0 || xx->ed->display_cutoffs)) {
	    int qbin = p->qual / 10;
	    if (qbin < 0) qbin = 0;
	    if (qbin > 9) qbin = 9;
	    i->sh |= sh_bg;
	    i->bg = xx->ed->qual_bg[qbin]->pixel;
	}

	/* Highlight disagreements */
	if (xx->ed->display_differences &&
	    (p->sclip == 0 || xx->ed->display_cutoffs)) {
	    char ubase = xx->ed->display_differences_case
		? p->base : toupper(p->base);

	    switch (xx->ed->display_differences) {
	    case 1:
		if (ubase == xx->displayedConsensus[c->ncols])
		    *s = '.';
		else if (p->qual < xx->ed->display_differences_qual)
		    *s = ':';
		break;

	    case 2:
		if (ubase != xx->displayedConsensus[c->ncols]) {
		    i->sh |= sh_fg;
		    if (p->qual >= xx->ed->display_differences_qual)
			i->fg = xx->ed->diff2_fg->pixel;
		    else
			i->fg = xx->ed->diff1_fg->pixel;
		}
		break;

	    case 3:
		if (ubase != xx->displayedConsensus[c->ncols]) {
		    i->sh |= sh_bg;
		    if (p->qual >= xx->ed->display_differences_qual)
			i->bg = xx->ed->diff2_bg->pixel;
		    else
			i->bg = xx->ed->diff1_bg->pixel;
		}
		break;
	    }
	}
    }
    //putchar('\n');

    c->ncols++;
    return 0;
}

static edscreen_t *
tk_redisplaySeqSequences_columns(edview *xx, rangec_t *r, int nr) {
    int ret;
    edscreen_t *c;
    int i, j, k, max_y = -1;

    /* Compute maximum Y dimension */
    for (i = 0; i < nr; i++) {
	if (max_y < xx->r[i].y)
	    max_y = xx->r[i].y;
    }
    if (max_y < xx->displayHeight)
	max_y = xx->displayHeight;
    max_y++;

    /* Allocate our screen cache */
    c = malloc(sizeof(*c));
    c->xx    = xx;
    c->ncols = 0;
    c->nrows = max_y;
    c->pos   = malloc((xx->displayWidth+9) * sizeof(int));
    c->nth   = calloc(xx->displayWidth+9, sizeof(int));
    c->str   = malloc(max_y * (sizeof(char *) + xx->displayWidth + 9));
    c->ink   = malloc((max_y * (sizeof(char *) + xx->displayWidth + 9))
		      * sizeof(XawSheetInk));
    c->rec   = malloc(max_y * sizeof(int));
    c->rpos  = malloc(nr * sizeof(int));
    c->cigar_ind = calloc(nr, sizeof(int));
    c->cigar_len = calloc(nr, sizeof(int));
    for (i = 0; i < xx->displayWidth+9; i++)
	c->pos[i] = INT_MIN;

    for (i = 0; i < max_y; i++) {
	c->str[i] = (char *)(&c->str[max_y])        + i * (xx->displayWidth+9);
	c->ink[i] = (XawSheetInk *)(&c->ink[max_y]) + i * (xx->displayWidth+9);
	c->rec[i] = -1;
    }
    for (i = 0; i < nr; i++)
	c->rpos[i] = INT_MAX;

    ret = pileup_loop(xx->io, r, nr,
		      xx->displayPos, 0,
		      xx->displayPos + xx->displayWidth + 9, 0,
		      tk_redisplaySeqSequences_callback, c);

    if (ret) {
	free(c);
	return NULL;
    }

    return c;
}

/*
 * Populates the cache of visible items in xx->r and xx->nr.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int edview_visible_items(edview *xx, int start, int end) {
    int i;
    int mode = xx->ed->stack_mode
	? CSIR_ALLOCATE_Y_MULTIPLE
	: CSIR_ALLOCATE_Y_SINGLE;

    /* sort... */
    mode |= CSIR_SORT_BY_SEQ_TECH;

    end += 9; /* So numbers always display correctly */

    /* Always reload for now as we can't spot edits yet */
    if ((xx->refresh_flags & ED_DISP_COLOUR) == 0 &&
	xx->r && xx->r_start == start && xx->r_end == end)
	return 0;

    /* Query sequences */
    if (xx->r)
	free(xx->r);
    xx->r_start = start;
    xx->r_end = end;
    xx->r = contig_items_in_range(xx->io, &xx->contig, start, end,
				  CSIR_SORT_BY_Y | mode, &xx->nr);
    if (!xx->r)
	return -1;

    if (xx->rec_hash) {
	HacheTableDestroy(xx->rec_hash, 0);
    }

    xx->rec_hash = HacheTableCreate(8192, HASH_DYNAMIC_SIZE);
    xx->rec_hash->name = "rec_hash";

    /* Work out Y dimension */
    xx->max_height = 0;
    for (i = 0; i < xx->nr; i++) {
	HacheData hd;
	tg_rec key = xx->r[i].rec;

	if (xx->max_height < xx->r[i].y)
	    xx->max_height = xx->r[i].y;

	hd.i = i;
	HacheTableAdd(xx->rec_hash, (char *)&key, sizeof(key), hd, NULL);
    }
    xx->max_height += 3; /* +1 for from 0, +2 for consensus+ruler */

    /* Fast map of annotations to sequences */
    if (xx->anno_hash)
	HacheTableDestroy(xx->anno_hash, 0);

    xx->anno_hash = HacheTableCreate(8192, HASH_DYNAMIC_SIZE |
				     HASH_ALLOW_DUP_KEYS);
    xx->anno_hash->name = "anno_hash";
    for (i = 0; i < xx->nr; i++) {
	tg_rec key = xx->r[i].pair_rec; /* aka obj_rec */
	HacheData hd;

	if ((xx->r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)
	    continue;

	/*
	 * Work around a bug (or design flaw?) in break_contig.
	 * When breaking a contig we need to reparent our tags to a new
	 * contig record. This isn't feasible as tags shouldn't "know"
	 * which contigs they're in. Instead we use the flag.
	 */
	if (!(xx->r[i].flags & GRANGE_FLAG_TAG_SEQ))
	    key = xx->cnum;

	hd.i = i;
	HacheTableAdd(xx->anno_hash, (char *)&key, sizeof(key), hd, NULL);
    }

    /*
     * Force consensus to be up to date as tk_redisplaySeqSequences_columns
     * makes use of this if we have highlight disagreements enabled
     */
    if (xx->ed->display_differences) {
	calculate_consensus(xx->io, xx->cnum, xx->displayPos,
			    xx->displayPos + xx->displayWidth - 1,
			    xx->cachedConsensus);
	for (i = 0; i < xx->displayWidth; i++) {
	    xx->displayedConsensus[i] = "ACGT*N"[xx->cachedConsensus[i].call];
	}
    }

    /* Update xx->screen cache. This holds data about sheet cells */
    if (xx->screen) {
	free(xx->screen->str);
	free(xx->screen->ink);
	free(xx->screen->pos);
	free(xx->screen->nth);
	free(xx->screen->rec);
	free(xx->screen->rpos);
	free(xx->screen->cigar_ind);
	free(xx->screen->cigar_len);
	free(xx->screen);
    }
    xx->screen = tk_redisplaySeqSequences_columns(xx, xx->r, xx->nr);

#if 0
    puts("");
    for (i = 0; i < xx->nr; i++) {
	tg_rec rec;

	if ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    anno_ele_t *a = (anno_ele_t *)cache_search(xx->io, GT_AnnoEle,
						       xx->r[i].rec);
	    rec = a->rec;
	} else {
	    seq_t *s = (seq_t *)cache_search(xx->io, GT_Seq,
					     xx->r[i].rec);
	    rec = s->rec;
	}
	printf("%d\t%d\t%s%d/%d\t%d\t%s\t%d..%d\n",
	       i, xx->r[i].y,
	       xx->r[i].rec == rec ? "" : "*", xx->r[i].rec, rec,
	       ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
	           ? xx->r[i].pair_rec : 0,
	       ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
	           ? "tag" : "seq",
	       xx->r[i].start, xx->r[i].end);
    }
#endif


    return 0;
}

/*
 * Given an X,Y coordinate return the reading id under this position.
 *
 * Returns record number on success
 *         -1 on failure.
 */
int edGetGelNumber(edview *xx, int x, int y) {
    if (xx->editorState == StateDown)
	return -1;

    if (y < 0 || y >= xx->displayHeight ||
	x < 0 || x >= xx->displayWidth)
	return -1;

    puts("edGetGelNumber unimplemented");
    return 0;
}

/*
 * Find the largest i such that r[i].y <= y.
 */
static int edview_binary_search_y(rangec_t *r, int nr, int y) {
    int i_start, i_end, i_mid;

    if (nr <= 0)
	return 0;

    i_start = 0;
    i_end = nr;

    while (i_start < i_end) {
	i_mid = (i_end - i_start) / 2 + i_start;

	if (r[i_mid].y < y) {
	    i_start = i_mid+1;
	} else {
	    i_end = i_mid;
	}
    }

    return i_mid;
}

static int ed_set_xslider_pos(edview *xx, int offset) {
    char buf[100];
    double len = contig_get_length(&xx->contig);

    offset -= contig_get_start(&xx->contig);

    sprintf(buf, " %.20f %.20f",
	    offset / len,
	    (offset + xx->displayWidth) / len);

    if (Tcl_VarEval(xx->interp, xx->ed->xScrollCmd, buf, NULL)
	!= TCL_OK) {
        Tcl_AddErrorInfo(xx->interp, "\n(xscrollcommand executed by Editor)");
        Tcl_BackgroundError(xx->interp);
	return -1;
    }

    return 0;
}

static int ed_set_yslider_pos(edview *xx, int offset, int size, int total) {
    char buf[100];
    double len = total;

    if (xx->ed->consensus_at_top) {
	sprintf(buf, " %.20f %.20f",
		offset / len,
		(offset + size) / len);
    } else {
	sprintf(buf, " %.20f %.20f",
		(len - offset - size) / len,
		(len - offset) / len);
    }

    if (Tcl_VarEval(xx->interp, xx->ed->yScrollCmd, buf, NULL)
	!= TCL_OK) {
        Tcl_AddErrorInfo(xx->interp, "\n(yscrollcommand executed by Editor)");
        Tcl_BackgroundError(xx->interp);
	return -1;
    }

    return 0;
}

/* Update X scrollbar of names display */
void ed_set_nslider_pos(edview *xx, int pos) {
    edNames *en = xx->names;
    char buf[1024];

    if (!en || xx->editorState == StateDown)
	return;

    if (en->xScrollCmd) {
	double fract1, fract2;

	if (xx->ed->stack_mode) {
	    fract1 = 0;
	    fract2 = 1;
	} else {
	    fract1 = pos / (double)MAX_NAME_LEN;
	    fract2 = (pos + en->sw.columns) / (double)MAX_NAME_LEN;
	}
	sprintf(buf, " %.20f %.20f", fract1, fract2);
	if (Tcl_VarEval(EDINTERP(en), en->xScrollCmd, buf, NULL) != TCL_OK) {
	    printf("Error in editor names scroll: %s\n", Tcl_GetStringResult(EDINTERP(en)));
	}
    }
}

static void tk_redisplaySeqTags(edview *xx, XawSheetInk *ink, seq_t *s,
				int sp, int left, int right) {
    char type[5];
    HacheItem *hi;
    tg_rec key = s ? s->rec : xx->cnum;
    int j;

    if (xx->ed->hide_annos)
	return;

    if (xx->nr == 0)
	return;

    for (hi = HacheTableSearch(xx->anno_hash,
			       (char *)&key, sizeof(key));
	 hi; hi = HacheTableNext(hi, (char *)&key, sizeof(key))) {
	int ai = hi->data.i;
	int db = idToIndex(type2str(xx->r[ai].mqual, type));
	
	printf("drawing tag %"PRIrec" %d+%d / %d+%d\n", xx->r[ai].rec,
	       xx->r[ai].start, xx->r[ai].start_nth,
	       xx->r[ai].end, xx->r[ai].end_nth);


	for (j = 0; j < xx->displayWidth; j++) {
	    if (xx->screen->pos[j] > xx->r[ai].start ||
		(xx->screen->pos[j] == xx->r[ai].start &&
		 xx->screen->nth[j] >= xx->r[ai].start_nth))
		break;
	}

	while (j < xx->displayWidth &&
	       xx->screen->pos[j] <= xx->r[ai].end &&
	       xx->screen->pos[j] != INT_MIN) {
	    if (tag_db[db].fg_colour!=NULL) {
		ink[j].sh|=sh_fg;
		ink[j].fg=tag_db[db].fg_pixel;
	    }
	    if (tag_db[db].bg_colour!=NULL) {
		ink[j].sh|=sh_bg;
		ink[j].bg=tag_db[db].bg_pixel;
	    }
		
	    if (xx->screen->pos[j] == xx->r[ai].end &&
		xx->screen->nth[j] >= xx->r[ai].end_nth)
		break;

	    j++;
	}
    }
}

/* Returns 1 if rec is in the global "readings" list, 0 if not */
static int seq_in_readings_list(edview *xx, tg_rec rec) {
    char srec[20];
    
    sprintf(srec, "#%"PRIrec, rec);
    return Tcl_GetVar2(xx->interp,
		       "NGList_read_hash",
		       srec, TCL_GLOBAL_ONLY) ? 1 : 0;
}

static void tk_redisplaySeqSequences(edview *xx, rangec_t *r, int nr) {
    int i, j, k, box_alt, seq_len;

    puts("");

    /*
    sheet_clear(&xx->ed->sw);
    sheet_clear(&xx->names->sw);
    */

    i = edview_binary_search_y(r, nr, xx->displayYPos);

    /* Work down the screen line by line */
    for (j = xx->y_seq_start;
	 j < xx->displayHeight - xx->y_seq_end && i < nr;
	 j++) {
	int sp, l;
	unsigned char seq_a[MAX_SEQ_LEN+1], *seq = seq_a;
	XawSheetInk ink[MAX_DISPLAY_WIDTH], nink[MAX_NAME_WIDTH];
	char line[MAX_DISPLAY_WIDTH+1], nline[MAX_NAME_WIDTH];
	int dir;
	int left, right;
	int seq_p;

	if (xx->refresh_flags & (ED_DISP_READS | ED_DISP_SEQ)) {
	    memset(ink, 0, MAX_DISPLAY_WIDTH * sizeof(*ink));
	    memset(line, ' ', MAX_DISPLAY_WIDTH);
	}
	if (xx->ed->stripe_mode) {
	    int n = xx->ed->stripe_mode;
	    for (k = n-xx->displayPos%n; k < xx->displayWidth; k+=n) {
		ink[k].sh |= sh_bg;
		ink[k].bg = xx->ed->stripe_bg->pixel;
	    }
	}

	if (xx->refresh_flags & (ED_DISP_NAMES | ED_DISP_NAME)) {
	    memset(nink, 0, MAX_NAME_WIDTH * sizeof(*nink));
	    memset(nline, ' ', MAX_NAME_WIDTH);
	}

	/* Iterate through all sequences on this line */
	box_alt = 0; /* alternating 0/1 */
	while (i < nr && xx->r[i].y - xx->displayYPos <= j - xx->y_seq_start) {
	    seq_t *s, *sorig;

	    if ((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO
#ifndef CACHED_CONS_VISIBLE
		|| (r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISCONS
#endif
		) {
		i++;
		continue;
	    }

	    if (xx->r[i].y - xx->displayYPos < j - xx->y_seq_start) {
		i++;
		continue;
	    }

	    s = sorig = get_seq(xx->io, r[i].rec);
	    sp = r[i].start;
	    l = s->len > 0 ? s->len : -s->len;
	    seq_p = 0;
	    dir = '>';

	    /* Optimisation for single sequence only */
	    if (xx->refresh_flags & ED_DISP_SEQ &&
		!(xx->refresh_flags & ED_DISP_SEQS)) {
		if (xx->refresh_seq != r[i].rec) {
		    continue;
		}
	    }
	
	    /* Complement data on-the-fly */
	    if ((s->len < 0) ^ r[i].comp) {
		dir = '<';
		s = dup_seq(s);
		complement_seq_t(s);
	    }

	    left = s->left;
	    right = s->right;

	    if (l > MAX_SEQ_LEN)
		seq = malloc(l);
	    memcpy(seq, s->seq, l);

	    if (sp < xx->displayPos) {
		seq_p += xx->displayPos - sp;
		//seq   += xx->displayPos - sp;
		//conf  += xx->displayPos - sp;
		l     -= xx->displayPos - sp;
		left  -= xx->displayPos - sp;
		right -= xx->displayPos - sp;
		sp = xx->displayPos;
	    }
	    if (l > xx->displayWidth - (sp-xx->displayPos))
		l = xx->displayWidth - (sp-xx->displayPos);

	    /* Sequence */
	    if (xx->refresh_flags & (ED_DISP_READS | ED_DISP_SEQ)) {
		int p, p2;

		/* Annotations */
		if (xx->anno_hash) {
		    tk_redisplaySeqTags(xx, xx->screen->ink[xx->r[i].y],
					s, sp, left, right);
		}

		memcpy(line, xx->screen->str[xx->r[i].y], xx->displayWidth);
		memcpy(ink,  xx->screen->ink[xx->r[i].y],
		       xx->displayWidth * sizeof(XawSheetInk));
	    }

	    /* Name */
	    if (xx->refresh_flags & (ED_DISP_NAMES | ED_DISP_NAME)) {
		int nl = s->name_len - xx->names_xPos;
		int ncol = xx->names->sw.columns;
		XColor **qual_bg = xx->ed->qual_bg;

		box_alt ^= 1;

		if (xx->ed->stack_mode) {
		    int p  = r[i].start - xx->displayPos;
		    int p2 = r[i].end   - xx->displayPos;
		    int bg = -1;
		    double nc = xx->names->sw.columns;
		    if (p < 0) p = 0;
		    p = p * (nc / xx->displayWidth);
		    if (p2 < 0) p2 = 0;
		    p2 = p2 * (nc / xx->displayWidth);
		    if (p2 > xx->names->sw.columns)
			p2 = xx->names->sw.columns;
		    while (nline[p] != ' ')
			p++;

		    if (seq_in_readings_list(xx, s->rec)) {
			int ptmp = p;
			do {
			    nink[ptmp++].sh |= box_alt ? sh_box : sh_box_alt;
			} while (ptmp < p2);
			qual_bg = xx->ed->qual_bg2;
			bg = xx->ed->qual_bg2[9]->pixel;
		    }

		    if (xx->ed->display_mapping_quality) {
			int qbin = s->mapping_qual / 10;
			if (qbin < 0) qbin = 0;
			if (qbin > 9) qbin = 9;
			bg = qual_bg[qbin]->pixel;
		    }

		    nline[p] = dir;
		    if (bg != -1) {
			nink[p].sh |= sh_bg;
			nink[p].bg = bg;
		    }
		    
		    for (++p; p < p2; p++) {
			nline[p] = '.';
			if (bg != -1) {
			    nink[p].sh |= sh_bg;
			    nink[p].bg = bg;
			}
		    }

		} else {
		    XColor **qual_bg = xx->ed->qual_bg;

		    nline[0] = dir;
		    if (nl > 0)
			memcpy(&nline[1], s->name + xx->names_xPos, nl);
		    nink[0].sh = sh_bg;
		    if (r[i].pair_rec) {
			nink[0].bg = xx->ed->qual_bg[9]->pixel;
		    } else {
			nink[0].bg = xx->ed->qual_bg[0]->pixel;
		    }
		    for (k = 1; k < ncol && k < 1024; k++) {
			nink[k].sh = sh_default;
		    }

		    if (seq_in_readings_list(xx, s->rec)) {
			qual_bg = xx->ed->qual_bg2;
			for (k = 1; k < ncol && k < MAX_NAME_WIDTH; k++) {
			    nink[k].sh |= sh_bg |
				(box_alt ? sh_box : sh_box_alt);
			    nink[k].bg = qual_bg[9]->pixel;
			}
		    }

		    if (xx->ed->display_mapping_quality) {
			int qbin = s->mapping_qual / 10;
			if (qbin < 0) qbin = 0;
			if (qbin > 9) qbin = 9;

			for (k = 1; k < ncol && k < MAX_NAME_WIDTH; k++) {
			    nink[k].sh |= sh_bg;
			    nink[k].bg = qual_bg[qbin]->pixel;
			}
		    }
		}
	    }

	    //cache_decr(xx->io, sorig);

	    if (s != sorig)
		free(s);

	    i++;

	    if (seq != seq_a) {
		free(seq);
		seq = seq_a;
	    }
	}

	if (xx->refresh_flags & (ED_DISP_READS | ED_DISP_SEQ))
	    XawSheetPutJazzyText(&xx->ed->sw, 0, j, xx->displayWidth,
				 line, ink);

	if (xx->refresh_flags & (ED_DISP_NAMES | ED_DISP_NAME))
	    XawSheetPutJazzyText(&xx->names->sw, 0, j, xx->names->sw.columns,
				 nline, nink);
    }

    /*
     * Clear any blank lines too.
     */
    if (xx->refresh_flags & ED_DISP_SEQS) {
	char line[MAX_DISPLAY_WIDTH];
	memset(line, ' ', MAX_DISPLAY_WIDTH);
	for (; j < xx->displayHeight - xx->y_seq_end; j++) {
	    XawSheetPutText(&xx->ed->sw, 0, j, xx->displayWidth, line);
	    XawSheetPutText(&xx->names->sw, 0, j, xx->names->sw.columns, line);
	}
    }
 
    //ed_set_xslider_pos(xx, xx->displayPos);
    //ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight, nr);

    return;
}


static void tk_redisplaySeqConsensus(edview *xx) {
    int pos = xx->displayPos;
    int wid =  xx->displayWidth;
    char name[] = " Consensus";
    int i;
    XawSheetInk ink[MAX_DISPLAY_WIDTH];

    /* Names panel */
    XawSheetPutText(&xx->names->sw, 0, xx->y_cons, strlen(name), name);

    /* Maybe already computed in edview_visible_items */
    if (!xx->ed->display_differences) {
	calculate_consensus(xx->io, xx->cnum, pos, pos+wid-1,
			    xx->cachedConsensus);
	for (i = 0; i < wid; i++) {
	    xx->displayedConsensus[i] = "ACGT*N"[xx->cachedConsensus[i].call];
	}
    }

    /* Editor panel */
    memset(ink, 0, MAX_DISPLAY_WIDTH * sizeof(*ink));
    if (xx->ed->display_quality) {
	int qbin;

	for (i = 0; i < wid; i++) {
	    qbin = xx->cachedConsensus[i].phred/10;
	    if (qbin < 0) qbin = 0;
	    if (qbin > 9) qbin = 9;
	    ink[i].sh |= sh_bg;
	    ink[i].bg = xx->ed->qual_bg[qbin]->pixel;
	}
    }

    /* Consensus annotations */
    if (xx->anno_hash)
	tk_redisplaySeqTags(xx, ink, NULL, 0, 0, 0);

    XawSheetPutJazzyText(&xx->ed->sw, 0, xx->y_cons, wid,
			 xx->displayedConsensus, ink);
}

/*
 * Returns the other editor in a join-editor pair.
 *         NULL if not joined.
 */
edview *linked_editor(edview *xx) {
    if (!xx->link)
	return NULL;
    return xx->link->xx[0] == xx ? xx->link->xx[1] : xx->link->xx[0];
}

static int tk_redisplaySeqDiffs(edview *xx) {
    char diff[MAX_DISPLAY_WIDTH+1];
    edview *xx0, *xx1;
    int i;

    xx0 = xx->link->xx[0];
    xx1 = xx->link->xx[1];

    for (i = 0; i < xx->displayWidth; i++) {
	diff[i] = xx0->displayedConsensus[i] == xx1->displayedConsensus[i]
	    ? ' ' : '!';
    }
    XawSheetPutText(&xx->link->diffs->sw, 0, 0,
		    xx->displayWidth, diff);
    return 0;
}

/*
 * Calculates the numbers for the contig editor ruler line.
 * The return value is the index into this ruler buffer to plot.
 */
static int generate_ruler(edview *xx, char *ruler, int pos, int width) {
    char *k = ruler;
    int j;

    //int padded_pos[MAX_DISPLAY_WIDTH+21];

    memset(ruler, ' ', MAX_DISPLAY_WIDTH+21);

#if 0
    if (DBI(xx)->reference_seq) {
	/* Number relative to a specific sequence number */
	char *bases = DBgetSeq(DBI(xx), DBI(xx)->reference_seq);
	int reflen = DB_Length(xx, DBI(xx)->reference_seq);
	int rp = DB_RelPos(xx, DBI(xx)->reference_seq);

	/* 
	 * Compute unpadded base positions for this window.
	 * WARNING: This code is full of warts. I carefully worked it out, but
	 * it was not quite perfect. Working out again gave a slightly
	 * different, but equally imperfect answer. So I applied the next
	 * logical step - I hope you like fudge!
	 *
	 * If the change the code - check for these cases and combinations:
	 *   Position from left end when uncomplemented
	 *   Position from right end when complemented
	 *   Negative offsets
	 *   When relpos of ref is 1 and when it is not, both when refseq
	 *     is complemented and when it is not.
	 *   Effect of pads in ref seq
	 *   That base numbers are never positioned above the pads
	 *   Circular sequences (both orientations)
	 *   Offset base numbers
	 */
	if ((DB_Comp(xx, DBI(xx)->reference_seq) == UNCOMPLEMENTED)) {
	    int unpadded = MIN(0, pos);

	    for (j = unpadded; j < pos + width + 9; j++) {
		if (j >= pos)
		    padded_pos[j-pos] = unpadded-rp;
		if (j-rp < 0 || j-rp >= reflen || bases[j-rp] != '*')
		    unpadded++;
	    }
	} else {
	    int unpadded = reflen - MAX(pos + width + 8, reflen);

	    for (j = MAX(pos + width + 8, reflen); j >= pos; j--) {
		if (j <= pos + width + 8)
		    padded_pos[j-pos] = unpadded+(rp-1);
		if (j-rp-1 < 0 || j-rp-1 >= reflen || bases[j-rp-1] != '*')
		    unpadded++;
	    }
	}

	for (j = pos; j < pos + width + 9; j++, k++) {
	    int unpadded2 = (padded_pos[j-pos] + DBI(xx)->reference_offset);

	    if (DBI(xx)->reference_len) {
		unpadded2 %= DBI(xx) -> reference_len;
		while (unpadded2 < 0) {
		    unpadded2 += DBI(xx)->reference_len;
		}
		if (unpadded2 == 0)
		    continue; /* Don't display 0 for circular seqs */
	    }
	    if (!(unpadded2 % 10)) {
		sprintf(k, "%10d", unpadded2);
	    }
	}

	/*
	 * Pads in the consensus can sometimes leave nulls in the buffer,
	 * which turn into square boxes on some fonts. Replace these with
	 * spaces.
	 */
	for (j = 0; j < MAX_DISPLAY_WIDTH+21; j++)
	    if (ruler[j] == '\0')
		ruler[j] = ' ';

	return 9;

    } else if (xx->unpadded_ruler) {
	/* Number by unpadded base calls, using the consensus */
	int unpadded;
	
	edUnpaddedBaseNumber(xx, pos, width + 9);
	for (j = pos; j < pos + width + 9; j++, k++) {
	    unpadded = edUnpaddedBaseNumber(xx, j, 0);
	    if (unpadded%10)
		continue;
	    sprintf(k, "%10d", unpadded);
	    if (k+10-ruler < MAX_DISPLAY_WIDTH+21)
		k[10]=' ';
	}
	edUnpaddedBaseNumber(xx, pos, -1);

	return 9;
    } /* else */
#endif

    /* Old padded numbering */
#if 0
    {
	/* Basic numbering */
	int lower,times;
	lower = (pos - pos%10);
	times = width/10 + 3;
	for (j=0;j<times;j++,k+=10,lower+=10)
	    sprintf(k,"%10d",lower);
	return 9+pos%10;
	
    }
#endif

    /* Unpadded data array, precomputed to xx->screen */
    if (xx->screen) {
	int last = INT_MAX;

	for (j = 0; j < xx->displayWidth+9; j++, k++) {
	    if (xx->screen->pos[j] == INT_MIN)
		continue;
	       
	    if (xx->screen->pos[j] != last && (xx->screen->pos[j]%10) == 0) {
		sprintf(k, "%10d", xx->screen->pos[j]);
		k[10] = ' ';
	    }
	    last = xx->screen->pos[j];
	}

	return 9;
    }

    return 0;
}

static void tk_redisplaySeqNumbers(edview *xx) {
    char ruler[MAX_DISPLAY_WIDTH+21];
    int off;

    off = generate_ruler(xx, ruler, xx->displayPos, xx->displayWidth);
    XawSheetPutText(&xx->ed->sw, 0, xx->y_numbers, xx->displayWidth,
			 &ruler[off]);
}


/* Handle scrolling changes */
static void tk_redisplaySeqScroll(edview *xx, rangec_t *r, int nr) {
    if (xx->refresh_flags & ED_DISP_YSCROLL) {
	ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight,
			   xx->max_height);

	/* No need to redraw the consenus/numbers */
	xx->refresh_flags |= ED_DISP_NAMES | ED_DISP_READS | ED_DISP_CURSOR |
	    ED_DISP_SELECTION;
    }

    if (xx->refresh_flags & ED_DISP_XSCROLL) {
	ed_set_xslider_pos(xx, xx->displayPos);

	/* Changing X may also change height */
	if (!(xx->refresh_flags & ED_DISP_YSCROLL))
	    ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight,
			       xx->max_height);

	xx->refresh_flags |= ED_DISP_ALL;
    }
}

/*
 * Forces the sheet widget to redisplay the cursor.
 */
static void tk_redisplayCursor(edview *xx, rangec_t *r, int nr) {
    int x, y;

    /* If visible, find the screen x/y coord */
    if (xx->cursor_rec == xx->cnum) {
	y = xx->y_cons;
    } else {
	tg_rec key;
	HacheItem *hi;
	
	key = xx->cursor_rec;
	if (!xx->rec_hash || !r)
	    return;

	hi = HacheTableSearch(xx->rec_hash, (char *)&key, sizeof(key));
	y = hi && hi->data.i < nr
	    ? r[hi->data.i].y + xx->y_seq_start - xx->displayYPos
	    : -1;

	if (y < xx->y_seq_start || y >= xx->displayHeight) {
	    XawSheetDisplayCursor(&xx->ed->sw, False);
	    return; /* not visible */
	}
    }

    //x = xx->cursor_apos - xx->displayPos;
    //XawSheetDisplayCursor(&xx->ed->sw, True);
    //XawSheetPositionCursor(&xx->ed->sw, x, y);

    if (xx->screen) {
	int i;
	for (i = 0; i < xx->displayWidth; i++) {
	    if (xx->screen->pos[i] == xx->cursor_apos) {
		XawSheetDisplayCursor(&xx->ed->sw, True);
		XawSheetPositionCursor(&xx->ed->sw, i+xx->cursor_anth, y);
		return;
	    }
	}
    }
    XawSheetDisplayCursor(&xx->ed->sw, False);
}


/*
 * Force the cursor to be visible. If x_safe or y_safe are true then
 * we omit some of the searching and assume there is no reason to check
 * that x or y is still visible.
 *
 * Returns 1 if redraw has taken place
 *         0 if not
 */
int showCursor(edview *xx, int x_safe, int y_safe) {
    int y_pos = 0;
    int do_x = 0;
    int do_y = 0;

    /* X position */
    if (!x_safe) {
	if (xx->screen) {
	    int w = xx->displayWidth > 10 ? 10 : xx->displayWidth;
	    if (xx->cursor_apos < xx->screen->pos[0]) {
		set_displayPos(xx, xx->cursor_apos -w);
		do_x = 1;
	    }
	    if (xx->cursor_apos > xx->screen->pos[xx->displayWidth-1] ||
		(xx->cursor_apos == xx->screen->pos[xx->displayWidth-1] &&
		 xx->cursor_anth > xx->screen->nth[xx->displayWidth-1])) {
		set_displayPos(xx, xx->cursor_apos - w);
		do_x = 1;
	    }
	} else {
	    do_x = 1;
	}
    }

    if (do_x)
	y_safe = 0;

    /* Y position */
    if (!y_safe && xx->cursor_type != GT_Contig) {
	int i;
	tg_rec key;
	int sheight = xx->displayHeight - xx->y_seq_end - xx->y_seq_start;
	HacheItem *hi;
	
	edview_visible_items(xx, xx->displayPos,
			     xx->displayPos + xx->displayWidth);

	key = xx->cursor_rec;
	if (!xx->rec_hash)
	    return 0;
	hi = HacheTableSearch(xx->rec_hash, (char *)&key, sizeof(key));
	if (!hi)
	    return 0;
	i = hi->data.i;

	y_pos = xx->r[i].y;
	if (y_pos == -1) {
	    y_pos = 0; /* tag on consensus */
	    xx->cursor_rec = xx->cnum;
	    xx->cursor_type = GT_Contig;
	}

	/* If above, scroll so this is the first row */
	if (y_pos < xx->displayYPos) {
	    xx->displayYPos = y_pos;
	    do_y = 1;
	}

	/* If below, scroll so this is the last row */
	if (y_pos >= xx->displayYPos + sheight) {
	    xx->displayYPos = y_pos - sheight + 1;
	    do_y = 1;
	}
    }

    tman_reposition_traces(xx, xx->cursor_apos, 0);

    if (do_x)
	ed_set_xslider_pos(xx, xx->displayPos);

    if (do_y)
	ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight,
			   xx->max_height);

    if (do_x || do_y) {
	xx->refresh_flags = ED_DISP_ALL;
	if (!do_x)
	    xx->refresh_flags &= ~(ED_DISP_CONS | ED_DISP_XSCROLL);
	edview_redraw(xx);
	return 1;
    }

    return 0;
}

/*
 * Sends out a notification of our cursor movement
 */
static void cursor_notify(edview *xx) {
    reg_cursor_notify cn;

    if (!xx->cursor)
	return;

    xx->cursor->seq = xx->cursor_rec;
    xx->cursor->pos = xx->cursor_pos;
    xx->cursor->nth = xx->cursor_nth;
    xx->cursor->abspos = xx->cursor_apos;
    //xx->cursor->absnth = xx->cursor_anth;
    xx->cursor->job = CURSOR_MOVE;
    xx->cursor->sent_by = xx->reg_id;
    cn.job = REG_CURSOR_NOTIFY;
    cn.cursor = xx->cursor;
    contig_notify(xx->io->base, xx->cnum, (reg_data *)&cn);
}

/*
 * Returns  1 if seq is visible in the current editor position.
 *          0 if not.
 *
 * It may also fill out the 'new_y' (if non-null) value to indicate
 * the preferred new displayYPos parameter. This is set to -1 when
 * there is no point in changing it as 'seq' will never be visible
 * (it's not overlapping this X coord).
 */
int edview_seq_visible(edview *xx, tg_rec seq, int *new_y) {
    int i, y_pos, vis = 0;
    int sheight = xx->displayHeight - xx->y_seq_end - xx->y_seq_start;
    HacheItem *hi;
	
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    if (new_y)
	*new_y = xx->displayYPos;

    if (!xx->rec_hash)
	return 0;
    hi = HacheTableSearch(xx->rec_hash, (char *)&seq, sizeof(seq));
    if (!hi || !xx->r)
	return 0;
    i = hi->data.i;

    y_pos = xx->r[i].y;
    vis = 1;
    if (y_pos == -1) {
	/* tag? */
	return 1;
    }

    if (!vis) {
	if (new_y)
	    *new_y = -1;
	return 0;
    }

    /* seq is above, scroll so this is the first row */
    if (y_pos < xx->displayYPos) {
	if (new_y)
	    *new_y = y_pos;
	return 0;
    }

    /* seq is below, scroll so this is the last row */
    if (y_pos >= xx->displayYPos + sheight) {
	if (new_y)
	    *new_y = y_pos - sheight + 1;;
	return 0;
    }

    /* Otherwise it's already on-screen */
    if (new_y)
	*new_y = y_pos;
    return 1;
}

/*
 * This is called by an editor xview method - ie the X scrollbar. It also
 * gets called programmatically on other conditions (such as making the
 * editing cursor visibile).
 *
 * We attempt to make sure that sequences that were previously visible on
 * screen are now also visible on screen (where possible). To achieve
 * this we may have to adjust the Y scrollbar position too.
 */
int set_displayPos(edview *xx, int pos) {
    char buf[100];
    int i, ret = 0;
    int delta = pos - xx->displayPos;
    edview *xx2[2];

    if (xx->link && xx->link->locked)
	xx = xx->link->xx[0];

    for (i = 0; i < 2; i++) {
	int new_y = -1, vis_pos;
	tg_rec vis_rec1, vis_rec2;
	int sheight;
	tg_rec vis_cur;

	xx2[i] = xx;

	if (!xx)
	    break;

	sheight = xx->displayHeight - xx->y_seq_end - xx->y_seq_start;

	edview_visible_items(xx, xx->displayPos,
			     xx->displayPos + xx->displayWidth);
	/*
	 * Pick an appropriate sequence to try and keep track of.
	 * This is one we'll attempt to keep visible before and after
	 * scrolling.
	 */
	vis_cur = edview_seq_visible(xx, xx->cursor_rec, NULL);

	edview_item_at_pos(xx, xx->y_seq_start,
			   0, 0, 0, 1, &vis_rec1, &vis_pos, NULL);
	edview_item_at_pos(xx, xx->displayHeight - xx->y_seq_end - 1,
			   0, 0, 0, 1, &vis_rec2, &vis_pos, NULL);
	    
	xx->displayPos += delta;

	sprintf(buf, "%d", pos);
	Tcl_SetVar2(xx->interp, xx->edname, "displayPos", buf,
		    TCL_GLOBAL_ONLY);

	xx->refresh_flags = ED_DISP_XSCROLL;
	if (i == 1)
	    xx->refresh_flags |= ED_DISP_NO_DIFFS;

	/* Try and ensure 'vis_rec' is still visible */
	if (vis_rec1 == -1 || !edview_seq_visible(xx, vis_rec1, &new_y)) {
	    if (new_y == -1 && vis_rec2 != -1) {
		if (edview_seq_visible(xx, vis_rec2, &new_y)) {
		    /* Already visible, but new_y is bottom loc */
		    new_y -= sheight-1;
		}
	    }

	    if (new_y != -1) {
		xx->displayYPos = new_y;
		xx->refresh_flags |= ED_DISP_YSCROLL;
	    }
	} else {
	    /* Still visible, but potentially changed Y */
	    if (new_y != -1 && new_y != xx->displayYPos) {
		xx->displayYPos = new_y;
		xx->refresh_flags |= ED_DISP_YSCROLL;
	    }
	}

	/* If editing cursor was visible, ensure it still is too */
	if (vis_cur) {
	    if (!edview_seq_visible(xx, xx->cursor_rec, &new_y)) {
		xx->displayYPos = new_y;
		xx->refresh_flags |= ED_DISP_YSCROLL;
	    }
	}

	if (xx->displayYPos + sheight > xx->nr) {
	    xx->displayYPos = xx->nr - sheight;
	    xx->refresh_flags |= ED_DISP_YSCROLL;
	}

	if (xx->displayYPos < 0) {
	    xx->displayYPos = 0;
	    xx->refresh_flags |= ED_DISP_YSCROLL;
	}

	xx = (xx->link && xx->link->locked)
	    ? xx->link->xx[1] : NULL;
    }

    if (xx2[0]->link)
	xx2[0]->link->lockOffset =
	    xx2[0]->link->xx[1]->displayPos - xx2[0]->link->xx[0]->displayPos;


    if (xx2[1]) ret |= edview_redraw(xx2[1]);
    if (xx2[0]) ret |= edview_redraw(xx2[0]);

    return ret;
}

/*
 * Resets the absolute position in the contig based on the contents of the
 * cursor_pos, cursor_nth and cursor_rec fields.
 *
 * cursor_pos/nth are the unpadded position and nth base at position in
 * either the consensus, sequence or anno. The only adjustment needed to make
 * this an absolute coordinate is to add the starting point for the data
 * item.
 */
void edSetApos(edview *xx) {
    switch (xx->cursor_type) {
    case GT_Contig:
	xx->cursor_apos = xx->cursor_pos;
	xx->cursor_anth = xx->cursor_nth;
	break;

    case GT_Seq: {
	tg_rec cnum;
	int cpos;
	sequence_get_position(xx->io, xx->cursor_rec, &cnum, &cpos, NULL,
			      NULL);
	xx->cursor_apos = cpos + xx->cursor_pos;
	xx->cursor_anth = xx->cursor_nth;
	break;
    }

    case GT_AnnoEle: {
	tg_rec cnum;
	range_t *r = anno_get_range(xx->io, xx->cursor_rec, &cnum, 0);
	xx->cursor_apos = r->start + xx->cursor_pos;
	xx->cursor_anth = xx->cursor_nth;
	break;
    }

    default:
	fprintf(stderr, "Unknown item type in edSetApos(): %d\n",
		xx->cursor_type);
    }

    /* Send a notification of cursor movement */
    cursor_notify(xx);
}

int edSetCursorPos(edview *xx, int type, tg_rec rec, int pos, int nth,
		   int visible) {
    if (!xx)
	return 0;

    if (type == GT_Seq) {
	seq_t *s = get_seq(xx->io, rec);
	int left = s->left-1;
	int right = s->right;
	int ulen = sequence_unpadded_len(s, NULL);

	if (xx->ed->display_cutoffs) {
	    left = 0;
	    right = sequence_unpadded_len(s, NULL);
	} else {
	    if (sequence_get_orient(xx->io, rec)) {
		s = get_seq(xx->io, rec);
		left  = ulen - (s->right-1) -1;
		right = ulen - (s->left-1);
	    }
	}

	/* If out of bounds, punt it to the consensus */
	if (pos < left || pos > right) {
	    if (visible) {
		if (pos < 0 || pos > sequence_unpadded_len(s, NULL))
		    return 0;
		xx->ed->display_cutoffs = 1;
	    } else {
		return 0;
	    }
	}
    }

    if (type != GT_Seq) {
	int ustart, uend;

	if (xx->ed->display_cutoffs) {
	    ustart = xx->contig->start;
	    uend   = xx->contig->end;
	} else {
	    char con;
	    calculate_consensus_simple(xx->io, xx->cnum, pos, pos, &con, NULL);
	    if (con != 'N') {
		/* Must be valid, so fake boundaries */
		ustart = pos;
		uend = pos;
	    } else {
		consensus_valid_range(xx->io, xx->contig->rec, &ustart, &uend);
	    }
	}

	uend++;
	if (pos < ustart)
	    pos = ustart;
	if (pos > uend)
	    pos = uend;
    }

    xx->cursor_type = type;
    xx->cursor_rec  = rec;
    xx->cursor_pos  = pos;
    xx->cursor_nth  = nth;

    edSetApos(xx);

    if (visible && !showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    } else {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

/*
 * Convert from a record and position to a window X,Y coordinate in sheet
 * units.
 *
 * Returns 0 on success and stores via x and y pointers.
 *        -1 on failure (rec/pos not visible).
 */
int edGetXY(edview *xx, int rec_type, tg_rec rec, int pos, int *x, int *y) {
    int i;

    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    if (xx->nr == 0)
	return -1;

    if (rec == xx->contig->rec) {
	int col = pos - xx->displayPos;

	if (col < 0 || col > xx->displayWidth)
	    return -1;
	
	*x = col;
	*y = 0;
	return 0;
    }

    for (i = 0; i < xx->nr; i++) {
	if (xx->r[i].rec == rec) {
	    int row, col;

	    row = xx->r[i].y + xx->y_seq_start - xx->displayYPos;
	    col = xx->r[i].start - xx->displayPos + pos;

	    if (col < 0 || col >= xx->displayWidth)
		return -1;

	    if (row < xx->y_seq_start ||
		row >= xx->displayHeight - xx->y_seq_end)
		return -1;

	    *x = col;
	    *y = row;
	    return 0;
	}
    }

    return -1;
}


int edCursorUp(edview *xx) {
    int j;
    int cpos = xx->cursor_apos;

    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    if (xx->nr == 0)
	return 0;

    /* Find the current sequence number */
    if (xx->cursor_type == GT_Contig) {
	j = xx->nr;
    } else {
	HacheItem *hi;
	tg_rec key = xx->cursor_rec;

	if (!xx->rec_hash)
	    return 0;
	hi = HacheTableSearch(xx->rec_hash, (char *)&key, sizeof(key));
	if (!hi)
	    return 0;
	j = hi->data.i;
    }

    /* Step up until we find something overlapping */
    for (j--; j >= 0; j--) {
	if (xx->r[j].start <= cpos && xx->r[j].end+1 >= cpos
#ifndef CACHED_CONS_VISIBLE
	    && ((xx->r[j].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISCONS)
#endif
	    && ((xx->r[j].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)) {
	    if (!xx->ed->display_cutoffs) {
		seq_t *s = get_seq(xx->io, xx->r[j].rec);
		int left = s->left;
		int right = s->right;
		if (sequence_get_orient(xx->io, xx->r[j].rec)) {
		    int ulen;
		    s = get_seq(xx->io, xx->r[j].rec);
		    ulen = sequence_unpadded_len(s, NULL);
		    left  = ulen - (s->right-1);
		    right = ulen - (s->left-1);
		}
		if (cpos - xx->r[j].start < left-1 ||
		    cpos - xx->r[j].start > right)
		    continue; /* Sequence present, but hidden */
	    }
	    xx->cursor_type = GT_Seq;
	    xx->cursor_pos = cpos - xx->r[j].start;
	    xx->cursor_rec = xx->r[j].rec;
	    break;
	}
    }

    /* Otherwise we've hit the consensus */
    if (j < 0) {
	xx->cursor_type = GT_Contig;
	xx->cursor_rec = xx->cnum;
	xx->cursor_pos = cpos;
    }

    cursor_notify(xx);
    if (!showCursor(xx, 1, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edCursorDown(edview *xx) {
    int j;
    int cpos = xx->cursor_apos;

    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    if (xx->nr == 0)
	return 0;

    /* Find the current sequence number */
    if (xx->cursor_type == GT_Contig) {
	cpos = xx->cursor_pos;
	j = -1;
    } else {
	HacheItem *hi;
	tg_rec key = xx->cursor_rec;

	if (!xx->rec_hash)
	    return 0;
	hi = HacheTableSearch(xx->rec_hash, (char *)&key, sizeof(key));
	if (!hi)
	    return 0;
	j = hi->data.i;
	cpos = xx->r[j].start + xx->cursor_pos;
    }

    /* Step up until we find something overlapping */
    for (j++; j < xx->nr; j++) {
	if (xx->r[j].start <= cpos && xx->r[j].end+1 >= cpos
#ifndef CACHED_CONS_VISIBLE
	    && ((xx->r[j].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISCONS)
#endif
	    && ((xx->r[j].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)) {
	    if (!xx->ed->display_cutoffs) {
		seq_t *s = get_seq(xx->io, xx->r[j].rec);
		int left = s->left;
		int right = s->right;
		if (sequence_get_orient(xx->io, xx->r[j].rec)) {
		    int ulen;
		    s = get_seq(xx->io, xx->r[j].rec);
		    ulen = sequence_unpadded_len(s, NULL);
		    left  = ulen - (s->right-1);
		    right = ulen - (s->left-1);
		}
		if (cpos - xx->r[j].start < left-1 ||
		    cpos - xx->r[j].start > right)
		    continue; /* Sequence present, but hidden */
	    }
	    xx->cursor_type = GT_Seq;
	    xx->cursor_pos = cpos - xx->r[j].start;
	    xx->cursor_rec = xx->r[j].rec;
	    break;
	}
    }

    /* Otherwise we've hit the consensus */
    if (j >= xx->nr) {
	xx->cursor_type = GT_Contig;
	xx->cursor_rec = xx->cnum;
	xx->cursor_pos = cpos;
    }

    cursor_notify(xx);
    if (!showCursor(xx, 1, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edCursorLeft(edview *xx) {
    int i, pd = -1, nd = -xx->cursor_anth;
    for (i = 1; i < xx->displayWidth; i++) {
	if (xx->screen->pos[i] == xx->cursor_apos &&
	    xx->screen->nth[i] == xx->cursor_anth) {
	    pd = xx->screen->pos[i-1] - xx->cursor_apos;
	    nd = xx->screen->nth[i-1] - xx->cursor_anth;
	    break;
	}
    }

    if (xx->cursor_type == GT_Seq) {
	if (xx->ed->display_cutoffs) {
	    if (xx->cursor_pos > 0 || xx->cursor_nth > 0) {
		//xx->cursor_pos--;
		//xx->cursor_apos--;
		xx->cursor_apos += pd;
		xx->cursor_anth += nd;
		xx->cursor_pos  += pd;
		xx->cursor_nth  += nd;
	    }
	} else {
	    seq_t *s = get_seq(xx->io, xx->cursor_rec);
	    int left = s->left;

	    if (sequence_get_orient(xx->io, xx->cursor_rec)) {
		s = get_seq(xx->io, xx->cursor_rec);
		left = sequence_unpadded_len(s, NULL) - (s->right-1);
	    }

	    if (xx->cursor_pos >= left || xx->cursor_nth > 0) {
		//xx->cursor_pos--;
		//xx->cursor_apos--;
		xx->cursor_apos += pd;
		xx->cursor_anth += nd;
		xx->cursor_pos  += pd;
		xx->cursor_nth  += nd;
	    }

	}
    } else {
	//xx->cursor_pos--;
	//xx->cursor_apos--;
	xx->cursor_apos += pd;
	xx->cursor_anth += nd;
	xx->cursor_pos  += pd;
	xx->cursor_nth  += nd;
    }

    cursor_notify(xx);
    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edCursorRight(edview *xx) {
    int i, pd = 1, nd = -xx->cursor_anth;
    for (i = 0; i < xx->displayWidth; i++) {
	if (xx->screen->pos[i] == xx->cursor_apos &&
	    xx->screen->nth[i] == xx->cursor_anth) {
	    pd = xx->screen->pos[i+1] - xx->cursor_apos;
	    nd = xx->screen->nth[i+1] - xx->cursor_anth;
	    break;
	}
    }

    if (xx->cursor_type == GT_Seq) {
	seq_t *s = get_seq(xx->io, xx->cursor_rec);

	if (xx->ed->display_cutoffs) {
	    if (xx->cursor_pos < sequence_unpadded_len(s, NULL)) {
		//xx->cursor_pos++;
		//xx->cursor_apos++;
		xx->cursor_apos += pd;
		xx->cursor_anth += nd;
		xx->cursor_pos  += pd;
		xx->cursor_nth  += nd;
	    }
	} else {
	    int right = s->right;

	    if (sequence_get_orient(xx->io, xx->cursor_rec)) {
		s = get_seq(xx->io, xx->cursor_rec);
		right = sequence_unpadded_len(s, NULL) - (s->left-1);
	    }

	    if (xx->cursor_pos < right) {
		//xx->cursor_pos++;
		//xx->cursor_apos++;
		xx->cursor_apos += pd;
		xx->cursor_anth += nd;
		xx->cursor_pos  += pd;
		xx->cursor_nth  += nd;
	    }
	}
    } else {
	//xx->cursor_pos++;
	//xx->cursor_apos++;
	xx->cursor_apos += pd;
	xx->cursor_anth += nd;
	xx->cursor_pos  += pd;
	xx->cursor_nth  += nd;
    }

    cursor_notify(xx);
    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edReadStart(edview *xx) {
    xx->cursor_nth = 0;

    if (xx->ed->display_cutoffs) {
	if (xx->cursor_type == GT_Seq) {
	    xx->cursor_pos = 0;
	} else {
	    xx->cursor_pos = xx->contig->start;
	}
    } else {
	if (xx->cursor_type == GT_Seq) {
	    seq_t *s = get_seq(xx->io, xx->cursor_rec);

	    xx->cursor_pos = s->left-1;
	    xx->cursor_nth = 0;

	    if (sequence_get_orient(xx->io, xx->cursor_rec)) {
		s = get_seq(xx->io, xx->cursor_rec);
		xx->cursor_pos = sequence_unpadded_len(s, NULL)
		    - (s->right-1) - 1;
	    }
	} else {
	    int ustart, uend;
	    consensus_valid_range(xx->io, xx->cursor_rec, &ustart, &uend);

	    xx->cursor_pos = ustart;
	}
    }

    edSetApos(xx);

    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edReadStart2(edview *xx) {
    return edReadStart(xx);
}

int edReadEnd(edview *xx) {
    xx->cursor_nth = 0;

    if (xx->ed->display_cutoffs) {
	if (xx->cursor_type == GT_Seq) {
	    seq_t *s = get_seq(xx->io, xx->cursor_rec);
	    xx->cursor_pos = sequence_unpadded_len(s, NULL);
	} else {
	    xx->cursor_pos = xx->contig->end;
	}
    } else {
	if (xx->cursor_type == GT_Seq) {
	    seq_t *s = get_seq(xx->io, xx->cursor_rec);

	    xx->cursor_pos = s->right;

	    if (sequence_get_orient(xx->io, xx->cursor_rec)) {
		s = get_seq(xx->io, xx->cursor_rec);
		xx->cursor_pos = sequence_unpadded_len(s, NULL) - (s->left-1);
	    }
	} else {
	    int ustart, uend;
	    consensus_valid_range(xx->io, xx->cursor_rec, &ustart, &uend);

	    xx->cursor_pos = uend+1;
	}
    } 

    edSetApos(xx);

    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edReadEnd2(edview *xx) {
    return edReadEnd(xx);
}

int edContigStart(edview *xx) {
    xx->cursor_pos = xx->contig->start;
    xx->cursor_nth = 0;
    xx->cursor_type = GT_Contig;
    xx->cursor_rec = xx->cnum;
    xx->cursor_apos = xx->cursor_pos;

    cursor_notify(xx);
    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edContigEnd(edview *xx) {
    xx->cursor_pos = xx->contig->end;
    xx->cursor_nth = 0;
    xx->cursor_type = GT_Contig;
    xx->cursor_rec = xx->cnum;
    xx->cursor_apos = xx->cursor_pos;

    cursor_notify(xx);
    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

/*
 * The main editor redraw function
 */
int edview_redraw(edview *xx) {
    if (!xx->ed || xx->editorState == StateDown)
	return -1;

    if (xx->refresh_flags & ED_DISP_SEQ && xx->refresh_seq == xx->cnum)
	xx->refresh_flags |= ED_DISP_CONS;

    if (xx->displayWidth > MAX_DISPLAY_WIDTH)
	xx->displayWidth = MAX_DISPLAY_WIDTH;

#if 0
    /* Work out the status line; this may control window height */
    cur_depth = xx->status_depth;
    if (xx->refresh_flags & (ED_DISP_STATUS | ED_DISP_SCROLL)) {
	tk_redisplaySeqStatusCompute(xx, xx->displayPos, xx->displayWidth);
    }

#endif

    /* Find out what's visible */
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    /* Deal with ED_DISP_XSCROLL and ED_DISP_YSCROLL events */
    if (xx->refresh_flags & (ED_DISP_XSCROLL | ED_DISP_YSCROLL))
	tk_redisplaySeqScroll(xx, xx->r, xx->nr);

    /* Consensus emacs-style edit-status ----, -%%-, -**- */
    //tk_redisplaySeqEditStatus(xx);

    /* Redraw the consensus and/or numbers */
    if (xx->refresh_flags & (ED_DISP_CONS | ED_DISP_XSCROLL)) {
	tk_redisplaySeqConsensus(xx);
    }
    if (xx->refresh_flags & ED_DISP_RULER) {
	tk_redisplaySeqNumbers(xx);
    }

    /* Redraw the main sequences or names section */
    if (xx->refresh_flags & (ED_DISP_SEQS  | ED_DISP_SEQ |
			     ED_DISP_NAMES | ED_DISP_NAME)) {
	tk_redisplaySeqSequences(xx, xx->r, xx->nr);
    }

    /* Editor cursor position */
    if (xx->refresh_flags & ED_DISP_CURSOR) {
	tk_redisplayCursor(xx, xx->r, xx->nr);
    }

#if 0
    /* We've already computed them, but now we actually display them */
    if (xx->refresh_flags & ED_DISP_STATUS) {
	tk_redisplaySeqStatusDisplay(xx);
    }
#endif

    /* Underlining for current selection */
    if (xx->refresh_flags & ED_DISP_SELECTION) {
	redisplaySelection(xx);
    }

    if (inJoinMode(xx) && !(xx->refresh_flags & ED_DISP_NO_DIFFS))
	tk_redisplaySeqDiffs(xx);

#if 0
    /* FIXME: only need to redraw here if major change => scrolling etc */
    if (xx->refresh_flags & (ED_DISP_SEQS))
	sheet_display(&xx->ed->sw);
    if (xx->refresh_flags & (ED_DISP_NAMES))
	sheet_display(&xx->names->sw);
#endif

    xx->refresh_flags = 0;
    xx->refresh_seq = 0;

    return 0;
}

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
		       int seq_only, tg_rec *rec, int *pos, int *nth) {
    int i;
    int type = -1;
    int best_delta = INT_MAX;
    char nline[MAX_NAME_WIDTH];
    //    int exact = (name && xx->ed->stack_mode) || !name;

    if (rec) *rec = -1;
    if (pos) *pos =  0;
    if (nth) *nth =  0;

    if (!xx->screen)
	return -1;

    if (row >= xx->y_seq_start) { 
	int row2 = row + xx->displayYPos - xx->y_seq_start;
	int rid = xx->screen->rec[row2];
	int sid = -1;

	//printf("rid = %d\n", rid);
	for (; rid >= 0 && rid < xx->nr && xx->r[rid].y == row2; rid++) {
	    //printf("Check #%"PRIrec"\n", xx->r[rid].rec);
	    if (col >= xx->screen->rpos[rid])
		sid = rid;
	}
	if (sid == -1)
	    return -1;

	//printf("Row %d Col %d Seq = %d/%"PRIrec" start %d\n",
	//       row, col, sid, xx->r[sid].rec, xx->screen->rpos[sid]);

	if (rec)
	    *rec = xx->r[sid].rec;
	if (pos) {
	    if (xx->screen->rpos[sid] < 0) {
		int sp, sn, rp, rn;
		seq_t *s = cache_search(xx->io, GT_Seq, xx->r[sid].rec);

		//int upos, nth;
		//sequence_get_upos(xx->io, &s, -xx->screen->rpos[sid],
		//		  &upos, &nth);
		//*pos = upos + xx->screen->pos[col] - xx->screen->pos[0];

		// Try 2
		sequence_cigar2pos(xx->io, s,
				   xx->screen->cigar_ind[sid],
				   xx->screen->cigar_len[sid],
				   &sp, &sn, &rp, &rn);
		sp--;
		rp--;
		*pos = rp + xx->screen->pos[col] - xx->screen->pos[0];
		
	    } else {
		*pos = 
		    xx->screen->pos[col] - xx->screen->pos[xx->screen->rpos[sid]];
	    }
	}
	if (nth)
	    *nth = xx->screen->nth[col];

	return GT_Seq;
    }

    /* Special case - the reserve row numbers */
    if (row == xx->y_cons) {
	*rec = xx->cnum;
	*pos = xx->screen->pos[col];
	*nth = xx->screen->nth[col];
	type = GT_Contig;

	if (xx->ed->hide_annos || seq_only)
	    return type;

	/* Look for consensus tags */
	for (i = 0; i < xx->nr; i++) {
	    if (xx->r[i].y != -1)
		break;

	    if ((xx->r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)
		continue;

	    if (col + xx->displayPos >= xx->r[i].start &&
		col + xx->displayPos <= xx->r[i].end) {
		*rec = xx->r[i].rec;
		*pos = col + xx->displayPos - xx->r[i].start;
		type = GT_AnnoEle;
	    }
	}

	return type;
    }

    if (row < xx->y_seq_start)
	return -1;

    return -1;
#if 0    
    /* A sequence, so find out what's visible */
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    /* Inefficient, but just a copy from tk_redisplaySeqSequences() */
    i = edview_binary_search_y(xx->r, xx->nr, xx->displayYPos);
    memset(nline, ' ', MAX_NAME_WIDTH);
    for (; i < xx->nr; i++) {
	if ((xx->ed->hide_annos || seq_only || name) &&
	    ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO))
	    continue;

#ifndef CACHED_CONS_VISIBLE
	if ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISCONS)
	    continue;
#endif

	if (xx->r[i].y + xx->y_seq_start - xx->displayYPos == row) {
	    int delta;

	    /* Find distance from object, in X */
	    if (xx->ed->stack_mode && name) {
		/* In names display during stacking mode */
		int p1 = xx->r[i].start - xx->displayPos;
		int p2 = xx->r[i].end   - xx->displayPos;
		double nc = xx->names->sw.columns;

		if (p1 < 0) p1 = 0;
		p1 = p1 * (nc / xx->displayWidth);
		if (p2 < 0) p2 = 0;
		p2 = p2 * (nc / xx->displayWidth);

		while (p1 < nc && nline[p1] != ' ')
		    p1++;

		if (col >= p1 && (col < p2 || col == p1))
		    delta = 0;
		else
		    delta = INT_MAX;

		if (p2 > nc)
		    p2 = nc;

		do {
		    nline[p1++] = '.';
		} while (p1 < p2);

	    } else {
		/* In sequence display, or only 1 seq per line */
		if (col + xx->displayPos < xx->r[i].start)
		    delta = xx->r[i].start - (col + xx->displayPos);
		else if (col + xx->displayPos > xx->r[i].end)
		    delta = col + xx->displayPos - xx->r[i].end;
		else {
		    delta = 0;
		}
	    }

	    /* And if this is closest match, use it */
	    if (best_delta >= delta) {
		best_delta =  delta;
		*rec = xx->r[i].rec;
		*pos = col + xx->displayPos - xx->r[i].start;
                type = (xx->r[i].flags & GRANGE_FLAG_ISMASK)
                       == GRANGE_FLAG_ISANNO
		    ? GT_AnnoEle
		    : GT_Seq;
	    }
	}
    }

    return !exact || best_delta == 0 ? type : -1;
#endif
}

/*
 * More or less the opposite of the above function. Given a record and type
 * this returns the row number on screen. If non-NULL xmin and xmax are the
 * X coordinates of the extends of this record. These may be beyond the bounds
 * of the window.
 *
 * Returns Y coordinate (and optionally min/max X coordinates) if found
 *        -1 if not (with xmin/xmax unset).
 */
int edview_row_for_item(edview *xx, tg_rec rec, int *xmin, int *xmax) {
    int i, r = -1;
    HacheItem *hi;

    if (rec == xx->cnum) {
	if (xmin) *xmin = -xx->displayPos;
	if (xmax) *xmax = -xx->displayPos;
	return 0;
    }

    /* A sequence, so find out what's visible */
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    /* And search for rec in this list */
    if (!xx->rec_hash)
	return -1;
    hi = HacheTableSearch(xx->rec_hash, (char *)&rec, sizeof(rec));
    if (!hi)
	return -1;

    i = hi->data.i;
    if (xmin) *xmin = xx->r[i].start - xx->displayPos;
    if (xmax) *xmax = xx->r[i].end   - xx->displayPos;
    r = xx->r[i].y + xx->y_seq_start - xx->displayYPos;
    
    return r >= xx->y_seq_start ? r : -1;
}


int inJoinMode(edview *xx) {
    return xx->link ? 1 : 0;
}

void edDisplayTrace(edview *xx) {
    seq_t *s;

    if (xx->cursor_type == GT_Seq) {
	/* Single sequence */
	s = get_seq(xx->io, xx->cursor_rec);
	tman_manage_trace("ANY", sequence_get_name(&s), xx->cursor_pos,
			  0, 0, /* left/right clips */
			  s->len < 0, /* complemented */
			  1, /* base spacing */
			  sequence_get_name(&s),
			  xx, xx->cursor_rec, 0, 0);
    } else if (xx->cursor_type == GT_Contig) {
	/* Consensus click */
	rangec_t *r;
	int nr, i;

	/* Shut down existing traces */
	tman_shutdown_traces(xx, 2);

	/* And display the new ones */
	puts("FIXME: reuse existing cache of items");
	r = contig_seqs_in_range(xx->io, &xx->contig,
				 xx->cursor_apos, xx->cursor_apos,
				 CSIR_SORT_BY_X, &nr);

	for (i = 0; i < nr; i++) {
	    s = get_seq(xx->io, r[i].rec);
	    tman_manage_trace("ANY", sequence_get_name(&s), xx->cursor_pos,
			      0, 0, /* left/right clips */
			      s->len < 0, /* complemented */
			      1, /* base spacing */
			      sequence_get_name(&s),
			      xx, r[i].rec, 0, 0);
	}
	free(r);
    }

    tman_reposition_traces(xx, xx->cursor_apos, 0);
}

/*
 * Given a sequence record number this identifies all other sequence
 * records from the same template. The returned array is malloced and should
 * be freed by the caller once finished with.
 *
 * Returns pointer to array of records of size *nrec on success
 *         NULL on failure (or zero found)
 */
tg_rec *edGetTemplateReads(edview *xx, tg_rec seqrec, int *nrec) {
    seq_t *s = get_seq(xx->io, seqrec);
    tg_rec *r = NULL, p;

    if (!s)
	return NULL;

    /* FIXME: support s->template_rec */

    /* Solexa data is easy: we have just one other end */
    p = sequence_get_pair(xx->io, s);
    if (p > 0) {
	*nrec = 1;
	r = malloc(sizeof(*r));
	*r = p;
    } else {
	*nrec = 0;
    }

    return r;
}


/* ---------------------------------------------------------------------- */
/* Selection aka cut'n'paste handling code */

/*
 * (Un)draws the selection - toggles it on or off
 */
static void toggle_select(edview *xx, tg_rec seq,
			  int from_pos, int from_nth,
			  int to_pos, int to_nth) {
    int row, xmin, xmax;

    if (from_pos > to_pos || (from_pos == to_pos && from_nth > to_nth)) {
	int temp;
	temp = from_pos; from_pos = to_pos; to_pos = temp;
	temp = from_nth; from_nth = to_nth; to_nth = temp;
    }

    /* Find out the X and Y coord, and exit now if it's not visible */
    if (-1 == (row = edview_row_for_item(xx, seq, &xmin, NULL)))
	return;

    if (from_pos < 0) {from_pos = 0; from_nth = 0;}
    if (to_pos   < 0) {to_pos   = 0; to_nth   = 0;}

    if (seq != xx->cnum) {
	HacheItem *hi;
	int i, p;
	tg_rec cnum;
	int cpos, orient;
	sequence_get_position(xx->io, seq, &cnum, &cpos, NULL, &orient);

	from_pos += cpos;
	to_pos += cpos;

	if (!(hi = HacheTableSearch(xx->rec_hash, (char *)&seq, sizeof(seq))))
	    return;
	assert(xx->r[hi->data.i].rec == seq);
	i = hi->data.i;

	/* xmin */
	p = xx->screen->rpos[i];
	if (p < 0) p = 0;

	while (p < xx->displayWidth) {
	    int p2 = xx->screen->pos[p];
	    if (p2 > from_pos) {
		/* Off left edge => xmin=0 */
		break;
	    }
	    if (p2 == from_pos && xx->screen->nth[p] == from_nth)
		break;
	    p++;
	}
	xmin = p;

	/* xmax */
	while (p < xx->displayWidth) {
	    int p2 = xx->screen->pos[p];
	    if (p2 > to_pos) {
		/* Off left edge => entirely invisible */
		return;
	    }
	    if (p2 == to_pos && xx->screen->nth[p] == to_nth)
		break;
	    p++;
	}
	xmax = p;
    } else {
	/* Consensus is selected */
	int p;

	xmin = 0; xmax = xx->displayWidth-1;
	for (p = 0; p < xx->displayWidth; p++) {
	    int p2 = xx->screen->pos[p];
	    if (p2 > from_pos)
		break;
	    if (p2 == from_pos && xx->screen->nth[p] == from_nth)
		break;
	}
	xmin = p;

	for (; p < xx->displayWidth; p++) {
	    int p2 = xx->screen->pos[p];
	    if (p2 > to_pos)
		return; /* off left edge */
	    if (p2 == to_pos && xx->screen->nth[p] == to_nth)
		break;
	}
	xmax = p;
    }
    
    /* Convert xmin/xmax to the region we wish to view */
    //    xmin += from_pos;
    //    xmax = xmin + to_pos - from_pos;

    /* clip to screen */
    if (xmin < 0) xmin = 0;
    if (xmax >= xx->displayWidth) xmax = xx->displayWidth-1;

    if (xmin > xmax)
	return;

    /* Toggle the line */
    XawSheetOpHilightText(&xx->ed->sw, xmin, row, xmax-xmin+1,
			  sh_select, HOP_TOG);
}

void redisplaySelection(edview *xx) {
    toggle_select(xx, xx->select_seq,
		  xx->select_start, xx->select_start_nth,
		  xx->select_end,   xx->select_end_nth);
}

/*
 * Callback from Tk_OwnSelection().
 */
static void edSelectionLost(ClientData cd) {
    edview *xx = (edview *)cd;

    /* Undisplay the selection */
    redisplaySelection(xx);

    xx->select_made = 0;
    xx->select_seq = 0;
    xx->select_start = 0;
    xx->select_start_nth = 0;
    xx->select_end = 0;
    xx->select_end_nth = 0;
}

int edSelectClear(edview *xx) {
    if (xx->select_made && EDTKWIN(xx->ed))
	Tk_ClearSelection(EDTKWIN(xx->ed), XA_PRIMARY);
    edSelectionLost((ClientData)xx);

    return 0;
}

void edSelectFrom(edview *xx, int pos, int nth) {
    /* Undisplay an old selection */
    if (xx->select_made)
	redisplaySelection(xx);
    else
	xx->select_made = 1;

    /* Set start/end */
    xx->select_seq = xx->cursor_rec;
    if (xx->select_seq != xx->cnum) {
 	tg_rec cnum;
	tg_rec key;
	HacheItem *hi;
	int i, p;
	
	key = xx->cursor_rec;
	if (!xx->rec_hash || !xx->r)
	    return;

	if (!(hi = HacheTableSearch(xx->rec_hash, (char *)&key, sizeof(key))))
	    return;
	assert(xx->r[hi->data.i].rec == xx->cursor_rec);

	i = hi->data.i;
	if ((p = xx->screen->rpos[i]) >= 0) {
	    /* p chars in from the left screen edge */
	    xx->select_start     = xx->screen->pos[pos] - xx->screen->pos[p];
	    xx->select_start_nth = xx->screen->nth[pos] - xx->screen->nth[p];
	} else {
	    /* -p bases into the raw sequence */
	    int sp, sn, rp, rn;
	    seq_t *s = get_seq(xx->io, xx->select_seq);
	    sequence_cigar2pos(xx->io, s,
			       xx->screen->cigar_ind[i],
			       xx->screen->cigar_len[i],
			       &sp, &sn, &rp, &rn);
	    sp--;
	    rp--;

	    xx->select_start     = rp + xx->screen->pos[pos] - 
		xx->screen->pos[0];
	    xx->select_start_nth = rn + xx->screen->nth[pos] -
		xx->screen->nth[0];
	}

	Tk_OwnSelection(EDTKWIN(xx->ed), XA_PRIMARY, edSelectionLost,
			(ClientData)xx);

	redisplaySelection(xx);
	return;
#if 0
	tg_rec cnum;
	int cpos, left, right, orient;
	seq_t *s = get_seq(xx->io, xx->select_seq);

	cache_incr(xx->io, s);
	sequence_get_position(xx->io, xx->select_seq,
			      &cnum, &cpos, NULL, &orient);
	pos -= cpos;

	if (xx->ed->display_cutoffs) {
	    left  = 0;
	    right = sequence_unpadded_len(s, NULL);
	} else {
	    if ((s->len < 0) ^ orient) {
		int ulen = sequence_unpadded_len(s, NULL);
		left  = ulen - (s->right-1) - 1;
		right = ulen - (s->left-1);
	    } else {
		left  = s->left - 1;
		right = s->right;
	    }
	}

	if (pos < left)
	    pos = left;
	if (pos > right+1)
	    pos = right+1;

	cache_decr(xx->io, s);
#endif
    } else {
	if (pos >= 0 && pos <= xx->displayWidth)
	    pos = xx->screen->pos[pos];
	else
	    pos += xx->displayPos; /* Random guess! */
    }
    xx->select_start = xx->select_end = pos;

    Tk_OwnSelection(EDTKWIN(xx->ed), XA_PRIMARY, edSelectionLost,
		    (ClientData)xx);

    /* Display new selection */
    redisplaySelection(xx);
}

void edSelectTo(edview *xx, int pos, int nth) {
    if (!xx->select_made) {
	edSelectFrom(xx, pos, nth);
    }

    /* Undisplay old selection */
    redisplaySelection(xx);

    /* Set start/end */
    if (xx->select_seq != xx->cnum) {
 	tg_rec cnum;
	tg_rec key;
	HacheItem *hi;
	int i, p, ulen, unth;
	seq_t *s = get_seq(xx->io, xx->select_seq);

	key = xx->cursor_rec;
	if (!xx->rec_hash || !xx->r)
	    return;

	if (!(hi = HacheTableSearch(xx->rec_hash, (char *)&key, sizeof(key))))
	    return;
	assert(xx->r[hi->data.i].rec == xx->cursor_rec);

	i = hi->data.i;
	if ((p = xx->screen->rpos[i]) >= 0) {
	    /* p chars in from the left screen edge */
	    xx->select_end     = xx->screen->pos[pos] - xx->screen->pos[p];
	    xx->select_end_nth = xx->screen->nth[pos] - xx->screen->nth[p];
	} else {
	    /* -p bases into the raw sequence */
	    int sp, sn, rp, rn;
	    sequence_cigar2pos(xx->io, s,
			       xx->screen->cigar_ind[i],
			       xx->screen->cigar_len[i],
			       &sp, &sn, &rp, &rn);
	    sp--;
	    rp--;

	    xx->select_end     = rp + xx->screen->pos[pos]-xx->screen->pos[0];
	    xx->select_end_nth = rn + xx->screen->nth[pos]-xx->screen->nth[0];
	}

	if (xx->select_end >= (ulen = sequence_unpadded_len(s, &unth))) {
	    xx->select_end = ulen-1;
	    xx->select_end_nth = unth;
	}

	redisplaySelection(xx);
	return;

#if 0
	int cpos, left, right, orient;
	seq_t *s = get_seq(xx->io, xx->select_seq);

	cache_incr(xx->io, s);
	sequence_get_position(xx->io, xx->select_seq,
			      &cnum, &cpos, NULL, &orient);
	pos -= cpos;

	if (xx->ed->display_cutoffs) {
	    left  = 0;
	    right = sequence_unpadded_len(s, NULL);
	} else {
	    if ((s->len < 0) ^ orient) {
		int ulen = sequence_unpadded_len(s, NULL);
		left  = ulen - (s->right-1) - 1;
		right = ulen - (s->left-1);
	    } else {
		left  = s->left - 1;
		right = s->right;
	    }
	}

	if (pos < left)
	    pos = left;
	if (pos > right-1)
	    pos = right-1;

	cache_decr(xx->io, s);
#endif
    } else {
	if (pos >= 0 && pos <= xx->displayWidth)
	    pos = xx->screen->pos[pos];
	else
	    pos += xx->displayPos;
    }
    xx->select_end = pos;

    /* Display new selection */
    redisplaySelection(xx);
}

void edSelectSet(edview *xx, tg_rec rec, int start, int start_nth, 
		 int end, int end_nth) {
    int do_x = 0;

    /* Undisplay an old selection */
    if (xx->select_made)
	redisplaySelection(xx);

    xx->select_made = 0;


    xx->select_seq       = rec;
    xx->select_start     = start;
    xx->select_start_nth = start_nth;
    xx->select_end       = end;
    xx->select_end_nth   = end_nth;
    xx->select_made      = 1;

    /* Scroll and redraw if appropriate */
    if (xx->select_end+2 >= xx->displayPos + xx->displayWidth) {
	set_displayPos(xx, xx->select_end+2 - xx->displayWidth);
	do_x = 1;
    }
    if (xx->select_start-1 <= xx->displayPos) {
	set_displayPos(xx, xx->select_start-1);
	do_x = 1;
    }

    if (do_x) {
	xx->refresh_flags = ED_DISP_ALL;
	ed_set_xslider_pos(xx, xx->displayPos);
    }

    xx->refresh_flags |= ED_DISP_SELECTION;
    edview_redraw(xx);
}

/*
 * Automatically called when X wishes to obtain a selection. We register
 * this procedure in the initialise code of tkEditor.c.
 *
 * Return codes expected:
 *    -1  Failure
 *    >0  Number of bytes
 */
int edGetSelection(ClientData clientData, int offset, char *buffer,
		   int bufsize) {
    Editor *ed = (Editor *)clientData;
    int start, end, len, st_exists, en_exists;
    edview *xx = ed->xx;
    seq_t *s, *sorig;
    int bin_comp = 0;
    HacheItem *hi;

    /* Do we have a selection? */
    if (!xx->select_made)
	return -1;

    if (xx->select_seq == xx->cnum) {
	/* Swap if start..range is backwards */
	if ((xx->select_start     > xx->select_end) || 
	    (xx->select_start    == xx->select_end &&
	     xx->select_start_nth > xx->select_end_nth)) {
	    int tmp;
	    tmp = xx->select_start;
	    xx->select_start = xx->select_end;
	    xx->select_end = tmp;
	    tmp = xx->select_start_nth;
	    xx->select_start_nth = xx->select_end_nth;
	    xx->select_end_nth = tmp;
	}
	start = xx->select_start;
	len = xx->select_end - xx->select_start + 1;
	calculate_consensus_simple(xx->io, xx->cnum,
				   xx->select_start, //xx->select_start_nth,
				   xx->select_end,   //xx->select_end_nth
				   buffer, NULL);
    } else {
	if (!(s = get_seq(xx->io, xx->select_seq)))
	    return -1;
	cache_incr(xx->io, s);
	start = sequence_get_spos(xx->io, &s,
				  xx->select_start, xx->select_start_nth,
				  &st_exists);
	end   = sequence_get_spos(xx->io, &s,
				  xx->select_end,   xx->select_end_nth,
				  &en_exists);
	start += offset;
	end   += offset;
	cache_decr(xx->io, s);
    
	hi = HacheTableSearch(xx->rec_hash,
			      (char *)&xx->select_seq,
			      sizeof(xx->select_seq));
	if (hi) {
	    assert(xx->r[hi->data.i].rec == xx->select_seq);
	    bin_comp = xx->r[hi->data.i].comp;
	}

	if (start == end && !st_exists && !en_exists)
	    return 0;

	if ((xx->select_start > xx->select_end ||
	     (xx->select_start == xx->select_end &&
	      xx->select_start_nth > xx->select_end_nth))
	    ^ (s->len < 0)
	    ^ bin_comp) {
	    if (!st_exists) start--;
	} else {
	    if (!st_exists) start++;
	}

	if (start > end) {
	    len = start;
	    start = end;
	    end = len;
	}

	len = end - start+1 > bufsize ? bufsize : end - start + 1;

	sorig = s = get_seq(xx->io, xx->select_seq);
	if (s->len < 0) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	    start = s->len - end-1;
	}
	memcpy(buffer, s->seq+start, len);
	if (s != sorig) {
	    free(s);
	}

	/* Reverse comp if bin is complemented too */
	if (bin_comp) {
	    complement_seq(buffer, len);
	}
    }

    return len;
}
