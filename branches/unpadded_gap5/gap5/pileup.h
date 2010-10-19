#ifndef _PILEUP_H_
#define _PILEUP_H_

#include <tg_gio.h>

typedef struct pileup {
    struct pileup *next;  // A link list, for active seqs

    rangec_t *r;          // The range and seq associated with struct
    seq_t *s;

    int pos;              // Current unpadded position in seq
    int nth;		  // nth base at unpadded position 'pos'
    int seq_offset;       // Current base position in s->seq[] array.

    int cigar_ind;        // Current location in s->alignment cigar str
    int cigar_op;         // Current cigar operation
    int cigar_len;        // Remaining length of this cigar op

    int  qual;            // Current qual (for active seq only)
    char base;		  // Current base (for active seq only)
    char sclip;		  // True if soft-clipped
    //char last_sclip;	  // True if last base emitted was soft-clipped
    char start;		  // True is this is a new sequence
} pileup_t;

/*
 * Loops through a set of supplied ranges producing columns of data.
 * When found, it calls func with clientdata as a callback. Func should
 * return 0 for success and non-zero for failure.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int pileup_loop(GapIO *io, rangec_t *r, int nr,
		int start, int start_nth, int end, int end_nth,
		int (*func)(void *client_data,
			    pileup_t *p,
			    int pos,
			    int nth),
		void *client_data);

#endif
