#include "pileup.h"

/*
 * Fetches the next base => the nth base at unpadded position pos. (Nth can
 * be greater than 0 if we have an insertion in this column). Do not call this
 * with pos/nth lower than the previous query, although higher is better.
 * (This allows it to be initialised at base 0.)
 *
 * Stores the result in base and also updates is_insert to indicate that
 * this sequence still has more bases in this position beyond the current
 * nth parameter.
 *
 * Returns 1 if a base was fetched
 *         0 if not (eg ran off the end of sequence)
 */
static int get_next_base(pileup_t *p, int pos, int nth, int of,
			 int *is_insert) {
    seq_t *s = p->s;
    int comp = (s->len < 0) ^ p->r->comp;
    int more = 1;
    int justify = 0;
    int ilen;

    p->base = '?';
    *is_insert = 0;

#if 1
    /* Find pos first */
    while (p->pos < pos) {
	p->nth = 0;

	if (p->cigar_len == 0) {
	    if (p->cigar_ind >= s->alignment_len)
		return 0;

	    if (comp) {
		int r_ind = s->alignment_len - p->cigar_ind - 2;
		p->cigar_len = s->alignment[r_ind];
		p->cigar_op  = s->alignment[r_ind+1];
	    } else {
		p->cigar_len = s->alignment[p->cigar_ind];
		p->cigar_op  = s->alignment[p->cigar_ind+1];
	    }
	    p->cigar_ind += 2;
	}

	switch (p->cigar_op) {
	case 'H':
	    p->cigar_len = 0;
	    break;

	case 'S':
	case 'M':
	    p->pos++;
	    p->seq_offset++;
	    p->cigar_len--;
	    break;

	case 'D':
	    p->pos++;
	    p->cigar_len--;
	    break;

	case 'P':
	    p->cigar_len = 0; /* An insert that doesn't consume seq */
	    break;

	case 'I':
	    p->seq_offset += p->cigar_len;
	    p->cigar_len = 0;
	    break;
	}
    }

    if (p->pos < pos)
	return 0;

//    /* Compute total size of insert due to this seq */
//    ilen = 0;
//    if (comp && nth /*&& p->in_size == 0*/) {
//	int tmp_ind = p->cigar_ind;
//	int tmp_op  = p->cigar_op;
//	int r_ind;
//
//	if (p->cigar_len == 0 && tmp_ind < s->alignment_len) {
//	    r_ind = s->alignment_len - tmp_ind - 2;
//	    tmp_op  = s->alignment[r_ind+1];
//	    tmp_ind += 2;
//	}
//
//	while (tmp_ind < s->alignment_len) {
//	    if (tmp_op == 'I' || tmp_op == 'P') {
//		r_ind = s->alignment_len - tmp_ind - 2;
//		ilen += s->alignment[r_ind+2];
//		tmp_op  = s->alignment[r_ind+1];
//		tmp_ind += 2;
//	    } else {
//		break;
//	    }
//	}
//    }

    /* Now at pos, find nth base */
    while (p->nth < nth && more) {
	if (p->cigar_len == 0) {
	    if (p->cigar_ind >= s->alignment_len)
		return 0; /* off end of seq */

	    if (comp) {
		int r_ind = s->alignment_len - p->cigar_ind - 2;
		p->cigar_len = s->alignment[r_ind];
		p->cigar_op  = s->alignment[r_ind+1];
	    } else {
		p->cigar_len = s->alignment[p->cigar_ind];
		p->cigar_op  = s->alignment[p->cigar_ind+1];
	    }
	    p->cigar_ind += 2;
	}

	switch (p->cigar_op) {
	case 'H':
	    p->cigar_len = 0;
	    break;

	case 'S':
	case 'M':
	case 'D':
	    more = 0;
	    break;

	case 'P':
	    p->cigar_len--;
	    p->nth++;
	    more = 1;
	    break;

	case 'I':
//	    if (comp && of - p->nth > ilen) {
//		p->nth++;
//		more = 0;
//		justify = 1;
//		break;
//	    }
	    p->cigar_len--;
	    p->nth++;
	    p->seq_offset++;
	    more = 1;
	    break;
	}
    }

    if (justify || (p->nth < nth && p->cigar_op != 'I')) {
	p->base = '-';
	p->qual = 0;
    } else {
	if (p->cigar_op == 'D') {
	    p->base = '*';
	    p->qual = 0;
	} else if (p->cigar_op == 'P') {
	    p->base = '+';
	    p->qual = 0;
	} else {
	    if (comp) {
		p->qual = s->conf[ABS(s->len) - (p->seq_offset+1)];
		p->base = complement_base(s->seq[ABS(s->len) - (p->seq_offset+1)]);
	    } else{
		p->qual = s->conf[p->seq_offset];
		p->base = s->seq[p->seq_offset];
	    }
	}
    }

    if (p->cigar_len == 0 && p->cigar_ind < s->alignment_len) {
	/*
	 * Query next operation to check for 'I', but don't modify our
	 * p->cigar_ fields as the callback function may need to query
	 * these to know where we are during CIGAR processing.
	 */

	int tmp_ind = p->cigar_ind;
	int tmp_len = p->cigar_len;
	int tmp_op  = p->cigar_op;

	do {
	    int r_ind;

	    r_ind = comp
		? s->alignment_len - tmp_ind - 2
		: tmp_ind;
	    tmp_len = s->alignment[r_ind];
	    tmp_op  = s->alignment[r_ind+1];

	    if (tmp_op == 'I' || tmp_op == 'P') {
		*is_insert += tmp_len;
	    } else {
		break;
	    }

	    tmp_ind += 2;
	} while (tmp_ind < s->alignment_len);

//	int tmp_len, tmp_op;
//	if (comp) {
//	    int r_ind = s->alignment_len - p->cigar_ind - 2;
//	    tmp_len = s->alignment[r_ind];
//	    tmp_op  = s->alignment[r_ind+1];
//	} else {
//	    tmp_len = s->alignment[p->cigar_ind];
//	    tmp_op  = s->alignment[p->cigar_ind+1];
//	}
//	if (tmp_op == 'I' || tmp_op == 'P')
//	    *is_insert = tmp_len;
    } else {
	int tmp_ind = p->cigar_ind;
	int tmp_len = p->cigar_len;
	int tmp_op  = p->cigar_op;

	while (tmp_op == 'I' || tmp_op == 'P') {
	    int r_ind;

	    *is_insert += tmp_len;
	    tmp_ind += 2;
	    if (tmp_ind >= s->alignment_len)
		break;

	    r_ind = comp
		? s->alignment_len - tmp_ind - 2
		: tmp_ind;
	    tmp_len = s->alignment[r_ind];
	    tmp_op  = s->alignment[r_ind+1];
	}
    }
    
    if (comp) {
	p->sclip = (p->seq_offset+1 <= ABS(s->len) - s->right ||
		    p->seq_offset   >  ABS(s->len) - s->left);
    } else {
	p->sclip = (p->seq_offset+1 < s->left || p->seq_offset >= s->right);
    }

    return 1;

#else

    /* May need to catch up to pos on the first time through */
    for (;;) {
	int tmp_pos;

	/* Fetch next base */
	if (p->cigar_len == 0) {
	    if (p->cigar_ind >= s->alignment_len)
		return 0;

	    if (comp) {
		int r_ind = s->alignment_len - p->cigar_ind - 2;
		p->cigar_len = s->alignment[r_ind];
		p->cigar_op  = s->alignment[r_ind+1];
	    } else {
		p->cigar_len = s->alignment[p->cigar_ind];
		p->cigar_op  = s->alignment[p->cigar_ind+1];
	    }
	    p->cigar_ind += 2;
	}

	/* Pos is an unpadded position, nth is the nth base at that pos */
	tmp_pos = p->pos;
	switch (p->cigar_op) {
	case 'S':
	case 'M':
	case 'P':
	    if (nth) {
		p->sclip = p->cigar_op == 'S' ? 1 : p->last_sclip;
		p->base = '-';
		nth--;
	    } else {
		p->sclip = p->cigar_op == 'S' ? 1 : 0;
		p->last_sclip = p->sclip;
		p->pos++;
		p->cigar_len--;
		p->qual = comp
		    ? s->conf[ABS(s->len) - (p->seq_offset+1)]
		    : s->conf[p->seq_offset];
		p->base = comp
		    ? complement_base(s->seq[ABS(s->len) - ++p->seq_offset])
		    : s->seq[p->seq_offset++];
	    }
	    break;

	case 'D':
	    p->sclip = p->last_sclip;
	    if (nth) {
		p->base = '-';
		p->qual = 0;
		nth--;
	    } else {
		p->pos++;
		p->cigar_len--;
		p->base = '*';
		p->qual = 0;
	    }
	    break;

	case 'I':
	    p->last_sclip = 0;
	    p->cigar_len--;
	    p->qual = comp
		? s->conf[ABS(s->len) - (p->seq_offset+1)]
		    : s->conf[p->seq_offset];
	    p->base = comp 
		? complement_base(s->seq[ABS(s->len) - ++p->seq_offset])
		: s->seq[p->seq_offset++];
	    break;

	case 'H':
	    p->cigar_len = 0; /* just skip it */
	    break;

	default:
	    fprintf(stderr, "Unrecognised cigar operation: %c\n",
		    p->cigar_op);
	    p->cigar_len = 0; /* just skip it */
	}

	/* indel next? */
	if (comp) {
	    int r_ind = s->alignment_len - p->cigar_ind - 2;
	    if (tmp_pos >= pos &&
		p->cigar_len == 0 &&
		p->cigar_ind < s->alignment_len &&
		(s->alignment[r_ind+1] == 'I' ||
		 s->alignment[r_ind+1] == 'P')) {
		p->cigar_len = s->alignment[r_ind];
		p->cigar_op  = 'I';
		p->cigar_ind += 2;

		*is_insert = p->cigar_len;
		return 1;
	    }
	} else {
	    if (tmp_pos >= pos &&
		p->cigar_len == 0 &&
		p->cigar_ind < s->alignment_len &&
		(s->alignment[p->cigar_ind+1] == 'I' ||
		 s->alignment[p->cigar_ind+1] == 'P') {
		p->cigar_len = s->alignment[p->cigar_ind];
		p->cigar_op  = 'I';
		p->cigar_ind += 2;

		*is_insert = p->cigar_len;
		return 1;
	    }
	}

	if (tmp_pos >= pos) {
	    *is_insert = (p->cigar_op == 'I' || p->cigar_op == 'P')
		? p->cigar_len : 0;
	    return 1;
	}
    }
#endif

    return 0;
}

static int pileup_sort(const void *a, const void *b) {
    const pileup_t *p1 = (const pileup_t *)a;
    const pileup_t *p2 = (const pileup_t *)b;

    return p1->r->start - p2->r->start;
}

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
		void *client_data) {
    int i, j, ret = -1;
    pileup_t *parray, *p, *next, *active = NULL;
    int is_insert, nth = 0, of = 0;
    int col = 0;

    /*
     * Allocate an array of parray objects, one per sequence, 
     * and sort them into positional order.
     *
     * While processing, we also thread a linked list through the parray[]
     * to keep track of current "active" sequences (ones that overlap
     * our X column).
     */
    parray = malloc(nr * sizeof(*parray));
    for (i = j = 0; i < nr; i++) {
	if ((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	    parray[j].r = &r[i];
	    parray[j].s = cache_search(io, GT_Seq, r[i].rec);
	    cache_incr(io, parray[j].s);
	    j++;
	}
    }
    nr = j;
    qsort(parray, nr, sizeof(*parray), pileup_sort);

    /*
     * Loop through all seqs in range producing padded versions.
     * For each column in start..end
     *   1. Add any new sequences to our 'active' linked list through parray.
     *   2. Foreach active sequence, fetch the next base and an indication
     *      of whether more bases are at this column.
     *   3. If at sequence EOF, remove from 'active' linked list, otherwise
     *      display the base.
     *   4. Loop. We may stall and loop multiple times on the same column
     *      while we've determined there are still more bases at this
     *      column (ie an insert).
     */
    j = 0;
    i = start;
    nth = start_nth;
    while (i <= end || (i == end && nth <= end_nth)) {
	int v;
	pileup_t *last;

	/* Add new seqs */
	/* FIXME: optimise this in case of very long seqs. May want to
	 * process CIGAR and jump ahead
	 */
	while (j < nr && parray[j].r->start <= i) {
	    if ((parray[j].r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
		p = &parray[j];
		p->start = 1;
		p->next = active;
		active = p;
		p->pos = parray[j].r->start-1;
		p->cigar_len = 0;
		p->cigar_ind = 0;
		p->cigar_op  = 'X';
		p->seq_offset = -1;
		//p->last_sclip = 1;
	    }
	    j++;
	}


	/* Pileup */
	is_insert = 0;
	for (p = active, last = NULL; p; p = next) {
	    int ins;

	    next = p->next;

	    if (get_next_base(p, i, nth, of, &ins)) {
		last = p;
	    } else {
		/* Remove sequence */
		if (last) {
		    last->next = p->next;
		} else {
		    active = p->next;
		}
	    }

	    if (is_insert < ins)
		is_insert = ins;
	}

	/* active is now a linked list, so pass into the callback */
	v = func(client_data, active, i, nth);
	if (v == 1)
	    break; /* early abort */

	if (v != 0)
	    goto error;


	if (is_insert) {
	    nth++;
	    if (of < is_insert)
		of = is_insert;
	} else {
	    nth = 0;
	    of = 0;
	    i++;
	}
    }

    ret = 0;
 error:

    /* Tidy up */
    for (i = 0; i < nr; i++) {
	cache_decr(io, parray[i].s);
	//free(parray[i].s);
    }
    free(parray);

    return ret;
}
