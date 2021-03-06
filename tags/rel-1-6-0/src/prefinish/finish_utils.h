#ifndef _FINISH_UTILS_H
#define _FINISH_UTILS_H

#include "IO.h"

int finish_next_expt_id(int reset);
int finish_next_group_id(int reset);

/*
 * Computes the average expected length for a sequence with position 'start'
 * to 'end' on a template with minimum and maximum size ranges from
 * 't_start1'..'t_end1' and 't_start2'..'t_end2'.
 * All ranges are inclusive.
 *
 * When the sequence lies entirely within the minimum template size then
 * the expected length is the length we've specified (start to end). Otherwise
 * there is a chance that the template may be short (> min, but < max) and so
 * the sequence may be terminated prematurely. We cannot possibly know when
 * this will happen, but we can guess at an average case.
 *
 * Returns:
 *	An average length as the return value. new_start and new_end are
 *	filled in with the new expected start and end values, which also
 *	takes into account template 'clipping'.
 */
int finish_avg_length(int start, int end, int dir,
		      int t_start1, int t_end1,
		      int t_start2, int t_end2,
		      int *new_start, int *new_end);


/*
 * Computes the chance (0 to 1) of all the sequence at position start..end
 * existing on a template spanning between t_start1..t_end1 (minimum) and
 * t_start2..t_end2 (maximum).
 *
 * Returns chance (1 == always, 0 == never, plus between values)
 */
double template_exists_chance(int start, int end,
			      int t_start1, int t_end1,
			      int t_start2, int t_end2);


/*
 * complement_seq_qual_mapping
 *
 * Complements a sequence, quality buffer, and integer array (which is
 * used externally for mapping unpadded to padded positions, but for the
 * purposes of this code it could contain anything).
 *
 * Arguments:
 *	len		Dimension of seq, qual and map arrays
 *	seq		The DNA sequence (does not need to be nul terminated)
 *	qual		The confidence values
 *	map		The integer mapping array.
 */
void complement_seq_qual_mapping(int len, char *seq, float *qual, int *map);

/*
 * Looks for cloning vector. It's useful to know if its there so that we
 * don't allow extending into vector. If svec_also is true then this also
 * treats SVEC as the same as CVEC. The reason for this is to treat sequencing
 * vector as the clone-end for EST projects.
 *
 * Returns:
 *	Updates contents of left and right to be 0 for no vector and 1 for
 * 	vector present.
 */
void find_cloning_vector(GapIO *io, int contig, int *left, int *right,
			 int svec_also, int min_length);

/*
 * Clip suggested start and end points for new sequences based on the known
 * vector found in an existing sequence.
 *
 * Arguments:
 *	io		Gap IO handle
 *	s_start		Current suggestion for new sequence start in contig
 *	s_end		Current suggestion for new sequence end in contig
 *	rnum		Reading number to compare against
 *	min_length	Minimum length of SVEC tag to consider it to be valid
 *
 * Returns:
 *	Void return, but updates s_start and s_end.
 */
void finish_clip_svec(GapIO *io, int *s_start, int *s_end, int rnum,
		      int min_length);

/*
 * Finds the readings covering a specific consensus base
 *
 * Arguments:
 *	io		Gap IO handle
 *	contig		Contig number
 *	pos		Contig position
 *
 * Returns:
 *	A malloc array of integer reading numbers. It is up to the caller
 *	to free this array with xfree.
 *	NULL for failure
 */
int *seqs_at_pos(GapIO *io, int contig, int pos);

/*
 * From a template (zero for none), finds a reading appropriate for
 * tagging (positions start to end).
 *
 * Ideally we choose a reading using that template.
 * If none exists then we choose any reading.
 * If this doesn't work then the tag gets truncated.
 *
 * Arguments:
 *	io		Gap IO handle
 *	contig		Contig number
 *	template	Template number (0 for none)
 *	start		First base in tag (offset into contig)
 *	end		Last base in tag
 *
 * Returns:
 *	Success: Template number. Also may modify start and end values.
 *	Failure: 0
 */
int tag_template(GapIO *io, int contig, int template,
		 int *start, int *end);

/*
 * Finds where a template has a duplicate or not.
 *
 * Returns 1 for yes,
 *         0 for no.
 */
int template_is_dup(finish_t *fin, int *templates_picked,
		    int num_picked,
		    int template);


/*
 * secondary_primer_match
 *
 * Identifies whether there is a high scoring match with a particular primer
 * elsewhere in a set of files, within the entire database sequence, or
 * within just a section of sequence.
 *
 * This uses fin->extern_seq buffer when check_external is set to true.
 * If check_contig is zero then do not compare against the consensus.
 * If check_contig >0 then compare against region of that contig (between
 * contig_start and contig_end).
 * If check_contig <0 then compare against all contigs in the database.
 *
 * If skip_self is set then the self_match (finding your own match in
 * the consensus) is ignored.
 *
 * Arguments:
 *	fin			'Finish' object
 *	check_contig		Contig to check against (0 => none, -1 => all)
 *	check_start		Start position in check_contig
 *	check_end		End position in check_contig
 *	self_match		Whether to ignore the self-match (true if so)
 *	self_strand		Which strand is self_match relevant to
 *	check_external		Should we check the external file list?
 *	prim			The primer to find matches to.
 *
 * Returns:
 * 	0	for no match (other than self)
 *	>=0	for match
 */
double secondary_primer_match(finish_t *fin,
			      int check_contig,
			      int check_start,
			      int check_end,
			      int self_match,
			      int self_strand,
			      int check_external,
			      char *prim);

#endif /* _FINISH_UTILS_H */
