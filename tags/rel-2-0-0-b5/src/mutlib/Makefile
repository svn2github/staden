#
# Makefile for Mutation Library
#
LIBS = mutlib
PROGS= $(LIBS)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

CXXFLAGS += $(SHLIB_CXXFLAGS) $(MISC_DLL)

INCLUDES_E += $(MISC_INC) $(MUT_INC) $(SEQUTILS_INC) $(IOLIB_INC)

DEFINES += -DMUTLIB_EXPORTS -DDYNMAT

#
# Objects
#
OBJS = \
	align.o \
	basecall.o \
	caller.o \
	caller_snr.o \
	caller_base.o \
	caller_level.o \
	dnatable.o \
	muttag.o \
	mutationtag.o \
	mutationtag_utils.o \
	pathutil.o \
	peakcall.o \
	stringlist.o \
	tagarray.o \
	trace.o \
	validate.o \
	read_matrix.o \
	sp_align_lib.o \
	sp_alignment.o \
	sp_hash_lib.o \
	sp_matrix.o \
	tracediff.o \
	tracediff_helper.o \
	tracediff_mutations.o \
	tracediff_parameters.o \
	tracediff_validate.o \
	tracealign.o \
	tracealign_cache.o \
	tracealign_helper.o \
	tracealign_interpolate.o \
	tracealign_preprocess.o \
	tracealign_validate.o \
	mutscan.o \
	mutscan_helper.o \
	mutscan_validate.o \
	mutscan_analysis.o \
	mutscan_preprocess.o \
	mutscan_parameters.o


MUT_DEP=\
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB) \
	$(MATH_LIB)

ifeq ($(MACHINE),macosx)
MUT_DEP += $(TK_LIB)
endif

#ifeq ($(MACHINE),alpha)
#MUT_DEP += -Bstatic -lgcc
#endif

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS)
	-cp $(MUTSRC)/$(LIBS).def $(DEF_FILE) 2>&-
	$(SHLIB_LDXX) $(SHLIB_LDXXFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(MUT_DEP)

DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.[ch]pp $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGLIBS) $(INSTALLLIB)

# DO NOT DELETE THIS LINE -- make depend depends on it

align.o: $(SRCROOT)/mutlib/align.hpp
align.o: $(SRCROOT)/mutlib/matrix.hpp
align.o: $(SRCROOT)/mutlib/sp_alignment.h
align.o: $(SRCROOT)/mutlib/sp_alignment_structs.h
basecall.o: $(SRCROOT)/mutlib/basecall.hpp
caller.o: $(SRCROOT)/Misc/misc.h
caller.o: $(SRCROOT)/Misc/os.h
caller.o: $(SRCROOT)/Misc/xalloc.h
caller.o: $(SRCROOT)/mutlib/caller.hpp
caller.o: $(SRCROOT)/mutlib/matrix.hpp
caller.o: $(SRCROOT)/mutlib/staden.h
caller.o: $(SRCROOT)/mutlib/trace.hpp
caller.o: $(SRCROOT)/seq_utils/dna_utils.h
caller.o: $(PWD)/staden_config.h
caller_base.o: $(SRCROOT)/Misc/misc.h
caller_base.o: $(SRCROOT)/Misc/os.h
caller_base.o: $(SRCROOT)/Misc/xalloc.h
caller_base.o: $(SRCROOT)/mutlib/caller.hpp
caller_base.o: $(SRCROOT)/mutlib/caller_base.hpp
caller_base.o: $(SRCROOT)/mutlib/dnatable.hpp
caller_base.o: $(SRCROOT)/mutlib/matrix.hpp
caller_base.o: $(SRCROOT)/mutlib/staden.h
caller_base.o: $(SRCROOT)/mutlib/trace.hpp
caller_base.o: $(SRCROOT)/seq_utils/dna_utils.h
caller_base.o: $(PWD)/staden_config.h
caller_level.o: $(SRCROOT)/Misc/misc.h
caller_level.o: $(SRCROOT)/Misc/os.h
caller_level.o: $(SRCROOT)/Misc/xalloc.h
caller_level.o: $(SRCROOT)/mutlib/caller.hpp
caller_level.o: $(SRCROOT)/mutlib/caller_level.hpp
caller_level.o: $(SRCROOT)/mutlib/dnatable.hpp
caller_level.o: $(SRCROOT)/mutlib/matrix.hpp
caller_level.o: $(SRCROOT)/mutlib/staden.h
caller_level.o: $(SRCROOT)/mutlib/trace.hpp
caller_level.o: $(SRCROOT)/seq_utils/dna_utils.h
caller_level.o: $(PWD)/staden_config.h
caller_snr.o: $(SRCROOT)/Misc/misc.h
caller_snr.o: $(SRCROOT)/Misc/os.h
caller_snr.o: $(SRCROOT)/Misc/xalloc.h
caller_snr.o: $(SRCROOT)/mutlib/caller.hpp
caller_snr.o: $(SRCROOT)/mutlib/caller_snr.hpp
caller_snr.o: $(SRCROOT)/mutlib/matrix.hpp
caller_snr.o: $(SRCROOT)/mutlib/staden.h
caller_snr.o: $(SRCROOT)/mutlib/trace.hpp
caller_snr.o: $(SRCROOT)/seq_utils/dna_utils.h
caller_snr.o: $(PWD)/staden_config.h
dnatable.o: $(SRCROOT)/mutlib/dnatable.hpp
mutationtag.o: $(SRCROOT)/Misc/misc.h
mutationtag.o: $(SRCROOT)/Misc/os.h
mutationtag.o: $(SRCROOT)/Misc/xalloc.h
mutationtag.o: $(SRCROOT)/mutlib/listitem.hpp
mutationtag.o: $(SRCROOT)/mutlib/mutationtag.hpp
mutationtag.o: $(SRCROOT)/mutlib/mutlib.h
mutationtag.o: $(SRCROOT)/mutlib/staden.h
mutationtag.o: $(SRCROOT)/seq_utils/dna_utils.h
mutationtag.o: $(PWD)/staden_config.h
mutationtag_utils.o: $(SRCROOT)/Misc/misc.h
mutationtag_utils.o: $(SRCROOT)/Misc/os.h
mutationtag_utils.o: $(SRCROOT)/Misc/xalloc.h
mutationtag_utils.o: $(SRCROOT)/mutlib/array.hpp
mutationtag_utils.o: $(SRCROOT)/mutlib/list.hpp
mutationtag_utils.o: $(SRCROOT)/mutlib/listitem.hpp
mutationtag_utils.o: $(SRCROOT)/mutlib/mutationtag.hpp
mutationtag_utils.o: $(SRCROOT)/mutlib/mutationtag_utils.hpp
mutationtag_utils.o: $(SRCROOT)/mutlib/mutlib.h
mutationtag_utils.o: $(SRCROOT)/mutlib/staden.h
mutationtag_utils.o: $(SRCROOT)/seq_utils/dna_utils.h
mutationtag_utils.o: $(PWD)/staden_config.h
mutscan.o: $(SRCROOT)/Misc/misc.h
mutscan.o: $(SRCROOT)/Misc/os.h
mutscan.o: $(SRCROOT)/Misc/xalloc.h
mutscan.o: $(SRCROOT)/mutlib/align.hpp
mutscan.o: $(SRCROOT)/mutlib/array.hpp
mutscan.o: $(SRCROOT)/mutlib/list.hpp
mutscan.o: $(SRCROOT)/mutlib/listitem.hpp
mutscan.o: $(SRCROOT)/mutlib/matrix.hpp
mutscan.o: $(SRCROOT)/mutlib/mutationtag.hpp
mutscan.o: $(SRCROOT)/mutlib/mutationtag_utils.hpp
mutscan.o: $(SRCROOT)/mutlib/mutlib.h
mutscan.o: $(SRCROOT)/mutlib/mutscan.hpp
mutscan.o: $(SRCROOT)/mutlib/mutscan_analysis.hpp
mutscan.o: $(SRCROOT)/mutlib/mutscan_parameters.hpp
mutscan.o: $(SRCROOT)/mutlib/mutscan_preprocess.hpp
mutscan.o: $(SRCROOT)/mutlib/parameter.hpp
mutscan.o: $(SRCROOT)/mutlib/sp_alignment.h
mutscan.o: $(SRCROOT)/mutlib/sp_alignment_structs.h
mutscan.o: $(SRCROOT)/mutlib/staden.h
mutscan.o: $(SRCROOT)/mutlib/trace.hpp
mutscan.o: $(SRCROOT)/seq_utils/dna_utils.h
mutscan.o: $(PWD)/staden_config.h
mutscan_analysis.o: $(SRCROOT)/Misc/misc.h
mutscan_analysis.o: $(SRCROOT)/Misc/os.h
mutscan_analysis.o: $(SRCROOT)/Misc/xalloc.h
mutscan_analysis.o: $(SRCROOT)/mutlib/array.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/caller.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/caller_base.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/caller_level.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/dnatable.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/list.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/listitem.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/matrix.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/mutationtag.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/mutlib.h
mutscan_analysis.o: $(SRCROOT)/mutlib/mutscan_analysis.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/mutscan_preprocess.hpp
mutscan_analysis.o: $(SRCROOT)/mutlib/staden.h
mutscan_analysis.o: $(SRCROOT)/mutlib/trace.hpp
mutscan_analysis.o: $(SRCROOT)/seq_utils/dna_utils.h
mutscan_analysis.o: $(PWD)/staden_config.h
mutscan_helper.o: $(SRCROOT)/Misc/misc.h
mutscan_helper.o: $(SRCROOT)/Misc/os.h
mutscan_helper.o: $(SRCROOT)/Misc/xalloc.h
mutscan_helper.o: $(SRCROOT)/mutlib/mutlib.h
mutscan_helper.o: $(SRCROOT)/mutlib/staden.h
mutscan_helper.o: $(SRCROOT)/seq_utils/dna_utils.h
mutscan_helper.o: $(PWD)/staden_config.h
mutscan_parameters.o: $(SRCROOT)/Misc/misc.h
mutscan_parameters.o: $(SRCROOT)/Misc/os.h
mutscan_parameters.o: $(SRCROOT)/Misc/xalloc.h
mutscan_parameters.o: $(SRCROOT)/mutlib/mutlib.h
mutscan_parameters.o: $(SRCROOT)/mutlib/mutscan_parameters.hpp
mutscan_parameters.o: $(SRCROOT)/mutlib/parameter.hpp
mutscan_parameters.o: $(SRCROOT)/mutlib/staden.h
mutscan_parameters.o: $(SRCROOT)/seq_utils/dna_utils.h
mutscan_parameters.o: $(PWD)/staden_config.h
mutscan_preprocess.o: $(SRCROOT)/Misc/misc.h
mutscan_preprocess.o: $(SRCROOT)/Misc/os.h
mutscan_preprocess.o: $(SRCROOT)/Misc/xalloc.h
mutscan_preprocess.o: $(SRCROOT)/mutlib/array.hpp
mutscan_preprocess.o: $(SRCROOT)/mutlib/matrix.hpp
mutscan_preprocess.o: $(SRCROOT)/mutlib/mutlib.h
mutscan_preprocess.o: $(SRCROOT)/mutlib/mutscan_preprocess.hpp
mutscan_preprocess.o: $(SRCROOT)/mutlib/staden.h
mutscan_preprocess.o: $(SRCROOT)/mutlib/trace.hpp
mutscan_preprocess.o: $(SRCROOT)/seq_utils/dna_utils.h
mutscan_preprocess.o: $(PWD)/staden_config.h
mutscan_validate.o: $(SRCROOT)/Misc/misc.h
mutscan_validate.o: $(SRCROOT)/Misc/os.h
mutscan_validate.o: $(SRCROOT)/Misc/xalloc.h
mutscan_validate.o: $(SRCROOT)/mutlib/mutlib.h
mutscan_validate.o: $(SRCROOT)/mutlib/mutscan_parameters.hpp
mutscan_validate.o: $(SRCROOT)/mutlib/parameter.hpp
mutscan_validate.o: $(SRCROOT)/mutlib/staden.h
mutscan_validate.o: $(SRCROOT)/mutlib/validate.hpp
mutscan_validate.o: $(SRCROOT)/seq_utils/dna_utils.h
mutscan_validate.o: $(PWD)/staden_config.h
muttag.o: $(SRCROOT)/Misc/misc.h
muttag.o: $(SRCROOT)/Misc/os.h
muttag.o: $(SRCROOT)/Misc/xalloc.h
muttag.o: $(SRCROOT)/mutlib/listitem.hpp
muttag.o: $(SRCROOT)/mutlib/mutlib.h
muttag.o: $(SRCROOT)/mutlib/muttag.hpp
muttag.o: $(SRCROOT)/mutlib/staden.h
muttag.o: $(SRCROOT)/seq_utils/dna_utils.h
muttag.o: $(PWD)/staden_config.h
peakcall.o: $(SRCROOT)/mutlib/listitem.hpp
peakcall.o: $(SRCROOT)/mutlib/peakcall.hpp
read_matrix.o: $(SRCROOT)/Misc/misc.h
read_matrix.o: $(SRCROOT)/Misc/os.h
read_matrix.o: $(SRCROOT)/Misc/xalloc.h
read_matrix.o: $(PWD)/staden_config.h
read_matrix.o: read_matrix.h
sp_align_lib.o: $(SRCROOT)/Misc/misc.h
sp_align_lib.o: $(SRCROOT)/Misc/os.h
sp_align_lib.o: $(SRCROOT)/Misc/xalloc.h
sp_align_lib.o: $(SRCROOT)/mutlib/sp_align_lib.h
sp_align_lib.o: $(SRCROOT)/mutlib/sp_alignment.h
sp_align_lib.o: $(SRCROOT)/mutlib/sp_alignment_structs.h
sp_align_lib.o: $(SRCROOT)/mutlib/staden.h
sp_align_lib.o: $(SRCROOT)/seq_utils/dna_utils.h
sp_align_lib.o: $(PWD)/staden_config.h
sp_align_lib.o: read_matrix.h
sp_alignment.o: $(SRCROOT)/Misc/misc.h
sp_alignment.o: $(SRCROOT)/Misc/os.h
sp_alignment.o: $(SRCROOT)/Misc/xalloc.h
sp_alignment.o: $(SRCROOT)/mutlib/sp_align_lib.h
sp_alignment.o: $(SRCROOT)/mutlib/sp_alignment.h
sp_alignment.o: $(SRCROOT)/mutlib/sp_alignment_structs.h
sp_alignment.o: $(SRCROOT)/mutlib/sp_hash_lib.h
sp_alignment.o: $(SRCROOT)/mutlib/staden.h
sp_alignment.o: $(SRCROOT)/seq_utils/dna_utils.h
sp_alignment.o: $(PWD)/staden_config.h
sp_alignment.o: read_matrix.h
sp_hash_lib.o: $(SRCROOT)/Misc/misc.h
sp_hash_lib.o: $(SRCROOT)/Misc/os.h
sp_hash_lib.o: $(SRCROOT)/Misc/xalloc.h
sp_hash_lib.o: $(SRCROOT)/mutlib/sp_align_lib.h
sp_hash_lib.o: $(SRCROOT)/mutlib/sp_alignment.h
sp_hash_lib.o: $(SRCROOT)/mutlib/sp_alignment_structs.h
sp_hash_lib.o: $(SRCROOT)/mutlib/sp_hash_lib.h
sp_hash_lib.o: $(SRCROOT)/mutlib/staden.h
sp_hash_lib.o: $(SRCROOT)/seq_utils/dna_utils.h
sp_hash_lib.o: $(PWD)/staden_config.h
sp_matrix.o: $(SRCROOT)/Misc/misc.h
sp_matrix.o: $(SRCROOT)/Misc/os.h
sp_matrix.o: $(SRCROOT)/Misc/xalloc.h
sp_matrix.o: $(SRCROOT)/mutlib/sp_matrix.h
sp_matrix.o: $(PWD)/staden_config.h
stringlist.o: $(SRCROOT)/mutlib/stringlist.hpp
tagarray.o: $(SRCROOT)/Misc/misc.h
tagarray.o: $(SRCROOT)/Misc/os.h
tagarray.o: $(SRCROOT)/Misc/xalloc.h
tagarray.o: $(SRCROOT)/mutlib/list.hpp
tagarray.o: $(SRCROOT)/mutlib/listitem.hpp
tagarray.o: $(SRCROOT)/mutlib/mutlib.h
tagarray.o: $(SRCROOT)/mutlib/muttag.hpp
tagarray.o: $(SRCROOT)/mutlib/staden.h
tagarray.o: $(SRCROOT)/mutlib/tagarray.hpp
tagarray.o: $(SRCROOT)/seq_utils/dna_utils.h
tagarray.o: $(PWD)/staden_config.h
trace.o: $(SRCROOT)/Misc/misc.h
trace.o: $(SRCROOT)/Misc/os.h
trace.o: $(SRCROOT)/Misc/xalloc.h
trace.o: $(SRCROOT)/mutlib/array.hpp
trace.o: $(SRCROOT)/mutlib/staden.h
trace.o: $(SRCROOT)/mutlib/trace.hpp
trace.o: $(SRCROOT)/seq_utils/dna_utils.h
trace.o: $(PWD)/staden_config.h
tracealign.o: $(SRCROOT)/Misc/misc.h
tracealign.o: $(SRCROOT)/Misc/os.h
tracealign.o: $(SRCROOT)/Misc/xalloc.h
tracealign.o: $(SRCROOT)/mutlib/align.hpp
tracealign.o: $(SRCROOT)/mutlib/array.hpp
tracealign.o: $(SRCROOT)/mutlib/list.hpp
tracealign.o: $(SRCROOT)/mutlib/listitem.hpp
tracealign.o: $(SRCROOT)/mutlib/matrix.hpp
tracealign.o: $(SRCROOT)/mutlib/mutlib.h
tracealign.o: $(SRCROOT)/mutlib/muttag.hpp
tracealign.o: $(SRCROOT)/mutlib/peakcall.hpp
tracealign.o: $(SRCROOT)/mutlib/sp_alignment.h
tracealign.o: $(SRCROOT)/mutlib/sp_alignment_structs.h
tracealign.o: $(SRCROOT)/mutlib/staden.h
tracealign.o: $(SRCROOT)/mutlib/tagarray.hpp
tracealign.o: $(SRCROOT)/mutlib/trace.hpp
tracealign.o: $(SRCROOT)/mutlib/tracealign.hpp
tracealign.o: $(SRCROOT)/mutlib/tracealign_cache.hpp
tracealign.o: $(SRCROOT)/mutlib/tracealign_preprocess.hpp
tracealign.o: $(SRCROOT)/seq_utils/dna_utils.h
tracealign.o: $(PWD)/staden_config.h
tracealign_cache.o: $(SRCROOT)/Misc/misc.h
tracealign_cache.o: $(SRCROOT)/Misc/os.h
tracealign_cache.o: $(SRCROOT)/Misc/xalloc.h
tracealign_cache.o: $(SRCROOT)/mutlib/array.hpp
tracealign_cache.o: $(SRCROOT)/mutlib/matrix.hpp
tracealign_cache.o: $(SRCROOT)/mutlib/staden.h
tracealign_cache.o: $(SRCROOT)/mutlib/trace.hpp
tracealign_cache.o: $(SRCROOT)/mutlib/tracealign_cache.hpp
tracealign_cache.o: $(SRCROOT)/mutlib/tracealign_preprocess.hpp
tracealign_cache.o: $(SRCROOT)/seq_utils/dna_utils.h
tracealign_cache.o: $(PWD)/staden_config.h
tracealign_helper.o: $(SRCROOT)/Misc/misc.h
tracealign_helper.o: $(SRCROOT)/Misc/os.h
tracealign_helper.o: $(SRCROOT)/Misc/xalloc.h
tracealign_helper.o: $(SRCROOT)/mutlib/array.hpp
tracealign_helper.o: $(SRCROOT)/mutlib/matrix.hpp
tracealign_helper.o: $(SRCROOT)/mutlib/mutlib.h
tracealign_helper.o: $(SRCROOT)/mutlib/staden.h
tracealign_helper.o: $(SRCROOT)/mutlib/trace.hpp
tracealign_helper.o: $(SRCROOT)/mutlib/tracealign_cache.hpp
tracealign_helper.o: $(SRCROOT)/mutlib/tracealign_preprocess.hpp
tracealign_helper.o: $(SRCROOT)/seq_utils/dna_utils.h
tracealign_helper.o: $(PWD)/staden_config.h
tracealign_interpolate.o: $(SRCROOT)/Misc/misc.h
tracealign_interpolate.o: $(SRCROOT)/Misc/os.h
tracealign_interpolate.o: $(SRCROOT)/Misc/xalloc.h
tracealign_interpolate.o: $(SRCROOT)/mutlib/align.hpp
tracealign_interpolate.o: $(SRCROOT)/mutlib/array.hpp
tracealign_interpolate.o: $(SRCROOT)/mutlib/matrix.hpp
tracealign_interpolate.o: $(SRCROOT)/mutlib/sp_alignment.h
tracealign_interpolate.o: $(SRCROOT)/mutlib/sp_alignment_structs.h
tracealign_interpolate.o: $(SRCROOT)/mutlib/staden.h
tracealign_interpolate.o: $(SRCROOT)/mutlib/trace.hpp
tracealign_interpolate.o: $(SRCROOT)/seq_utils/dna_utils.h
tracealign_interpolate.o: $(PWD)/staden_config.h
tracealign_preprocess.o: $(SRCROOT)/Misc/misc.h
tracealign_preprocess.o: $(SRCROOT)/Misc/os.h
tracealign_preprocess.o: $(SRCROOT)/Misc/xalloc.h
tracealign_preprocess.o: $(SRCROOT)/mutlib/array.hpp
tracealign_preprocess.o: $(SRCROOT)/mutlib/basecall.hpp
tracealign_preprocess.o: $(SRCROOT)/mutlib/staden.h
tracealign_preprocess.o: $(SRCROOT)/mutlib/trace.hpp
tracealign_preprocess.o: $(SRCROOT)/mutlib/tracealign_preprocess.hpp
tracealign_preprocess.o: $(SRCROOT)/seq_utils/dna_utils.h
tracealign_preprocess.o: $(PWD)/staden_config.h
tracealign_validate.o: $(SRCROOT)/Misc/misc.h
tracealign_validate.o: $(SRCROOT)/Misc/os.h
tracealign_validate.o: $(SRCROOT)/Misc/xalloc.h
tracealign_validate.o: $(SRCROOT)/mutlib/mutlib.h
tracealign_validate.o: $(SRCROOT)/mutlib/staden.h
tracealign_validate.o: $(SRCROOT)/mutlib/validate.hpp
tracealign_validate.o: $(SRCROOT)/seq_utils/dna_utils.h
tracealign_validate.o: $(PWD)/staden_config.h
tracediff.o: $(SRCROOT)/Misc/misc.h
tracediff.o: $(SRCROOT)/Misc/os.h
tracediff.o: $(SRCROOT)/Misc/xalloc.h
tracediff.o: $(SRCROOT)/mutlib/align.hpp
tracediff.o: $(SRCROOT)/mutlib/array.hpp
tracediff.o: $(SRCROOT)/mutlib/list.hpp
tracediff.o: $(SRCROOT)/mutlib/listitem.hpp
tracediff.o: $(SRCROOT)/mutlib/matrix.hpp
tracediff.o: $(SRCROOT)/mutlib/mutlib.h
tracediff.o: $(SRCROOT)/mutlib/muttag.hpp
tracediff.o: $(SRCROOT)/mutlib/parameter.hpp
tracediff.o: $(SRCROOT)/mutlib/sp_alignment.h
tracediff.o: $(SRCROOT)/mutlib/sp_alignment_structs.h
tracediff.o: $(SRCROOT)/mutlib/staden.h
tracediff.o: $(SRCROOT)/mutlib/tagarray.hpp
tracediff.o: $(SRCROOT)/mutlib/trace.hpp
tracediff.o: $(SRCROOT)/mutlib/tracealign.hpp
tracediff.o: $(SRCROOT)/mutlib/tracediff.hpp
tracediff.o: $(SRCROOT)/mutlib/tracediff_parameters.hpp
tracediff.o: $(SRCROOT)/seq_utils/dna_utils.h
tracediff.o: $(PWD)/staden_config.h
tracediff_helper.o: $(SRCROOT)/Misc/misc.h
tracediff_helper.o: $(SRCROOT)/Misc/os.h
tracediff_helper.o: $(SRCROOT)/Misc/xalloc.h
tracediff_helper.o: $(SRCROOT)/mutlib/align.hpp
tracediff_helper.o: $(SRCROOT)/mutlib/array.hpp
tracediff_helper.o: $(SRCROOT)/mutlib/listitem.hpp
tracediff_helper.o: $(SRCROOT)/mutlib/matrix.hpp
tracediff_helper.o: $(SRCROOT)/mutlib/mutlib.h
tracediff_helper.o: $(SRCROOT)/mutlib/muttag.hpp
tracediff_helper.o: $(SRCROOT)/mutlib/sp_alignment.h
tracediff_helper.o: $(SRCROOT)/mutlib/sp_alignment_structs.h
tracediff_helper.o: $(SRCROOT)/mutlib/staden.h
tracediff_helper.o: $(SRCROOT)/mutlib/trace.hpp
tracediff_helper.o: $(SRCROOT)/mutlib/tracealign.hpp
tracediff_helper.o: $(SRCROOT)/seq_utils/dna_utils.h
tracediff_helper.o: $(PWD)/staden_config.h
tracediff_mutations.o: $(SRCROOT)/Misc/misc.h
tracediff_mutations.o: $(SRCROOT)/Misc/os.h
tracediff_mutations.o: $(SRCROOT)/Misc/xalloc.h
tracediff_mutations.o: $(SRCROOT)/mutlib/array.hpp
tracediff_mutations.o: $(SRCROOT)/mutlib/list.hpp
tracediff_mutations.o: $(SRCROOT)/mutlib/listitem.hpp
tracediff_mutations.o: $(SRCROOT)/mutlib/mutlib.h
tracediff_mutations.o: $(SRCROOT)/mutlib/muttag.hpp
tracediff_mutations.o: $(SRCROOT)/mutlib/parameter.hpp
tracediff_mutations.o: $(SRCROOT)/mutlib/peakcall.hpp
tracediff_mutations.o: $(SRCROOT)/mutlib/staden.h
tracediff_mutations.o: $(SRCROOT)/mutlib/trace.hpp
tracediff_mutations.o: $(SRCROOT)/mutlib/tracediff_parameters.hpp
tracediff_mutations.o: $(SRCROOT)/seq_utils/dna_utils.h
tracediff_mutations.o: $(PWD)/staden_config.h
tracediff_parameters.o: $(SRCROOT)/Misc/misc.h
tracediff_parameters.o: $(SRCROOT)/Misc/os.h
tracediff_parameters.o: $(SRCROOT)/Misc/xalloc.h
tracediff_parameters.o: $(SRCROOT)/mutlib/mutlib.h
tracediff_parameters.o: $(SRCROOT)/mutlib/parameter.hpp
tracediff_parameters.o: $(SRCROOT)/mutlib/staden.h
tracediff_parameters.o: $(SRCROOT)/mutlib/tracediff_parameters.hpp
tracediff_parameters.o: $(SRCROOT)/seq_utils/dna_utils.h
tracediff_parameters.o: $(PWD)/staden_config.h
tracediff_validate.o: $(SRCROOT)/Misc/misc.h
tracediff_validate.o: $(SRCROOT)/Misc/os.h
tracediff_validate.o: $(SRCROOT)/Misc/xalloc.h
tracediff_validate.o: $(SRCROOT)/mutlib/mutlib.h
tracediff_validate.o: $(SRCROOT)/mutlib/parameter.hpp
tracediff_validate.o: $(SRCROOT)/mutlib/staden.h
tracediff_validate.o: $(SRCROOT)/mutlib/tracediff_parameters.hpp
tracediff_validate.o: $(SRCROOT)/mutlib/validate.hpp
tracediff_validate.o: $(SRCROOT)/seq_utils/dna_utils.h
tracediff_validate.o: $(PWD)/staden_config.h
validate.o: $(SRCROOT)/Misc/misc.h
validate.o: $(SRCROOT)/Misc/os.h
validate.o: $(SRCROOT)/Misc/xalloc.h
validate.o: $(SRCROOT)/mutlib/mutlib.h
validate.o: $(SRCROOT)/mutlib/staden.h
validate.o: $(SRCROOT)/seq_utils/dna_utils.h
validate.o: $(PWD)/staden_config.h
