#
# Makefile for polyA_clip
#

PROGS= 	polyA_clip

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(IOLIB_INC) $(SEQUTILS_INC) $(TKUTILS_INC) $(MISC_INC)

OBJ=\
	polyA_clip.o\
	seqInfo.o

POLYA_CLIP_LIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB)

DEPS=\
	$(IOLIB_DEP) \
	$(MISC_DEP)


polyA_clip: $(OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJ) $(POLYA_CLIP_LIBS) $(LIBSC) $(DEPS)

DEPEND_OBJ = $(OBJ)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE THIS LINE -- make depend depends on it

polyA_clip.o: $(SRCROOT)/Misc/xalloc.h
polyA_clip.o: $(SRCROOT)/seq_utils/dna_utils.h
polyA_clip.o: seqInfo.h
seqInfo.o: seqInfo.h
