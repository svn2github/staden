#
# Makefile for polyA_clip
#

PROGS= 	$(O)/polyA_clip

SRCROOT=..
#SRCROOT=/nfs/arran/home5/pubseq/share/src

include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC) $(SEQUTILS_INC) $(TKUTILS_INC)

POLYA_CLIPSRC=$(O)

all: $(O)/polyA_clip

OBJ=\
	$(POLYA_CLIPSRC)/polyA_clip.o\
	$(POLYA_CLIPSRC)/seqInfo.o

POLYA_CLIP_LIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB)

DEPS=\
	$(IOLIB_DEP) \
	$(MISC_DEP)


$(O)/polyA_clip: $(OBJ)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(OBJ) $(POLYA_CLIP_LIBS) $(LIBSC) $(DEPS)

DEPEND_OBJ = $(OBJ)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
