#
# Makefile for find_renz
#

PROGS = find_renz

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(SEQUTILS_INC) $(MISC_INC)

OBJS	= find_renz.o

LIBS=\
	$(SEQUTILS_LIB) \
	$(TEXTUTILS_LIB) \
	$(MISC_LIB)

find_renz: $(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(LIBS) $(LIBSC) 


DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE THIS LINE -- make depend depends on it.

find_renz.o: $(SRCROOT)/Misc/getfile.h
find_renz.o: $(SRCROOT)/Misc/misc.h
find_renz.o: $(SRCROOT)/Misc/os.h
find_renz.o: $(SRCROOT)/Misc/xalloc.h
find_renz.o: $(SRCROOT)/seq_utils/dna_utils.h
find_renz.o: $(SRCROOT)/seq_utils/genetic_code.h
find_renz.o: $(SRCROOT)/seq_utils/renz_utils.h
find_renz.o: $(SRCROOT)/seq_utils/sequence_formats.h
find_renz.o: $(PWD)/staden_config.h
