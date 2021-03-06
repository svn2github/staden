#
# Makefile for tk_utils routines
#

LIBS = tk_utils
ifneq ($(MACHINE),windows)
PROGS= $(LIBS) $(O)/stash
else
PROGS= $(LIBS)
endif
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(TKUTILS_INC) $(IOLIB_INC) $(TK_INC) $(SEQUTILS_INC) \
	$(BIOLIMS_INC) $(PNG_INC)

CFLAGS += $(SHLIB_CFLAGS) $(TK_UTILS_DLL)

TKUTILSBIN=$(O)

#DEFINES += -DCHECK_LICENCE -DTCL_LOCAL_APPINIT=Stash_AppInit
DEFINES += -DTCL_LOCAL_APPINIT=Stash_AppInit

#
# Objects
#
OBJS = \
	$(TKUTILSBIN)/cli_arg.o \
	$(TKUTILSBIN)/tclXkeylist.o \
	$(TKUTILSBIN)/tclXutil.o \
	$(TKUTILSBIN)/tcl_utils.o \
	$(TKUTILSBIN)/tcl_debug.o \
	$(TKUTILSBIN)/misc.o \
	$(TKUTILSBIN)/init.o \
	$(TKUTILSBIN)/text_output.o \
	$(TKUTILSBIN)/tkRaster.o \
	$(TKUTILSBIN)/tkRasterBuiltIn.o \
	$(TKUTILSBIN)/sheet.o \
	$(TKUTILSBIN)/tkSheet.o \
	$(TKUTILSBIN)/tkSheet_common.o \
	$(TKUTILSBIN)/trace_print.o \
	$(TKUTILSBIN)/postscript.o \
	$(TKUTILSBIN)/split.o \
	$(TKUTILSBIN)/tkTrace.o \
	$(TKUTILSBIN)/tkTraceComp.o \
	$(TKUTILSBIN)/tkTraceIO.o \
	$(TKUTILSBIN)/tkTraceDisp.o \
	$(TKUTILSBIN)/capture.o \
	$(TKUTILSBIN)/canvas_box.o \
	$(TKUTILSBIN)/ruler_tick.o \
	$(TKUTILSBIN)/restriction_enzyme_map.o \
	$(TKUTILSBIN)/check_licence.o\
	$(TKUTILSBIN)/licence_utils.o\
	$(TKUTILSBIN)/valid_seq.o\
	$(TKUTILSBIN)/md52.o\
	$(TKUTILSBIN)/container_cmds.o\
	$(TKUTILSBIN)/container.o\
	$(TKUTILSBIN)/container_ruler.o\
	$(TKUTILSBIN)/element_canvas.o\
	$(TKUTILSBIN)/tclCanvGraph.o\
	$(TKUTILSBIN)/tkCanvGraph.o\
	$(TKUTILSBIN)/matrix.o\
	$(TKUTILS_EXTRAS)

#
# Main dependency
#
$(LIBS) : .links $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

TKU_LIBS = \
	$(IOLIB_LIB) \
	$(TK_LIB) \
	$(TCL_LIB) \
	$(MISC_LIB) \
	$(PNG_LIB) \
	$(ZLIB_LIB) \
	$(SOCKET)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(TKU_LIBS)


$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

# Stash interpreter executable
# Uses the main tcl interpreter, but overrides the appinit to use a local
# function.
ifneq ($(MACHINE),windows)
STASH_OBJS = \
	$(TKUTILSBIN)/stash.o \
	$(TKUTILSBIN)/tclAppInit.o

STASH_LIBS = \
	$(TKUTILS_LIB) \
	$(TK_LIB)

$(TKUTILSBIN)/stash: $(STASH_OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMWIN) $(STASH_OBJS) $(STASH_LIBS) $(LIBSC)
endif

.PRECIOUS:	.links
.links:
	ln -s ../licence/boxes.h .
	ln -s ../licence/licence.h .
	ln -s ../licence/check_licence.c .
	ln -s ../licence/licence_utils.c .
	ln -s ../licence/licence_utils.h .
	ln -s ../licence/md52.c .
	ln -s ../licence/md52.h .
	ln -s ../licence/valid_seq.c .
	ln -s ../tcl8.4.6/unix/tclAppInit.c .
	ln -s ../tk8.4.6/winMain.c .
	touch .links
	$(MAKE) $(MFLAGS)

# Copy the tk_utils source
distsrc: distsrc_dirs
	-cp -R *.[ch] *.tcl tclIndex Makefile dependencies \
		$(DIRNAME)
	-rm $(DIRNAME)/resource.h

install:
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)
	-mkdir $(INSTALLLIB)/tk_utils
	cp *.tcl tclIndex $(INSTALLLIB)/tk_utils
ifneq ($(MACHINE),windows)
	cp $(O)/stash $(INSTALLBIN)
endif

include dependencies
# DO NOT DELETE
