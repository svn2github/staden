SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

distsrc: distsrc_dirs
	cp $(S)/*.tcl $(S)/trevrc $(S)/trev.in $(S)/trev.bat $(S)/Makefile $(DIRNAME)

install:
	$(INSTALL) trev $(INSTALLBIN)
	-mkdir $(INSTALLTCL)/trev
	cp $(S)/trev.tcl $(INSTALLTCL)/trev
	cp $(S)/trevrc $(INSTALLETC)
