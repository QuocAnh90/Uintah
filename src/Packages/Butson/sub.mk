#Makefile fragment for the Packages/Butson directory

include $(OBJTOP_ABS)/scripts/largeso_prologue.mk

SRCDIR := Packages/Butson
SUBDIRS := \
	$(SRCDIR)/Dataflow \
#	$(SRCDIR)/Core \

include $(OBJTOP_ABS)/scripts/recurse.mk

PSELIBS := 
LIBS := $(TK_LIBRARY) $(GL_LIBS) -lm

include $(OBJTOP_ABS)/scripts/largeso_epilogue.mk
