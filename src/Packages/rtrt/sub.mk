#Makefile fragment for the Packages/rtrt directory

include $(SCIRUN_SCRIPTS)/largeso_prologue.mk

SRCDIR := Packages/rtrt

ifeq ($(findstring -n32, $(C_FLAGS)),-n32)
#ifneq ($(USE_SOUND),no)
    SOUNDDIR := $(SRCDIR)/Sound
endif

SUBDIRS := \
	$(SRCDIR)/Core         \
	$(SRCDIR)/Dataflow     \
	$(SOUNDDIR)            \
	$(SRCDIR)/visinfo

include $(SCIRUN_SCRIPTS)/recurse.mk

PSELIBS := Core
LIBS := $(OOGL_LIBRARY) $(GL_LIBRARY) $(FASTM_LIBRARY) $(M_LIBRARY) $(THREAD_LIBRARY) $(PERFEX_LIBRARY) 

include $(SCIRUN_SCRIPTS)/largeso_epilogue.mk

SUBDIRS := \
	$(SRCDIR)/StandAlone
include $(SCIRUN_SCRIPTS)/recurse.mk

#CFLAGS := $(CFLAGS) -OPT:IEEE_arithmetic=3 
#CFLAGS := $(CFLAGS) -OPT:Olimit=16383
#CFLAGS := $(CFLAGS) -OPT:Olimit_opt=on
