# *** NOTE ***
#
# Do not remove or modify the comment line:
#
# #[INSERT NEW ?????? HERE]
#
# It is required by the Component Wizard to properly edit this file.
# if you want to edit this file by hand, see the "Create A New Component"
# documentation on how to do it correctly.

SRCDIR := Packages/VS/Dataflow/GUI

SRCDIR := Packages/VS/Dataflow/GUI

SRCS := \
	$(SRCDIR)/HotBox.tcl          \
	$(SRCDIR)/ICUMonitor.tcl \
	$(SRCDIR)/ExecutiveState.tcl

include $(SCIRUN_SCRIPTS)/tclIndex.mk


