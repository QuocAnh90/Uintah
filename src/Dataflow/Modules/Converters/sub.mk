# *** NOTE ***
#
# Do not remove or modify the comment line:
#
# #[INSERT NEW ?????? HERE]
#
# It is required by the Component Wizard to properly edit this file.
# if you want to edit this file by hand, see the "Create A New Component"
# documentation on how to do it correctly.

include $(SCIRUN_SCRIPTS)/smallso_prologue.mk

SRCDIR   := Dataflow/Modules/Converters

SRCS     += \
	$(SRCDIR)/ConvertMatrixToString.cc\
	$(SRCDIR)/ConvertMatrixToField.cc\
	$(SRCDIR)/ConvertNrrdToField.cc\
#[INSERT NEW CODE FILE HERE]

PSELIBS := Core/Datatypes Dataflow/Network \
        Core/Persistent Core/Containers Core/Util \
        Core/Exceptions Core/Thread Dataflow/GuiInterface \
        Core/Geom Core/Datatypes Core/Geometry \
        Core/GeomInterface Dataflow/TkExtensions \
        Core/Algorithms/Converter \


LIBS := $(TK_LIBRARY) $(GL_LIBRARY) $(M_LIBRARY)

include $(SCIRUN_SCRIPTS)/smallso_epilogue.mk


