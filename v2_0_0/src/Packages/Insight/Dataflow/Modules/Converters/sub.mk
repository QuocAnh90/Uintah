#
#  The contents of this file are subject to the University of Utah Public
#  License (the "License"); you may not use this file except in compliance
#  with the License.
#  
#  Software distributed under the License is distributed on an "AS IS"
#  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#  License for the specific language governing rights and limitations under
#  the License.
#  
#  The Original Source Code is SCIRun, released March 12, 2001.
#  
#  The Original Source Code was developed by the University of Utah.
#  Portions created by UNIVERSITY are Copyright (C) 2001, 1994 
#  University of Utah. All Rights Reserved.
#


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

SRCDIR   := Packages/Insight/Dataflow/Modules/Converters

SRCS     += \
	$(SRCDIR)/FieldToImage.cc\
	$(SRCDIR)/ImageToField.cc\
	$(SRCDIR)/UCharToDouble.cc\
	$(SRCDIR)/UCharToFloat.cc\
	$(SRCDIR)/CharToUChar.cc\
	$(SRCDIR)/CharToDouble.cc\
	$(SRCDIR)/CharToFloat.cc\
	$(SRCDIR)/UShortToUChar.cc\
	$(SRCDIR)/UShortToDouble.cc\
	$(SRCDIR)/UShortToFloat.cc\
	$(SRCDIR)/ShortToUChar.cc\
	$(SRCDIR)/ShortToDouble.cc\
	$(SRCDIR)/ShortToFloat.cc\
	$(SRCDIR)/DoubleToUChar.cc\
	$(SRCDIR)/DoubleToFloat.cc\
	$(SRCDIR)/FloatToUChar.cc\
	$(SRCDIR)/FloatToDouble.cc\
	$(SRCDIR)/IntToUChar.cc\
	$(SRCDIR)/IntToDouble.cc\
	$(SRCDIR)/IntToFloat.cc\
	$(SRCDIR)/ULongToUChar.cc\
	$(SRCDIR)/ULongToDouble.cc\
	$(SRCDIR)/ULongToFloat.cc\
	$(SRCDIR)/RGBPixelToVector.cc\
	$(SRCDIR)/Image2DToImage3D.cc\
	$(SRCDIR)/Image3DToImage2D.cc\
#[INSERT NEW CODE FILE HERE]

PSELIBS := Packages/Insight/Core/Datatypes \
	Core/Datatypes Dataflow/Network Dataflow/Ports \
        Core/Persistent Core/Containers Core/Util \
        Core/Exceptions Core/Thread Core/GuiInterface \
        Core/GeomInterface Core/Geom Core/Datatypes Core/Geometry \
        Core/TkExtensions
LIBS := $(TK_LIBRARY) $(GL_LIBRARY) $(M_LIBRARY) $(INSIGHT_LIBRARY)

include $(SCIRUN_SCRIPTS)/smallso_epilogue.mk


