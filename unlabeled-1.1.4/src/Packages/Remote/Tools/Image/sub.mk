#
# Makefile fragment for this subdirectory
# $Id$
#

SRCDIR   := Remote/Tools/Image

SRCS     += \
	$(SRCDIR)/Bmp.cc \
	$(SRCDIR)/Filter.cc \
	$(SRCDIR)/Gif.cc \
	$(SRCDIR)/Image.cc \
	$(SRCDIR)/Quant.cc

#
# $Log$
# Revision 1.1.4.1  2000/10/26 23:24:02  moulding
# merge HEAD into FIELD_REDESIGN
#
# Revision 1.1  2000/07/10 20:39:00  dahart
# initial add
#
# Revision 1.1  2000/06/06 15:29:31  dahart
# - Added a new package / directory tree for the remote visualization
# framework.
# - Added a new Salmon-derived module with a small server-daemon,
# network support, and a couple of very basic remote vis. services.
#
# Revision 1.1  2000/03/17 09:27:18  sparker
# New makefile scheme: sub.mk instead of Makefile.in
# Use XML-based files for module repository
# Plus many other changes to make these two things work
#
#
