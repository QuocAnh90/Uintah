//static char *id="@(#) $Id$";

/*
 *  Geom.cc: Displayable Geometry
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <SCICore/Geom/GeomObj.h>
#include <SCICore/Geometry/Vector.h>

namespace SCICore {
namespace GeomSpace {

PersistentTypeID GeomObj::type_id("GeomObj", "Persistent", 0);

GeomObj::GeomObj(int id) : id(id)
{
}

GeomObj::GeomObj(const GeomObj&)
{
}

GeomObj::~GeomObj()
{
}

void GeomObj::reset_bbox()
{
    // Nothing to do, by default.
}

void GeomObj::io(Piostream&)
{
    // Nothing for now...
}

void Pio( Piostream & stream, GeomObj *& obj )
{
    Persistent* tmp=obj;
    stream.io(tmp, GeomObj::type_id);
    if(stream.reading())
	obj=(GeomObj*)tmp;
}

} // End namespace GeomSpace
} // End namespace SCICore

//
// $Log$
// Revision 1.5  2000/01/03 20:12:36  kuzimmer
//  Forgot to check in these files for picking spheres
//
// Revision 1.4  1999/09/08 02:26:50  sparker
// Various #include cleanups
//
// Revision 1.3  1999/08/17 23:50:22  sparker
// Removed all traces of the old Raytracer and X11 renderers.
// Also removed a .o and .d file
//
// Revision 1.2  1999/08/17 06:39:09  sparker
// Merged in modifications from PSECore to make this the new "blessed"
// version of SCIRun/Uintah.
//
// Revision 1.1  1999/07/27 16:56:40  mcq
// Initial commit
//
// Revision 1.1.1.1  1999/04/24 23:12:22  dav
// Import sources
//
//
