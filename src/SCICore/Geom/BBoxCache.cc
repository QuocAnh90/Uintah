//static char *id="@(#) $Id$"

/*
 *  BBoxCache.cc: ?
 *
 *  Written by:
 *   Author?
 *   Department of Computer Science
 *   University of Utah
 *   Date?
 *
 *  Copyright (C) 199? SCI Group
 */

#include <SCICore/Geom/BBoxCache.h>
#include <SCICore/Malloc/Allocator.h>
#include <iostream.h>

namespace SCICore {
namespace GeomSpace {

Persistent* make_GeomBBoxCache()
{
    return scinew GeomBBoxCache(0);
}

PersistentTypeID GeomBBoxCache::type_id("GeomBBoxCache", "GeomObj",
					make_GeomBBoxCache);


GeomBBoxCache::GeomBBoxCache(GeomObj* obj)
:child(obj),bbox_cached(0)
{

}

GeomBBoxCache::GeomBBoxCache(GeomObj* obj, BBox &box)
:child(obj),bbox_cached(1)
{
  bbox.extend( box );
}

GeomBBoxCache::~GeomBBoxCache()
{
    if(child)
	delete child;
}

GeomObj* GeomBBoxCache::clone()
{
    cerr << "GeomBBoxCache::clone not implemented!\n";
    return 0;
}

void GeomBBoxCache::reset_bbox()
{
    bbox_cached = 0;
}

void GeomBBoxCache::get_bounds(BBox& box)
{
    if (!bbox_cached) {
	bbox.reset();
	child->get_bounds(bbox);
	bbox_cached = 1;
    }

    box.extend( bbox );
}

#define GEOMBBOXCACHE_VERSION 2

void GeomBBoxCache::io(Piostream& stream)
{
    using SCICore::PersistentSpace::Pio;

    int version=stream.begin_class("GeomBBoxCache", GEOMBBOXCACHE_VERSION);
    Pio(stream, bbox_cached);
    if(version < 2){
	int bsphere_cached;
	Pio(stream, bsphere_cached);
    }
    Pio(stream, bbox);
    if(version < 2){
	// Old BSphere stuff...
	stream.begin_cheap_delim();
	int have_some;
	Pio(stream, have_some);
	Point cen;
	Pio(stream, cen);
	double rad;
	Pio(stream, rad);
	stream.end_cheap_delim();
    }
    Pio(stream, child);
    stream.end_class();
}

bool GeomBBoxCache::saveobj(ostream& out, const clString& format,
			    GeomSave* saveinfo)
{
    return child->saveobj(out, format, saveinfo);
}

} // End namespace GeomSpace
} // End namespace SCICore

//
// $Log$
// Revision 1.5  1999/08/28 17:54:38  sparker
// Integrated new Thread library
//
// Revision 1.4  1999/08/19 23:18:05  sparker
// Removed a bunch of #include <SCICore/Util/NotFinished.h> statements
// from files that did not need them.
//
// Revision 1.3  1999/08/17 23:50:17  sparker
// Removed all traces of the old Raytracer and X11 renderers.
// Also removed a .o and .d file
//
// Revision 1.2  1999/08/17 06:39:03  sparker
// Merged in modifications from PSECore to make this the new "blessed"
// version of SCIRun/Uintah.
//
// Revision 1.1  1999/07/27 16:56:36  mcq
// Initial commit
//
// Revision 1.2  1999/07/07 21:10:48  dav
// added beginnings of support for g++ compilation
//
// Revision 1.1.1.1  1999/04/24 23:12:20  dav
// Import sources
//
//

