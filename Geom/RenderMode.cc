
/*
 * RenderMode.cc: RenderMode objects
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   October 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <Geom/RenderMode.h>
#include <Classlib/NotFinished.h>
#include <Geom/Tri.h>
#include <Geometry/BBox.h>
#include <Malloc/Allocator.h>
#include <Math/TrigTable.h>
#include <Math/Trig.h>

GeomRenderMode::GeomRenderMode(DrawType drawtype, GeomObj* child)
: GeomContainer(child), drawtype(drawtype)
{
}

GeomRenderMode::GeomRenderMode(const GeomRenderMode& copy)
: GeomContainer(copy), drawtype(copy.drawtype)
{
}

GeomRenderMode::~GeomRenderMode()
{
    if(child)
	delete child;
}

GeomObj* GeomRenderMode::clone()
{
    return scinew GeomRenderMode(*this);
}

void GeomRenderMode::make_prims(Array1<GeomObj*>& free,
				Array1<GeomObj*>& dontfree)
{
    if(child)
	child->make_prims(free, dontfree);
}

void GeomRenderMode::intersect(const Ray&, Material*,
			       Hit&)
{
    NOT_FINISHED("GeomRenderMode::intersect");
}
