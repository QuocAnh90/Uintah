/*
 *  ClipField.cc:  Unfinished modules
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   February 1995
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <Classlib/NotFinished.h>
#include <Dataflow/Module.h>
#include <Datatypes/ScalarFieldPort.h>
#include <Datatypes/SurfacePort.h>
#include <Datatypes/ScalarField.h>
#include <Datatypes/ScalarFieldRG.h>
#include <Datatypes/ScalarFieldRGBase.h>
#include <Datatypes/ScalarFieldRGdouble.h>
#include <Datatypes/ScalarFieldRGfloat.h>
#include <Datatypes/ScalarFieldRGint.h>
#include <Datatypes/ScalarFieldRGchar.h>
#include <Geometry/Point.h>
#include <TCL/TCLvar.h>
#include <stdio.h>

class ClipField : public Module {
    ScalarFieldIPort* ifield;
    ScalarFieldOPort* ofield;
public:
    TCLint x_min;
    TCLint x_max;
    TCLint y_min;
    TCLint y_max;
    TCLint z_min;
    TCLint z_max;
    TCLint sameInput;
    int first_time;
    int last_x_min;
    int last_y_min;
    int last_z_min;
    int last_x_max;
    int last_y_max;
    int last_z_max;
    ClipField(const clString& id);
    ClipField(const ClipField&, int deep);
    virtual ~ClipField();
    virtual Module* clone(int deep);
    virtual void execute();
    ScalarFieldHandle fldHandle;
    ScalarFieldRGBase* osf;
};

extern "C" {
Module* make_ClipField(const clString& id)
{
    return new ClipField(id);
}
};

ClipField::ClipField(const clString& id)
: Module("ClipField", id, Filter), 
  x_min("x_min", id, this),y_min("y_min", id, this),z_min("z_min", id, this), 
  x_max("x_max", id, this),y_max("y_max", id, this),z_max("z_max", id, this),
  sameInput("sameInput", id, this)
{
    ifield=new ScalarFieldIPort(this, "Geometry", ScalarFieldIPort::Atomic);
    add_iport(ifield);
    // Create the output port
    ofield=new ScalarFieldOPort(this, "Geometry", ScalarFieldIPort::Atomic);
    add_oport(ofield);
    first_time=1;
}

ClipField::ClipField(const ClipField& copy, int deep)
: Module(copy, deep),
  x_min("x_min", id, this),y_min("y_min", id, this),z_min("z_min", id, this), 
  x_max("x_max", id, this),y_max("y_max", id, this),z_max("z_max", id, this),
  sameInput("sameInput", id, this)
{
}

ClipField::~ClipField()
{
}

Module* ClipField::clone(int deep)
{
    return new ClipField(*this, deep);
}

void ClipField::execute()
{
    ScalarFieldHandle ifh;
    if(!ifield->get(ifh))
	return;
    ScalarFieldRGBase* isf=ifh->getRGBase();
    if(!isf){
	error("ClipField can't deal with unstructured grids!");
	return;
    }

    ScalarFieldRGdouble *ifd=isf->getRGDouble();
    ScalarFieldRGfloat *iff=isf->getRGFloat();
    ScalarFieldRGint *ifi=isf->getRGInt();
    ScalarFieldRGchar *ifc=isf->getRGChar();
    
    int mxx, mxy, mxz, mnx, mny, mnz;
    mxx=x_max.get()-1;
    mxy=y_max.get()-1;
    mxz=z_max.get()-1;
    mnx=x_min.get()-1;
    mny=y_min.get()-1;
    mnz=z_min.get()-1;
    if (!sameInput.get() || first_time || 
	mxx!=last_x_max || mxy!=last_y_max || mxz!=last_z_max ||
	mnx!=last_x_min || mny!=last_y_min || mnz != last_z_min) {
	first_time=0;
	if (ifd) {
	    ScalarFieldRGdouble *of;
	    fldHandle = of = 0;
	    fldHandle = of = new ScalarFieldRGdouble;
	    of->resize(mxx-mnx+1, mxy-mny+1, mxz-mnz+1);
	    for (int i=0; i<=mxx-mnx; i++) {
		for (int j=0; j<=mxy-mny; j++) {
		    for (int k=0; k<=mxz-mnz; k++) {
			of->grid(i,j,k)=ifd->grid(i+mnx, j+mny, k+mnz);
		    }
		}
	    }
	    of->compute_minmax();
	    osf=of;
	} else if (iff) {
	    ScalarFieldRGfloat *of;
	    fldHandle = of = 0;
	    fldHandle = of = new ScalarFieldRGfloat;
	    of->resize(mxx-mnx+1, mxy-mny+1, mxz-mnz+1);
	    for (int i=0; i<=mxx-mnx; i++) {
		for (int j=0; j<=mxy-mny; j++) {
		    for (int k=0; k<=mxz-mnz; k++) {
			of->grid(i,j,k)=iff->grid(i+mnx, j+mny, k+mnz);
		    }
		}
	    }
	    of->compute_minmax();
	    osf=of;
	} else if (ifi) {
	    ScalarFieldRGint *of;
	    fldHandle = of = 0;
	    fldHandle = of = new ScalarFieldRGint;
	    of->resize(mxx-mnx+1, mxy-mny+1, mxz-mnz+1);
	    for (int i=0; i<=mxx-mnx; i++) {
		for (int j=0; j<=mxy-mny; j++) {
		    for (int k=0; k<=mxz-mnz; k++) {
			of->grid(i,j,k)=ifi->grid(i+mnx, j+mny, k+mnz);
		    }
		}
	    }
	    of->compute_minmax();
	    osf=of;
	} else {
	    ScalarFieldRGchar *of;
	    fldHandle = of = 0;
	    fldHandle = of = new ScalarFieldRGchar;
	    of->resize(mxx-mnx+1, mxy-mny+1, mxz-mnz+1);
	    for (int i=0; i<=mxx-mnx; i++) {
		for (int j=0; j<=mxy-mny; j++) {
		    for (int k=0; k<=mxz-mnz; k++) {
			of->grid(i,j,k)=ifc->grid(i+mnx, j+mny, k+mnz);
		    }
		}
	    }
	    of->compute_minmax();
	    osf=of;
	}
	osf->set_bounds(Point(mnx-1, mny-1, mnz-1), 
			Point(mxx-1, mxy-1, mxz-1));
	last_x_max=mxx;
	last_y_max=mxy;
	last_z_max=mxz;
	last_x_min=mnx;
	last_y_min=mny;
	last_z_min=mnz;
    }
    ofield->send(osf);
}
