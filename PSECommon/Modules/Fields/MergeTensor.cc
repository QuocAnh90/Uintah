//static char *id="@(#) $Id$";

/*
 *  Gradient.cc:  Unfinished modules
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <SCICore/Util/NotFinished.h>
#include <PSECore/Dataflow/Module.h>
#include <PSECore/CommonDatatypes/ScalarFieldPort.h>
#include <PSECore/CommonDatatypes/SurfacePort.h>
#include <SCICore/CoreDatatypes/ScalarFieldRG.h>
#include <PSECore/CommonDatatypes/VectorFieldPort.h>
#include <SCICore/CoreDatatypes/VectorFieldRG.h>
#include <SCICore/Geometry/Point.h>
#include <SCICore/TclInterface/TCLvar.h>

namespace PSECommon {
namespace Modules {

using namespace PSECore::Dataflow;
using namespace PSECore::CommonDatatypes;
using namespace SCICore::TclInterface;
using namespace SCICore::GeomSpace;

class MergeTensor : public Module {
    ScalarFieldIPort* infield;

    VectorFieldIPort* tensora;
    VectorFieldIPort* tensorb;

    VectorFieldOPort* outfield;
    VectorFieldOPort* outJ;
public:
    MergeTensor(const clString& id);
    MergeTensor(const MergeTensor&, int deep);
    virtual ~MergeTensor();
    virtual Module* clone(int deep);
    virtual void execute();
};

Module* make_MergeTensor(const clString& id) {
  return new MergeTensor(id);
}

MergeTensor::MergeTensor(const clString& id)
: Module("MergeTensor", id, Filter)
{
    infield=new ScalarFieldIPort(this, "Geometry", ScalarFieldIPort::Atomic);
    add_iport(infield);

    tensora = new 
      VectorFieldIPort(this,"Tensor Diagonal",VectorFieldIPort::Atomic);
    add_iport(tensora);

    tensorb = new 
      VectorFieldIPort(this,"Tensor Off Diagonal",VectorFieldIPort::Atomic);
    add_iport(tensorb);

    // Create the output ports
    outfield=new VectorFieldOPort(this, "Vec Field", VectorFieldIPort::Atomic);
    add_oport(outfield);

    outJ=new VectorFieldOPort(this, "J", VectorFieldIPort::Atomic);
    add_oport(outJ);
}

MergeTensor::MergeTensor(const MergeTensor& copy, int deep)
: Module(copy, deep)
{
}

MergeTensor::~MergeTensor()
{
}

Module* MergeTensor::clone(int deep)
{
    return new MergeTensor(*this, deep);
}

void MergeTensor::execute()
{
  ScalarFieldHandle iff;
  if(!infield->get(iff))
    return;
  ScalarFieldRG* sfield=iff->getRG();
  if(!sfield){
    error("Gradient can't deal with this field");
    return;
  }

  // this is all that is needed to create vector field from potentials

  VectorFieldRG* vfield;

  vfield = new VectorFieldRG();

  vfield->resize(sfield->nx,sfield->ny,sfield->nz);
  Point bmin,bmax;
  sfield->get_bounds(bmin,bmax);
  vfield->set_bounds(bmin,bmax);

  int x;
  for(x=0;x<sfield->nx;x++)
    for(int y=0;y<sfield->ny;y++)
      for(int z=0;z<sfield->nz;z++) {
	vfield->grid(x,y,z) = sfield->gradient(x,y,z);
      }
  
  outfield->send(VectorFieldHandle(vfield));

  VectorFieldRG* jfield,*arg,*brg; // J = sigma*grad, where grad is vfield...

  VectorFieldHandle ta,tb;

  if (!tensora->get(ta))
    return;
  if (!tensorb->get(tb))
    return;

  arg = ta->getRG();
  brg = tb->getRG();

  if (!arg || !brg)
    return; // tensors aren't regular grids...

  jfield = new VectorFieldRG();

  jfield->resize(sfield->nx,sfield->ny,sfield->nz);
  sfield->get_bounds(bmin,bmax);
  jfield->set_bounds(bmin,bmax);
  
  // ok, now run through tensor and vectors, multiply and store...

  for(x=0;x<sfield->nx;x++)
    for(int y=0;y<sfield->ny;y++)
      for(int z=0;z<sfield->nz;z++) {
	Vector G = vfield->grid(x,y,z);
	Vector TD = arg->grid(x,y,z);
	Vector TOD = brg->grid(x,y,z);

	// now do matrix vector multiply...

	Vector nvec(TD.x()*G.x()+TOD.x()*G.y()+TOD.y()*G.z(),
		    TOD.x()*G.x() + TD.y()*G.y() + TOD.z()*G.z(),
		    TOD.y()*G.x() + TOD.z()*G.y() + TD.z()*G.z());

	jfield->grid(x,y,z) = nvec;
      }
  outJ->send(VectorFieldHandle(jfield));

}

} // End namespace Modules
} // End namespace PSECommon

//
// $Log$
// Revision 1.2  1999/08/17 06:37:29  sparker
// Merged in modifications from PSECore to make this the new "blessed"
// version of SCIRun/Uintah.
//
// Revision 1.1  1999/07/27 16:57:43  mcq
// Initial commit
//
// Revision 1.2  1999/04/28 20:51:12  dav
// deleted some files that are dependent on DaveW files
//
// Revision 1.1.1.1  1999/04/24 23:12:32  dav
// Import sources
//
//
