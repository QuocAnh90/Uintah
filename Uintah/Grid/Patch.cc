/* REFERENCED */
static char *id="@(#) $Id$";

#include <Uintah/Grid/Patch.h>
#include <Uintah/Exceptions/InvalidGrid.h>
#include <Uintah/Grid/CellIterator.h>
#include <Uintah/Grid/NodeIterator.h>
#include <Uintah/Grid/NodeSubIterator.h>
#include <Uintah/Grid/SubPatch.h>
#include <Uintah/Math/Primes.h>

#include <SCICore/Exceptions/InternalError.h>
#include <SCICore/Thread/AtomicCounter.h>

#include <values.h>
#include <iostream>
#include <sstream>

using namespace Uintah;
using namespace SCICore::Geometry;
using namespace std;
using SCICore::Exceptions::InternalError;
using SCICore::Math::Floor;
static SCICore::Thread::AtomicCounter ids("Patch ID counter");
using SCICore::Geometry::Min;
using SCICore::Geometry::Max;

Patch::Patch(const Level* level,
	     const IntVector& lowIndex, const IntVector& highIndex,
	     int id)
    : d_level(level), d_lowIndex(lowIndex), d_highIndex(highIndex),
      d_id( id )
{
   if(d_id == -1)
      d_id = ids++;
}

Patch::~Patch()
{
}

bool Patch::findCell(const Point& pos, IntVector& ci) const
{
   ci = d_level->getCellIndex(pos);
   return containsCell(ci);
}

void Patch::findCellsFromNode( const IntVector& nodeIndex,
                               IntVector cellIndex[8]) const
{
   int ix = nodeIndex.x();
   int iy = nodeIndex.y();
   int iz = nodeIndex.z();

   cellIndex[0] = IntVector(ix, iy, iz);
   cellIndex[1] = IntVector(ix, iy, iz-1);
   cellIndex[2] = IntVector(ix, iy-1, iz);
   cellIndex[3] = IntVector(ix, iy-1, iz-1);
   cellIndex[4] = IntVector(ix-1, iy, iz);
   cellIndex[5] = IntVector(ix-1, iy, iz-1);
   cellIndex[6] = IntVector(ix-1, iy-1, iz);
   cellIndex[7] = IntVector(ix-1, iy-1, iz-1);
}

void Patch::findNodesFromCell( const IntVector& cellIndex,
                               IntVector nodeIndex[8]) const
{
   int ix = cellIndex.x();
   int iy = cellIndex.y();
   int iz = cellIndex.z();

   nodeIndex[0] = IntVector(ix, iy, iz);
   nodeIndex[1] = IntVector(ix, iy, iz+1);
   nodeIndex[2] = IntVector(ix, iy+1, iz);
   nodeIndex[3] = IntVector(ix, iy+1, iz+1);
   nodeIndex[4] = IntVector(ix+1, iy, iz);
   nodeIndex[5] = IntVector(ix+1, iy, iz+1);
   nodeIndex[6] = IntVector(ix+1, iy+1, iz);
   nodeIndex[7] = IntVector(ix+1, iy+1, iz+1);
}

bool Patch::findCellAndWeights(const Point& pos,
				IntVector ni[8], double S[8]) const
{
   Point cellpos = d_level->positionToIndex(pos);
   int ix = Floor(cellpos.x());
   int iy = Floor(cellpos.y());
   int iz = Floor(cellpos.z());
   ni[0] = IntVector(ix, iy, iz);
   ni[1] = IntVector(ix, iy, iz+1);
   ni[2] = IntVector(ix, iy+1, iz);
   ni[3] = IntVector(ix, iy+1, iz+1);
   ni[4] = IntVector(ix+1, iy, iz);
   ni[5] = IntVector(ix+1, iy, iz+1);
   ni[6] = IntVector(ix+1, iy+1, iz);
   ni[7] = IntVector(ix+1, iy+1, iz+1);
   double fx = cellpos.x() - ix;
   double fy = cellpos.y() - iy;
   double fz = cellpos.z() - iz;
   double fx1 = 1-fx;
   double fy1 = 1-fy;
   double fz1 = 1-fz;
   S[0] = fx1 * fy1 * fz1;
   S[1] = fx1 * fy1 * fz;
   S[2] = fx1 * fy * fz1;
   S[3] = fx1 * fy * fz;
   S[4] = fx * fy1 * fz1;
   S[5] = fx * fy1 * fz;
   S[6] = fx * fy * fz1;
   S[7] = fx * fy * fz;
   return ix>= d_lowIndex.x()-1 && iy>=d_lowIndex.y()-1 && iz>=d_lowIndex.z()-1 && ix<d_highIndex.x() && iy<d_highIndex.y() && iz<d_highIndex.z();
}


bool Patch::findCellAndShapeDerivatives(const Point& pos,
					 IntVector ni[8],
					 Vector d_S[8]) const
{
   Point cellpos = d_level->positionToIndex(pos);
   int ix = Floor(cellpos.x());
   int iy = Floor(cellpos.y());
   int iz = Floor(cellpos.z());
   ni[0] = IntVector(ix, iy, iz);
   ni[1] = IntVector(ix, iy, iz+1);
   ni[2] = IntVector(ix, iy+1, iz);
   ni[3] = IntVector(ix, iy+1, iz+1);
   ni[4] = IntVector(ix+1, iy, iz);
   ni[5] = IntVector(ix+1, iy, iz+1);
   ni[6] = IntVector(ix+1, iy+1, iz);
   ni[7] = IntVector(ix+1, iy+1, iz+1);
   double fx = cellpos.x() - ix;
   double fy = cellpos.y() - iy;
   double fz = cellpos.z() - iz;
   double fx1 = 1-fx;
   double fy1 = 1-fy;
   double fz1 = 1-fz;
   d_S[0] = Vector(- fy1 * fz1, -fx1 * fz1, -fx1 * fy1);
   d_S[1] = Vector(- fy1 * fz,  -fx1 * fz,   fx1 * fy1);
   d_S[2] = Vector(- fy  * fz1,  fx1 * fz1, -fx1 * fy);
   d_S[3] = Vector(- fy  * fz,   fx1 * fz,   fx1 * fy);
   d_S[4] = Vector(  fy1 * fz1, -fx  * fz1, -fx  * fy1);
   d_S[5] = Vector(  fy1 * fz,  -fx  * fz,   fx  * fy1);
   d_S[6] = Vector(  fy  * fz1,  fx  * fz1, -fx  * fy);
   d_S[7] = Vector(  fy  * fz,   fx  * fz,   fx  * fy);
   return ix>= d_lowIndex.x()-1 && iy>=d_lowIndex.y()-1 && iz>=d_lowIndex.z()-1 && ix<d_highIndex.x() && iy<d_highIndex.y() && iz<d_highIndex.z();
}

ostream& operator<<(ostream& out, const Patch* r)
{
   out << "(Patch: box=" << r->getBox() << ", lowIndex=" << r->getCellLowIndex() << ", highIndex=" << r->getCellHighIndex() << ")";
  return out;
}

long Patch::totalCells() const
{
   IntVector res(d_highIndex-d_lowIndex);
   return res.x()*res.y()*res.z();
}

void Patch::performConsistencyCheck() const
{
   IntVector res(d_highIndex-d_lowIndex);
   if(res.x() < 1 || res.y() < 1 || res.z() < 1) {
      ostringstream msg;
      msg << "Degenerate patch: " << toString() << " (resolution=" << res << ")";
      throw InvalidGrid( msg.str() );
  }
}

Patch::BCType 
Patch::getBCType(Patch::FaceType face) const
{
  return d_bctypes[face];
}

void
Patch::setBCType(Patch::FaceType face, BCType newbc)
{
   d_bctypes[face]=newbc;
}

void
Patch::getFace(FaceType face, int offset, IntVector& l, IntVector& h) const
{
   l=getCellLowIndex();
   h=getCellHighIndex();
   switch(face){
   case xminus:
      l.x(l.x()-offset);
      h.x(l.x()+1-offset);
      break;
   case xplus:
      l.x(h.x()-1+offset);
      h.x(h.x()+offset);
      break;
   case yminus:
      l.y(l.y()-offset);
      h.y(l.y()+1-offset);
      break;
   case yplus:
      l.y(h.y()-1+offset);
      h.y(h.y()+offset);
      break;
   case zminus:
      l.z(l.z()-offset);
      h.z(l.z()+1-offset);
      break;
   case zplus:
      l.z(h.z()-1+offset);
      h.z(h.z()+offset);
      break;
   }
}

string
Patch::toString() const
{
  char str[ 1024 ];

  Box box(getBox());
  sprintf( str, "[ [%2.2lf, %2.2lf, %2.2lf] [%2.2lf, %2.2lf, %2.2lf] ]",
	   box.lower().x(), box.lower().y(), box.lower().z(),
	   box.upper().x(), box.upper().y(), box.upper().z() );

  return string( str );
}

CellIterator
Patch::getCellIterator(const Box& b) const
{
   Point l = d_level->positionToIndex(b.lower());
   Point u = d_level->positionToIndex(b.upper());
   IntVector low((int)l.x(), (int)l.y(), (int)l.z());
   IntVector high(RoundUp(u.x()), RoundUp(u.y()), RoundUp(u.z()));
   low = SCICore::Geometry::Max(low, getCellLowIndex());
   high = SCICore::Geometry::Min(high, getCellHighIndex());
   return CellIterator(low, high);
}

Box Patch::getGhostBox(const IntVector& lowOffset,
		       const IntVector& highOffset) const
{
   return Box(d_level->getNodePosition(d_lowIndex+lowOffset),
	      d_level->getNodePosition(d_highIndex+highOffset));
}

NodeIterator Patch::getNodeIterator() const
{
   return NodeIterator(getNodeLowIndex(), getNodeHighIndex());
}

IntVector Patch::getNodeHighIndex() const
{
   IntVector h(d_highIndex+
	       IntVector(getBCType(xplus) == Neighbor?0:1,
			 getBCType(yplus) == Neighbor?0:1,
			 getBCType(zplus) == Neighbor?0:1));
   return h;
}

//
// $Log$
// Revision 1.11  2000/06/15 21:57:19  sparker
// Added multi-patch support (bugzilla #107)
// Changed interface to datawarehouse for particle data
// Particles now move from patch to patch
//
// Revision 1.10  2000/06/14 19:58:03  guilkey
// Added a different version of findCell.
//
// Revision 1.9  2000/06/13 21:28:30  jas
// Added missing TypeUtils.h for fun_forgottherestofname and copy constructor
// was wrong for CellIterator.
//
// Revision 1.8  2000/06/08 17:47:47  dav
// longer error message
//
// Revision 1.7  2000/06/07 18:31:00  tan
// Requirement for getHighGhostCellIndex() and getLowGhostCellIndex()
// cancelled.
//
// Revision 1.6  2000/06/05 19:25:14  tan
// I need the following two functions,
// (1) IntVector getHighGhostCellIndex() const;
// (2) IntVector getLowGhostCellIndex() const;
// The temporary empty functions are created.
//
// Revision 1.5  2000/06/04 04:36:07  tan
// Added function findNodesFromCell() to find the 8 neighboring node indexes
// according to a given cell index.
//
// Revision 1.4  2000/06/02 20:44:56  tan
// Corrected a mistake in function findCellsFromNode().
//
// Revision 1.3  2000/06/02 19:58:01  tan
// Added function findCellsFromNode() to find the 8 neighboring cell
// indexes according to a given node index.
//
// Revision 1.2  2000/06/01 22:14:06  tan
// Added findCell(const Point& pos).
//
// Revision 1.1  2000/05/30 20:19:31  sparker
// Changed new to scinew to help track down memory leaks
// Changed region to patch
//
// Revision 1.20  2000/05/28 17:25:06  dav
// adding mpi stuff
//
// Revision 1.19  2000/05/20 08:09:26  sparker
// Improved TypeDescription
// Finished I/O
// Use new XML utility libraries
//
// Revision 1.18  2000/05/15 19:39:49  sparker
// Implemented initial version of DataArchive (output only so far)
// Other misc. cleanups
//
// Revision 1.17  2000/05/10 20:03:02  sparker
// Added support for ghost cells on node variables and particle variables
//  (work for 1 patch but not debugged for multiple)
// Do not schedule fracture tasks if fracture not enabled
// Added fracture directory to MPM sub.mk
// Be more uniform about using IntVector
// Made patches have a single uniform index space - still needs work
//
// Revision 1.16  2000/05/09 03:24:39  jas
// Added some enums for grid boundary conditions.
//
// Revision 1.15  2000/05/07 06:02:12  sparker
// Added beginnings of multiple patch support and real dependencies
//  for the scheduler
//
// Revision 1.14  2000/05/05 06:42:45  dav
// Added some _hopefully_ good code mods as I work to get the MPI stuff to work.
//
// Revision 1.13  2000/05/04 19:06:48  guilkey
// Added the beginnings of grid boundary conditions.  Functions still
// need to be filled in.
//
// Revision 1.12  2000/05/02 20:30:59  jas
// Fixed the findCellAndShapeDerivatives.
//
// Revision 1.11  2000/05/02 20:13:05  sparker
// Implemented findCellAndWeights
//
// Revision 1.10  2000/05/02 06:07:23  sparker
// Implemented more of DataWarehouse and SerialMPM
//
// Revision 1.9  2000/04/28 20:24:44  jas
// Moved some private copy constructors to public for linux.  Velocity
// field is now set from the input file.  Simulation state now correctly
// determines number of velocity fields.
//
// Revision 1.8  2000/04/28 03:58:20  sparker
// Fixed countParticles
// Implemented createParticles, which doesn't quite work yet because the
//   data warehouse isn't there yet.
// Reduced the number of particles in the bar problem so that it will run
//   quickly during development cycles
//
// Revision 1.7  2000/04/27 23:18:50  sparker
// Added problem initialization for MPM
//
// Revision 1.6  2000/04/26 06:48:54  sparker
// Streamlined namespaces
//
// Revision 1.5  2000/04/13 06:51:01  sparker
// More implementation to get this to work
//
// Revision 1.4  2000/04/12 23:00:49  sparker
// Starting problem setup code
// Other compilation fixes
//
// Revision 1.3  2000/03/16 22:08:01  dav
// Added the beginnings of cocoon docs.  Added namespaces.  Did a few other coding standards updates too
//
//
