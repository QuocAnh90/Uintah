#include "Lattice.h"

#include "Cell.h"
#include <Uintah/Grid/Patch.h>

namespace Uintah {
namespace MPM {

Lattice::
Lattice(const Patch* patch,const ParticleVariable<Point>& pX)
: Array3<Cell>(patch->getLowGhostCellIndex(),patch->getHighGhostCellIndex()),
  d_patch(patch), d_pX(pX)
{
  ParticleSubset* pset = pX.getParticleSubset();

  for(ParticleSubset::iterator part_iter = pset->begin();
      part_iter != pset->end(); part_iter++)
  {
    (*this)[ patch->findCell(pX[*part_iter]) ].insert(*part_iter);
  }
}

void Lattice::fingNeighborWithout(particleIndex pIndex)
{
//  (*this)[ d_patch->findCell( d_pX[pIndex] ) ];
}

void Lattice::fingNeighborWith(particleIndex pIndex)
{
}

const Patch* Lattice::getPatch()
{
  return d_patch;
}

const ParticleVariable<Point>& Lattice::getParticlesPosition()
{
  return d_pX;
}

  
} //namespace MPM
} //namespace Uintah

// $Log$
// Revision 1.3  2000/06/05 22:32:38  tan
// Added function to find neighbor for a given particle index.
//
// Revision 1.2  2000/06/05 19:49:01  tan
// Added d_lattice which is a Array3 data of cells in a given patch.
//
// Revision 1.1  2000/06/05 17:22:05  tan
// Lattice class will be designed to make it easier to handle the grid/particle
// relationship in a given patch and a given velocity field.
//
