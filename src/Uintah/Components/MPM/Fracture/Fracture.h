#ifndef __FRACTURE_H__
#define __FRACTURE_H__

#include <Uintah/Grid/SimulationState.h>
#include <Uintah/Grid/SimulationStateP.h>
#include <Uintah/Interface/ProblemSpec.h>
#include <Uintah/Interface/ProblemSpecP.h>

#include <Uintah/Interface/DataWarehouseP.h>

namespace Uintah {

   class VarLabel;
   class ProcessorContext;
   class Region;

namespace MPM {

class Fracture {
public:
  void   materialDefectsInitialize(const Region* region,
                                   DataWarehouseP& new_dw);
  
  void   initializeFracture(const Region* region,
                           DataWarehouseP& new_dw);
  
  void   updateSurfaceNormalOfBoundaryParticle(
           const ProcessorContext*,
           const Region* region,
           DataWarehouseP& old_dw,
           DataWarehouseP& new_dw);
  
  void   labelSelfContactNodesAndCells (
           const ProcessorContext*,
           const Region* region,
           DataWarehouseP& old_dw,
           DataWarehouseP& new_dw);

  void   updateParticleInformationInContactCells (
           const ProcessorContext*,
           const Region* region,
           DataWarehouseP& old_dw,
           DataWarehouseP& new_dw);

  void   updateNodeInformationInContactCells (
           const ProcessorContext*,
           const Region* region,
           DataWarehouseP& old_dw,
           DataWarehouseP& new_dw);

  void   crackGrow(
           const ProcessorContext*,
           const Region* region,
           DataWarehouseP& old_dw,
           DataWarehouseP& new_dw);

         Fracture();
	 Fracture(ProblemSpecP& ps, SimulationStateP& d_sS);
                
private:
  double           d_averageMicrocrackLength;
  double           d_materialToughness;
  SimulationStateP d_sharedState;
};

} //namespace MPM
} //namespace Uintah

#endif //__FRACTURE_H__

// $Log$
// Revision 1.9  2000/05/30 04:36:46  tan
// Using MPMLabel instead of VarLabel.
//
// Revision 1.8  2000/05/15 18:58:53  tan
// Initialized NCVariables and CCVaribles for Fracture.
//
// Revision 1.7  2000/05/12 01:46:21  tan
// Added initializeFracture linked to SerialMPM's actuallyInitailize.
//
// Revision 1.6  2000/05/11 20:10:18  dav
// adding MPI stuff.  The biggest change is that old_dws cannot be const and so a large number of declarations had to change.
//
// Revision 1.5  2000/05/10 18:32:22  tan
// Added member funtion to label self-contact cells.
//
// Revision 1.4  2000/05/10 05:04:39  tan
// Basic structure of fracture class.
//
