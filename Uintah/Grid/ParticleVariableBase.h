
#ifndef UINTAH_HOMEBREW_ParticleVariableBase_H
#define UINTAH_HOMEBREW_ParticleVariableBase_H

#include <vector>
#include <Uintah/Grid/ParticleSubset.h>

namespace Uintah {
   class OutputContext;
   class ParticleSubset;
   class ParticleSet;
   class Patch;

/**************************************

CLASS
   ParticleVariableBase
   
   Short description...

GENERAL INFORMATION

   ParticleVariableBase.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   Copyright (C) 2000 SCI Group

KEYWORDS
   ParticleVariableBase

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/

   class ParticleVariableBase {
   public:
      
      virtual ~ParticleVariableBase();
      virtual void copyPointer(const ParticleVariableBase&) = 0;
      
      //////////
      // Insert Documentation Here:
      virtual ParticleVariableBase* clone() const = 0;

      virtual void allocate(ParticleSubset*) = 0;
      virtual void gather(ParticleSubset* dest,
			  std::vector<ParticleSubset*> subsets,
			  std::vector<ParticleVariableBase*> srcs) = 0;
      virtual void emit(OutputContext&) = 0;

      //////////
      // Insert Documentation Here:
      ParticleSubset* getParticleSubset() const {
	 return d_pset;
      }

      //////////
      // Insert Documentation Here:
      ParticleSet* getParticleSet() const {
	 return d_pset->getParticleSet();
      }
      
   protected:
      ParticleVariableBase(const ParticleVariableBase&);
      ParticleVariableBase(ParticleSubset* pset);
      ParticleVariableBase& operator=(const ParticleVariableBase&);
      
      ParticleSubset*  d_pset;

   private:
   };
   
} // end namespace Uintah

//
// $Log$
// Revision 1.8  2000/06/15 21:57:19  sparker
// Added multi-patch support (bugzilla #107)
// Changed interface to datawarehouse for particle data
// Particles now move from patch to patch
//
// Revision 1.7  2000/05/30 20:19:31  sparker
// Changed new to scinew to help track down memory leaks
// Changed region to patch
//
// Revision 1.6  2000/05/15 19:39:48  sparker
// Implemented initial version of DataArchive (output only so far)
// Other misc. cleanups
//
// Revision 1.5  2000/05/10 20:03:02  sparker
// Added support for ghost cells on node variables and particle variables
//  (work for 1 patch but not debugged for multiple)
// Do not schedule fracture tasks if fracture not enabled
// Added fracture directory to MPM sub.mk
// Be more uniform about using IntVector
// Made patches have a single uniform index space - still needs work
//
// Revision 1.4  2000/05/01 16:18:18  sparker
// Completed more of datawarehouse
// Initial more of MPM data
// Changed constitutive model for bar
//
// Revision 1.3  2000/04/28 07:35:37  sparker
// Started implementation of DataWarehouse
// MPM particle initialization now works
//
// Revision 1.2  2000/04/26 06:48:52  sparker
// Streamlined namespaces
//
// Revision 1.1  2000/04/20 20:09:22  jas
// I don't know what these do, but Steve says we need them.
//
// Revision 1.1  2000/04/19 05:26:15  sparker
// Implemented new problemSetup/initialization phases
// Simplified DataWarehouse interface (not finished yet)
// Made MPM get through problemSetup, but still not finished
//
//

#endif
