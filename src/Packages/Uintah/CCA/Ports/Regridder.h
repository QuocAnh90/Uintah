#ifndef UINTAH_HOMEBREW_REGRIDDER_H
#define UINTAH_HOMEBREW_REGRIDDER_H

#include <Packages/Uintah/Core/Parallel/UintahParallelPort.h>
#include <Packages/Uintah/Core/Grid/GridP.h>
#include <Packages/Uintah/Core/Grid/LevelP.h>
#include <Packages/Uintah/Core/Grid/SimulationStateP.h>
#include <Packages/Uintah/Core/ProblemSpec/ProblemSpecP.h>
#include <Packages/Uintah/CCA/Ports/SchedulerP.h>

#include <Packages/Uintah/CCA/Ports/share.h>

namespace Uintah {

/**************************************

CLASS
   Regridder
   
   Short description...

GENERAL INFORMATION

   Regridder.h

   Bryan Worthen
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   Copyright (C) 2000 SCI Group

KEYWORDS
   Regridder

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/

  //! Takes care of AMR Regridding.
  class SCISHARE Regridder : public UintahParallelPort {
  public:
    Regridder();
    virtual ~Regridder();

    //! Initialize with regridding parameters from ups file
    virtual void problemSetup(const ProblemSpecP& params, const GridP&,
			      const SimulationStateP& state) = 0;

    //! Asks if we need to recompile the task graph.
    virtual bool needRecompile(double time, double delt,
			       const GridP& grid) = 0;

    //! Do we need to regrid this timestep?
    virtual bool needsToReGrid() = 0;

    //! Asks if we are going to do regridding
    virtual bool isAdaptive() = 0;

    //! Schedules task to initialize the error flags to 0
    virtual void scheduleInitializeErrorEstimate(const LevelP& level) = 0;

    //! Schedules task to dilate existing error flags
    virtual void scheduleDilation(const LevelP& level) = 0;

    //! Asks if we are going to do regridding
    virtual bool flaggedCellsOnFinestLevel(const GridP& grid) = 0;

    //! Returns the max number of levels this regridder will store
    virtual int maxLevels() = 0;

    //! Create a new Grid
    virtual Grid* regrid(Grid* oldGrid) = 0;

    //! If the Regridder set up the load balance in the process of Regridding
    virtual bool isLoadBalanced() { return false; }
  private:
    Regridder(const Regridder&);
    Regridder& operator=(const Regridder&);
  };

} // End namespace Uintah

#endif
