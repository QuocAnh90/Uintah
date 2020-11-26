/*
 * The MIT License
 *
 * Copyright (c) 1997-2020 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef UINTAH_HOMEBREW_SEGEMPM_H
#define UINTAH_HOMEBREW_SEGEMPM_H

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/SwitchingCriteria.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Variables/ComputeSet.h>
// put here to avoid template problems
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Materials/Contact/Contact.h>
#include <CCA/Components/MPM/MPMCommon.h>
#include <Core/Geometry/Vector.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBC.h>
#include <CCA/Components/MPM/PhysicalBC/LoadCurve.h>
#include <CCA/Components/OnTheFlyAnalysis/AnalysisModule.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <CCA/Components/MPM/Materials/MPMMaterial.h>
#include <CCA/Components/MPM/SerialMPM.h>
#include <CCA/Components/MPM/SegeMPM.h>


namespace Uintah {
class AnalysisModule;
class MPM;
class MPMLabel;
class Output;

/**************************************

CLASS
   SerialMPM
   
   Short description...

GENERAL INFORMATION

   SerialMPM.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   
KEYWORDS
   SerialMPM

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/

  class SegeMPM : public MPMCommon {
   public:
    SegeMPM(const ProcessorGroup* myworld,
              const MaterialManagerP materialManager);
    virtual ~SegeMPM();

  //////////
  // Insert Documentation Here:

  virtual void problemSetup(const ProblemSpecP& params, 
                            const ProblemSpecP& restart_prob_spec,
                            GridP&);

  virtual void outputProblemSpec(ProblemSpecP& ps);

  virtual void scheduleInitialize(const LevelP& level,
                                  SchedulerP&);
                                  
  virtual void scheduleDeleteGeometryObjects(const LevelP& level,
                                             SchedulerP& sched);

  virtual void scheduleRestartInitialize(const LevelP& level,
                                         SchedulerP& sched);

  void schedulePrintParticleCount(const LevelP& level, SchedulerP& sched);
  
  void scheduleTotalParticleCount(SchedulerP& sched,
                                 const PatchSet* patches,
                                 const MaterialSet* matls);
  //////////
  // Insert Documentation Here:
  virtual void scheduleComputeStableTimeStep(const LevelP& level, SchedulerP&);

  //////////
  // Insert Documentation Here:
  virtual void scheduleTimeAdvance(const LevelP& level, SchedulerP&);

   protected:
  //////////
  // Insert Documentation Here:
  friend class MPMICE;
  friend class MPM;
  friend class MPMICE2;
  friend class MPMArches;

  MPMFlags* flags;

  SerialMPM* d_mpm;
  MPMLabel* Mlb;

  void scheduleInitializePressureBCs(const LevelP& level, SchedulerP&);

  //__________________________________
  // refinement criteria threshold knobs
  struct thresholdVar {
    std::string name;
    int matl;
    double value;
  };
  std::vector<thresholdVar> d_thresholdVars;
 
  double           d_nextOutputTime;
  double           d_SMALL_NUM_MPM;
  int              NGP;      // Number of ghost particles needed.
  int              NGN;      // Number of ghost nodes     needed.
  int              d_8or27;

  std::vector<AnalysisModule*> d_analysisModules;

   private:

  SegeMPM(const SegeMPM&);
  SegeMPM& operator=(const SegeMPM&);
};
      
} // end namespace Uintah

#endif
