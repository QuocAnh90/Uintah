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
#include <CCA/Components/MPM/SerialMPM.h>
#include <CCA/Components/MPM/SegeMPM.h>

#include <CCA/Components/MPM/Core/MPMDiffusionLabel.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Materials/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/Materials/ConstitutiveModel/PlasticityModels/DamageModel.h>
#include <CCA/Components/MPM/Materials/ConstitutiveModel/PlasticityModels/ErosionModel.h>
#include <CCA/Components/MPM/Materials/MPMMaterial.h>
#include <CCA/Components/MPM/Materials/Contact/Contact.h>
#include <CCA/Components/MPM/Materials/Contact/ContactFactory.h>
#include <CCA/Components/MPM/CohesiveZone/CZMaterial.h>
#include <CCA/Components/MPM/HeatConduction/HeatConduction.h>
#include <CCA/Components/MPM/Materials/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/MMS/MMS.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContactFactory.h>
#include <CCA/Components/OnTheFlyAnalysis/AnalysisModuleFactory.h>
#include <CCA/Components/MPM/SerialMPM.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/PerPatchVars.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Point.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/Matrix3.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/DOUT.hpp>

// Diffusion includes
#include <CCA/Components/MPM/Materials/Diffusion/DiffusionInterfaces/SDInterfaceModel.h>
#include <CCA/Components/MPM/Materials/Diffusion/DiffusionModels/ScalarDiffusionModel.h>
#include <CCA/Components/MPM/Materials/Diffusion/SDInterfaceModelFactory.h>
#include <CCA/Components/MPM/PhysicalBC/FluxBCModelFactory.h>
#include <CCA/Components/MPM/PhysicalBC/ScalarFluxBC.h>
#include <CCA/Components/MPM/PhysicalBC/FluxBCModel.h>


#include <iostream>
#include <fstream>
#include <cmath>

using namespace Uintah;
using namespace std;

static DebugStream cout_doing("MPM", false);
static DebugStream cout_dbg("SegeMPM", false);

static Vector face_norm(Patch::FaceType f)
{
  switch(f) {
  case Patch::xminus: return Vector(-1,0,0);
  case Patch::xplus:  return Vector( 1,0,0);
  case Patch::yminus: return Vector(0,-1,0);
  case Patch::yplus:  return Vector(0, 1,0);
  case Patch::zminus: return Vector(0,0,-1);
  case Patch::zplus:  return Vector(0,0, 1);
  default:            return Vector(0,0,0); // oops !
  }
}

SegeMPM::SegeMPM( const ProcessorGroup* myworld,
                      const MaterialManagerP materialManager) :
  MPMCommon( myworld, materialManager )
{
  flags = scinew MPMFlags(myworld);
  d_mpm = scinew SerialMPM(myworld, m_materialManager);
  Mlb = d_mpm->lb;

  d_nextOutputTime=0.;
  d_SMALL_NUM_MPM=1e-200;
  NGP     = 1;
  NGN     = 1;
}

SegeMPM::~SegeMPM()
{
  delete d_mpm;
  delete flags;

  if (flags->d_doScalarDiffusion) {
    delete d_sdInterfaceModel;
  }

  MPMPhysicalBCFactory::clean();

  if(d_analysisModules.size() != 0){
    vector<AnalysisModule*>::iterator iter;
    for( iter  = d_analysisModules.begin();
         iter != d_analysisModules.end(); iter++){
      AnalysisModule* am = *iter;
      am->releaseComponents();
      delete am;
    }
  }

  if(d_switchCriteria) {
    delete d_switchCriteria;
  }
}

void SegeMPM::problemSetup(const ProblemSpecP& prob_spec,
                             const ProblemSpecP& restart_prob_spec,
                             GridP& grid)
{
  cout_doing<<"Doing SegeMPM::problemSetup\t\t\t\t\t MPM"<<endl;

  d_mpm->setComponents(this);
  dynamic_cast<ApplicationCommon*>(d_mpm)->problemSetup(prob_spec);

  d_mpm->problemSetup(prob_spec, restart_prob_spec, grid);

  d_8or27 = flags->d_8or27;
  if (d_8or27 == 8) {
      NGN = 1;
  }
  else {
      NGN = 2;
  }

  d_mpm->setParticleGhostLayer(Ghost::AroundNodes, NGP);

  //__________________________________
  //  create analysis modules
  // call problemSetup
  if(!flags->d_with_ice && !flags->d_with_arches){ // mpmice or mpmarches handles this
      d_analysisModules = AnalysisModuleFactory::create(d_myworld,
                                                      m_materialManager,
                                                      prob_spec);

    if(d_analysisModules.size() != 0){
      vector<AnalysisModule*>::iterator iter;
      for( iter  = d_analysisModules.begin();
           iter != d_analysisModules.end(); iter++) {
        AnalysisModule* am = *iter;
        am->setComponents( dynamic_cast<ApplicationInterface*>( this ) );
        am->problemSetup(prob_spec,restart_prob_spec, grid,
            d_mpm->d_mpmd_particleState, d_mpm->d_particleState_preReloc);
      }
    }
  }

  //__________________________________
  //  create the switching criteria port
  d_mpm->d_switchCriteria = dynamic_cast<SwitchingCriteria*>(getPort("switch_criteria"));

  if (d_switchCriteria) {
    d_switchCriteria->problemSetup(restart_mat_ps,
                                   restart_prob_spec, m_materialManager);
  }
}
//______________________________________________________________________
//
void SegeMPM::outputProblemSpec(ProblemSpecP& root_ps)
{

    d_mpm->outputProblemSpec(root_ps);

  //__________________________________
  //  output data analysis modules. Mpmice or mpmarches handles this
  if(!flags->d_with_ice && !flags->d_with_arches && d_analysisModules.size() != 0){

    vector<AnalysisModule*>::iterator iter;
    for( iter  = d_analysisModules.begin();
         iter != d_analysisModules.end(); iter++){
      AnalysisModule* am = *iter;

      am->outputProblemSpec( root_ps );
    }
  }

}

void SegeMPM::scheduleInitialize(const LevelP& level,
                                   SchedulerP& sched)
{
    printSchedule(level, cout_doing, "SegeMPM::scheduleInitialize");

    d_mpm->scheduleInitialize(level, sched);

}
//______________________________________________________________________
//
void SegeMPM::scheduleRestartInitialize(const LevelP& level,
                                          SchedulerP& sched)
{
}

//______________________________________________________________________
void SegeMPM::schedulePrintParticleCount(const LevelP& level,
                                           SchedulerP& sched)
{
    printSchedule(level, cout_doing, "SegeMPM::schedulePrintParticleCount");
    d_mpm->schedulePrintParticleCount(level, sched);

}

//__________________________________
//  Diagnostic task: compute the total number of particles
void SegeMPM::scheduleTotalParticleCount(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
    d_mpm->scheduleTotalParticleCount(level, sched);
}
//__________________________________
//  Diagnostic task: compute the total number of particles

void SegeMPM::scheduleInitializePressureBCs(const LevelP& level,
                                              SchedulerP& sched)
{
    printSchedule(level, cout_doing, "SegeMPM::scheduleInitializePressureBCs");
    d_mpm->scheduleInitializePressureBCs(level, sched);
}

void SegeMPM::scheduleDeleteGeometryObjects(const LevelP& level,
                                              SchedulerP& sched)
{
    printSchedule(level, cout_doing, "SegeMPM::scheduleDeleteGeometryObjects");
    d_mpm->scheduleDeleteGeometryObjects(level, sched);
}

void SegeMPM::scheduleComputeStableTimeStep(const LevelP& level,
                                              SchedulerP& sched)
{
    printSchedule(level, cout_doing, "SegeMPM::scheduleComputeStableTimeStep");
    d_mpm->scheduleComputeStableTimeStep(level, sched);
}

void
SegeMPM::scheduleTimeAdvance(const LevelP & level,
                               SchedulerP   & sched)
{
  if (!flags->doMPMOnLevel(level->getIndex(), level->getGrid()->numLevels()))
    return;

  const PatchSet* patches = level->eachPatch();
  const MaterialSet* matls = m_materialManager->allMaterials( "MPM" );
  const MaterialSet* cz_matls = m_materialManager->allMaterials( "CZ" );
  const MaterialSet* all_matls = m_materialManager->allMaterials();

  const MaterialSubset* mpm_matls_sub = (   matls ?    matls->getUnion() : nullptr);;
  const MaterialSubset*  cz_matls_sub = (cz_matls ? cz_matls->getUnion() : nullptr);

  d_mpm->scheduleComputeCurrentParticleSize(     sched, patches, matls);
  d_mpm->scheduleApplyExternalLoads(             sched, patches, matls);
  if(flags->d_doScalarDiffusion) {
      d_mpm->d_fluxBC->scheduleApplyExternalScalarFlux(sched, patches, matls);
  }
  d_mpm->scheduleInterpolateParticlesToGrid(     sched, patches, matls);
  if(flags->d_computeNormals){
      d_mpm->scheduleComputeNormals(               sched, patches, matls);
  }
  if(flags->d_useLogisticRegression){
      d_mpm->scheduleFindSurfaceParticles(         sched, patches, matls);
      d_mpm->scheduleComputeLogisticRegression(    sched, patches, matls);
  }
  d_mpm->scheduleExMomInterpolated(              sched, patches, matls);
  if(flags->d_doScalarDiffusion) {
      d_mpm->scheduleConcInterpolated(             sched, patches, matls);
  }
  if(flags->d_useCohesiveZones){
      d_mpm->scheduleUpdateCohesiveZones(          sched, patches, mpm_matls_sub,
                                                          cz_matls_sub,
                                                          all_matls);
      d_mpm->scheduleAddCohesiveZoneForces(        sched, patches, mpm_matls_sub,
                                                          cz_matls_sub,
                                                          all_matls);
  }
  if(d_mpm->d_bndy_traction_faces.size()>0) {
      d_mpm->scheduleComputeContactArea(           sched, patches, matls);
  }
  d_mpm->scheduleComputeInternalForce(           sched, patches, matls);
  if (flags->d_doScalarDiffusion) {
      d_mpm->scheduleComputeFlux(                  sched, patches, matls);
      d_mpm->scheduleComputeDivergence(            sched, patches, matls);
      d_mpm->scheduleDiffusionInterfaceDiv(        sched, patches, matls);
  }

  d_mpm->scheduleComputeAndIntegrateAcceleration(sched, patches, matls);
  if (flags->d_doScalarDiffusion) {
      d_mpm->scheduleComputeAndIntegrateDiffusion( sched, patches, matls);
  }
  d_mpm->scheduleExMomIntegrated(                sched, patches, matls);
  d_mpm->scheduleSetGridBoundaryConditions(      sched, patches, matls);
  if (flags->d_prescribeDeformation){
      d_mpm->scheduleSetPrescribedMotion(          sched, patches, matls);
  }
  if(flags->d_XPIC2){
      d_mpm->scheduleComputeSSPlusVp(              sched, patches, matls);
      d_mpm->scheduleComputeSPlusSSPlusVp(         sched, patches, matls);
  }
  if(flags->d_doExplicitHeatConduction){
      d_mpm->scheduleComputeHeatExchange(          sched, patches, matls);
      d_mpm->scheduleComputeInternalHeatRate(      sched, patches, matls);
      d_mpm->scheduleComputeNodalHeatFlux(         sched, patches, matls);
      d_mpm->scheduleSolveHeatEquations(           sched, patches, matls);
      d_mpm->scheduleIntegrateTemperatureRate(     sched, patches, matls);
  }
  d_mpm->scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);
  d_mpm->scheduleComputeParticleGradients(       sched, patches, matls);
  d_mpm->scheduleComputeStressTensor(            sched, patches, matls);
  d_mpm->scheduleFinalParticleUpdate(            sched, patches, matls);
  d_mpm->scheduleInsertParticles(                    sched, patches, matls);
  if(flags->d_computeScaleFactor){
      d_mpm->scheduleComputeParticleScaleFactor(       sched, patches, matls);
  }

  if(d_analysisModules.size() != 0){
    vector<AnalysisModule*>::iterator iter;
    for( iter  = d_analysisModules.begin();
         iter != d_analysisModules.end(); iter++){
      AnalysisModule* am = *iter;
      am->scheduleDoAnalysis_preReloc( sched, level);
    }
  }

  sched->scheduleParticleRelocation(level, Mlb->pXLabel_preReloc,
                                    d_particleState_preReloc,
      Mlb->pXLabel,
                                    d_particleState,
      Mlb->pParticleIDLabel, matls, 1);

 if(flags->d_useCohesiveZones){
  sched->scheduleParticleRelocation(level, Mlb->pXLabel_preReloc,
                                    d_cohesiveZoneState_preReloc,
      Mlb->pXLabel,
                                    d_cohesiveZoneState,
      Mlb->czIDLabel, cz_matls,2);
  }

  //__________________________________
  //  on the fly analysis
  if(d_analysisModules.size() != 0){
    vector<AnalysisModule*>::iterator iter;
    for( iter  = d_analysisModules.begin();
         iter != d_analysisModules.end(); iter++){
      AnalysisModule* am = *iter;
      am->scheduleDoAnalysis( sched, level);
    }
  }
}