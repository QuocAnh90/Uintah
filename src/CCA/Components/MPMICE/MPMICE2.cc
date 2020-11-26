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

// MPMICE2.cc
#include <CCA/Components/MPMICE/MPMICE2.h>
#include <CCA/Components/MPMICE/Core/MPMICELabel.h>

#include <CCA/Components/ICE/AMRICE.h>
#include <CCA/Components/ICE/Core/ICELabel.h>
#include <CCA/Components/ICE/CustomBCs/BoundaryCond.h>
#include <CCA/Components/ICE/EOS/EquationOfState.h>
#include <CCA/Components/ICE/ICE.h>
#include <CCA/Components/MPM/Materials/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/RigidMPM.h>
#include <CCA/Components/MPM/SerialMPM.h>
#include <CCA/Components/MPM/ShellMPM.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>
#include <CCA/Components/OnTheFlyAnalysis/AnalysisModuleFactory.h>

#include <CCA/Ports/Output.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SolverInterface.h>

#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/AMR_CoarsenRefine.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/Utils.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Math/MiscMath.h>

#include <cfloat>
#include <cstdio>
#include <Core/Util/DebugStream.h>

#include <iomanip>
#include <errno.h>
#include <fenv.h>


using namespace Uintah;
using namespace std;
//__________________________________
//  To turn on normal output
//  setenv SCI_DEBUG "MPMICE_NORMAL_COUT:+,MPMICE_DOING_COUT".....
//  MPMICE_NORMAL_COUT:  dumps out during problemSetup 
//  MPMICE_DOING_COUT:   dumps when tasks are scheduled and performed
//  default is OFF


static DebugStream cout_norm("MPMICE2_NORMAL_COUT", false);  
static DebugStream cout_doing("MPMICE2_DOING_COUT", false);
static DebugStream ds_EqPress("DBG_EqPress",false);

MPMICE2::MPMICE2(const ProcessorGroup* myworld,
               const MaterialManagerP materialManager,
               MPMType2 mpmtype, const bool doAMR)
  : ApplicationCommon(myworld, materialManager)
{
  MIlb = scinew MPMICELabel();
 
  d_rigidMPM = false;
  d_testForNegTemps_mpm = true;

  switch(mpmtype) {
      
  case RIGID_MPMICE2:
    d_mpm = scinew RigidMPM(myworld, m_materialManager);
    d_rigidMPM = true;
    break;
  case SHELL_MPMICE2:
    d_mpm = scinew ShellMPM(myworld, m_materialManager);
    break;
    
  default:
    d_mpm = scinew SerialMPM(myworld, m_materialManager);
  }
  
  d_mpmice = scinew MPMICE(myworld, m_materialManager);

  // Don't do AMRICE with MPMICE for now...
  if (doAMR) {
    d_ice  = scinew AMRICE(myworld, m_materialManager);
  }
  else {
    d_ice  = scinew ICE(myworld, m_materialManager);
  }

  d_exchModel = d_ice->d_exchModel;

  Ilb = d_ice->lb;
  Mlb = d_mpm->lb;

  d_SMALL_NUM = d_ice->d_SMALL_NUM;
  d_TINY_RHO  = 1.e-12;  // Note, within MPMICE2, d_TINY_RHO is only applied
                         // to MPM materials, for which its value is hardcoded,
                         // unlike the situation for ice materials

  d_switchCriteria = 0;
}
//______________________________________________________________________
//
MPMICE2::~MPMICE2()
{
  d_mpm->releaseComponents();
  d_ice->releaseComponents();

  delete MIlb;
  delete d_mpm;
  delete d_ice;
  
  if(d_analysisModules.size() != 0){
    vector<AnalysisModule*>::iterator iter;
    for( iter  = d_analysisModules.begin();
         iter != d_analysisModules.end(); iter++){
      AnalysisModule* am = *iter;
      am->releaseComponents();
      delete am;
    }
  }
}

//__________________________________
//    For recomputing timesteps
double MPMICE2::recomputeDelT(const double delT)
{
  return delT / 2.0;
} 
//______________________________________________________________________
//
void MPMICE2::problemSetup(const ProblemSpecP& prob_spec, 
                          const ProblemSpecP& restart_prob_spec, 
                          GridP& grid)
{
  cout_doing << "Doing MPMICE2::problemSetup " << endl;

  //__________________________________
  //  M P M
  d_mpm->setComponents( this );
  dynamic_cast<ApplicationCommon*>(d_mpm)->problemSetup( prob_spec );

  d_mpm->setWithICE();
  d_mpm->problemSetup(prob_spec, restart_prob_spec,grid);

  d_8or27 = d_mpm->flags->d_8or27; 
  if(d_8or27==8){
    NGN=1;
  } else{
    NGN=2;
  }

  Ghost::GhostType gp;
  int ngc_p;
  d_mpm->getParticleGhostLayer(gp, ngc_p);

  //__________________________________
  //  I C E
  d_ice->setComponents( this );
  dynamic_cast<ApplicationCommon*>(d_ice)->problemSetup( prob_spec );

  d_ice->setWithMPM();

  // Communicate the particle ghost from MPM to ICE. Used only in the
  // HEChem/Unsteady_Burn model.
  d_ice->setParticleGhostLayer(gp, ngc_p);

  if(d_rigidMPM){
   d_ice->setWithRigidMPM();
  }

  d_switchCriteria = dynamic_cast<SwitchingCriteria*>
    (getPort("switch_criteria"));
  
  if (d_switchCriteria) {
    d_switchCriteria->problemSetup(prob_spec,restart_prob_spec,m_materialManager);
  }

  d_ice->problemSetup(prob_spec, restart_prob_spec,grid);

  ProblemSpecP mpm_ps = 0;
  mpm_ps = prob_spec->findBlock("MPM");
  
  if(!mpm_ps){
    mpm_ps = restart_prob_spec->findBlock("MPM");
  }
  mpm_ps->get("testForNegTemps_mpm",d_testForNegTemps_mpm);
  
  //__________________________________
  //  bulletproofing
  if(isAMR() && !isLockstepAMR()){
    ostringstream msg;
    msg << "\n ERROR: You must add \n"
        << " <useLockStep> true </useLockStep> \n"
        << " inside of the <AMR> section for MPMICE2 and AMR. \n"; 
    throw ProblemSetupException(msg.str(),__FILE__, __LINE__);
  }
    
  if (cout_norm.active()) {
    cout_norm << "Done with problemSetup \t\t\t MPMICE2" <<endl;
    cout_norm << "--------------------------------\n"<<endl;
  }
  
  //__________________________________
  //  Set up data analysis modules
  d_analysisModules = AnalysisModuleFactory::create(d_myworld,
                                                    m_materialManager,
                                                    prob_spec);

  if(d_analysisModules.size() != 0){
    vector<AnalysisModule*>::iterator iter;
    for( iter  = d_analysisModules.begin();
         iter != d_analysisModules.end(); iter++){
      AnalysisModule* am = *iter;
      am->setComponents( dynamic_cast<ApplicationInterface*>( this ) );
      am->problemSetup(prob_spec, restart_prob_spec, grid,
                       d_mpm->d_particleState, d_mpm->d_particleState_preReloc);
    }
  }  
}

//______________________________________________________________________
//
void MPMICE2::outputProblemSpec(ProblemSpecP& root_ps)
{
  d_mpm->outputProblemSpec(root_ps);
  d_ice->outputProblemSpec(root_ps);
  
  // Global flags required by MPMICE2
  ProblemSpecP mpm_ps = root_ps->findBlock("MPM");
  mpm_ps->appendElement("testForNegTemps_mpm", d_testForNegTemps_mpm);  
  
  //__________________________________
  //  output data analysis modules
  if( d_analysisModules.size() != 0 ){

    vector<AnalysisModule*>::iterator iter;
    for( iter  = d_analysisModules.begin();
         iter != d_analysisModules.end(); iter++){
      AnalysisModule* am = *iter;

      am->outputProblemSpec( root_ps );
    }
  } 
}

//______________________________________________________________________
//
void MPMICE2::scheduleInitialize(const LevelP& level,
                            SchedulerP& sched)
{
  printSchedule(level,cout_doing,"MPMICE2::scheduleInitialize");

  d_mpm->scheduleInitialize(level, sched);
  d_ice->scheduleInitialize(level, sched);

  //__________________________________
  //  What isn't initialized in either ice or mpm
  Task* t = scinew Task("MPMICE2::actuallyInitialize",
                  this, &MPMICE2::actuallyInitialize);
                  
  // Get the material subsets
  const MaterialSubset* ice_matls = m_materialManager->allMaterials( "ICE" )->getUnion();
  const MaterialSubset* mpm_matls = m_materialManager->allMaterials( "MPM" )->getUnion();

  t->requires(Task::NewDW, Ilb->timeStepLabel);
  // These values are calculated for ICE materials in d_ice->actuallyInitialize(...)
  //  so they are only needed for MPM
  t->computes(MIlb->vel_CCLabel,       mpm_matls);
  t->computes(Ilb->rho_CCLabel,        mpm_matls); 
  t->computes(Ilb->temp_CCLabel,       mpm_matls);
  t->computes(Ilb->sp_vol_CCLabel,     mpm_matls);
  t->computes(Ilb->speedSound_CCLabel, mpm_matls); 
  t->computes(Mlb->heatRate_CCLabel,   mpm_matls);

  // This is compute in d_ice->actuallyInitalize(...), and it is needed in 
  //  MPMICE2's actuallyInitialize()
  t->requires(Task::NewDW, Ilb->vol_frac_CCLabel, ice_matls, Ghost::None, 0);

  if (d_switchCriteria) {
    d_switchCriteria->scheduleInitialize(level,sched);
  }
  
  //__________________________________
  // dataAnalysis 
  if(d_analysisModules.size() != 0){
    vector<AnalysisModule*>::iterator iter;
    for( iter  = d_analysisModules.begin();
         iter != d_analysisModules.end(); iter++){
      AnalysisModule* am = *iter;
      am->scheduleInitialize( sched, level);
    }
  }
    
  sched->addTask(t, level->eachPatch(), m_materialManager->allMaterials());
}

//______________________________________________________________________
//       A C T U A L   S T E P S :
//______________________________________________________________________
void MPMICE2::actuallyInitialize(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse*,
    DataWarehouse* new_dw)
{
    timeStep_vartype timeStep;
    new_dw->get(timeStep, VarLabel::find(timeStep_name));

    bool isNotInitialTimeStep = (timeStep > 0);

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing actuallyInitialize ");
        //__________________________________
        //output material indices
        if (patch->getID() == 0) {
            cout << "Materials Indicies:   MPM [" << *(m_materialManager->allMaterials("MPM")) << "] "
                << "ICE[" << *(m_materialManager->allMaterials("ICE")) << "]" << endl;

            cout << "Material Names:";
            unsigned int numAllMatls = m_materialManager->getNumMatls();
            for (unsigned int m = 0; m < numAllMatls; m++) {
                Material* matl = m_materialManager->getMaterial(m);
                cout << " " << matl->getDWIndex() << ") " << matl->getName();
            }
            cout << "\n";
        }

        // Sum variable for testing that the volume fractions sum to 1
        CCVariable<double> vol_frac_sum;
        new_dw->allocateTemporary(vol_frac_sum, patch);
        vol_frac_sum.initialize(0.0);

        //__________________________________
        //  Initialize CCVaribles for MPM Materials
        //  Even if mass = 0 in a cell you still need
        //  CC Variables defined.
        double junk = -9, tmp;
        unsigned int numMPM_matls = m_materialManager->getNumMatls("MPM");
        double p_ref = d_ice->getRefPress();
        for (unsigned int m = 0; m < numMPM_matls; m++) {
            CCVariable<double> rho_micro, sp_vol_CC, rho_CC, Temp_CC, speedSound, vol_frac_CC;
            CCVariable<Vector> vel_CC;
            MPMMaterial* mpm_matl = (MPMMaterial*)m_materialManager->getMaterial("MPM", m);
            int indx = mpm_matl->getDWIndex();
            new_dw->allocateTemporary(rho_micro, patch);
            // Allocate volume fraction for use in intializeCCVariables
            new_dw->allocateTemporary(vol_frac_CC, patch);
            new_dw->allocateAndPut(sp_vol_CC, Ilb->sp_vol_CCLabel, indx, patch);
            new_dw->allocateAndPut(rho_CC, Ilb->rho_CCLabel, indx, patch);
            new_dw->allocateAndPut(speedSound, Ilb->speedSound_CCLabel, indx, patch);
            new_dw->allocateAndPut(Temp_CC, MIlb->temp_CCLabel, indx, patch);
            new_dw->allocateAndPut(vel_CC, MIlb->vel_CCLabel, indx, patch);


            CCVariable<double> heatFlux;
            new_dw->allocateAndPut(heatFlux, Mlb->heatRate_CCLabel, indx, patch);
            heatFlux.initialize(0.0);

            mpm_matl->initializeCCVariables(rho_micro, rho_CC,
                Temp_CC, vel_CC,
                vol_frac_CC, patch);

            setBC(rho_CC, "Density", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            setBC(rho_micro, "Density", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            setBC(Temp_CC, "Temperature", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            setBC(vel_CC, "Velocity", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            for (CellIterator iter = patch->getExtraCellIterator();
                !iter.done(); iter++) {
                IntVector c = *iter;
                sp_vol_CC[c] = 1.0 / rho_micro[c];

                mpm_matl->getConstitutiveModel()->
                    computePressEOSCM(rho_micro[c], junk, p_ref, junk, tmp, mpm_matl, Temp_CC[c]);
                speedSound[c] = sqrt(tmp);

                // sum volume fraction
                vol_frac_sum[c] += vol_frac_CC[c];
            }

            //__________________________________
            //    B U L L E T   P R O O F I N G
            IntVector neg_cell;
            ostringstream warn;
            if (!areAllValuesPositive(rho_CC, neg_cell)) {
                Point pt = patch->getCellPosition(neg_cell);
                warn << "ERROR MPMICE2::actuallyInitialize, mat " << indx << " cell "
                    << neg_cell << " position " << pt << " rho_CC is negative\n";
                throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
            }
            if (!areAllValuesPositive(Temp_CC, neg_cell)) {
                Point pt = patch->getCellPosition(neg_cell);

                warn << "ERROR MPMICE2::actuallyInitialize, mat " << indx << " cell "
                    << neg_cell << " position " << pt << " Temp_CC is negative\n";
                throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
            }
            if (!areAllValuesPositive(sp_vol_CC, neg_cell)) {
                Point pt = patch->getCellPosition(neg_cell);

                warn << "ERROR MPMICE2::actuallyInitialize, mat " << indx << " cell "
                    << neg_cell << " position " << pt << " sp_vol_CC is negative\n";
                throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
            }
            if (!areAllValuesNumbers(speedSound, neg_cell)) {
                Point pt = patch->getCellPosition(neg_cell);

                warn << "ERROR MPMICE2::actuallyInitialize, mat " << indx << " cell "
                    << neg_cell << " position " << pt << " speedSound is nan\n";
                warn << "speedSound = " << speedSound[neg_cell] << " sp_vol_CC = " << sp_vol_CC[neg_cell]
                    << " rho_micro = " << rho_micro[neg_cell] << " Temp_CC = " << Temp_CC[neg_cell] << endl;
                throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
            }
        }  // num_MPM_matls loop 

        //___________________________________
        //   B U L L E T  P R O O F I N G
        // Verify volume fractions sum to 1.0
        // Loop through ICE materials to get their contribution to volume fraction
        unsigned int numICE_matls = m_materialManager->getNumMatls("ICE");
        for (unsigned int m = 0; m < numICE_matls; m++) {
            constCCVariable<double> vol_frac;
            ICEMaterial* ice_matl = (ICEMaterial*)m_materialManager->getMaterial("ICE", m);
            int indx = ice_matl->getDWIndex();

            // Get the Volume Fraction computed in ICE's actuallyInitialize(...)
            new_dw->get(vol_frac, Ilb->vol_frac_CCLabel, indx, patch, Ghost::None, 0);

            for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                vol_frac_sum[c] += vol_frac[c];
            }
        }  // num_ICE_matls loop

        double errorThresholdTop = 1.0e0 + 1.0e-10;
        double errorThresholdBottom = 1.0e0 - 1.0e-10;

        for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
            IntVector c = *iter;
            Point pt = patch->getCellPosition(c);

            if (!(vol_frac_sum[c] <= errorThresholdTop && vol_frac_sum[c] >= errorThresholdBottom)) {
                \
                    ostringstream warn;
                warn << "ERROR MPMICE2::actuallyInitialize cell " << *iter << " position" << pt << "\n\n"
                    << "volume fraction (" << std::setprecision(13) << vol_frac_sum[*iter] << ") does not sum to 1.0 +- 1e-10.\n"
                    << "Verify that this region of the domain contains at least 1 geometry object.  If you're using the optional\n"
                    << "'volumeFraction' tags verify that they're correctly specified.\n";
                throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
            }
        } // cell iterator for volume fraction
    } // Patch loop
}

//______________________________________________________________________
//
void MPMICE2::scheduleRestartInitialize(const LevelP& level,
                                       SchedulerP& sched)
{
  printSchedule(level,cout_doing,"MPMICE2::scheduleInitialize");

  d_mpm->scheduleRestartInitialize(level, sched);
  d_ice->scheduleRestartInitialize(level, sched);
  
  //__________________________________
  // dataAnalysis 
  if(d_analysisModules.size() != 0){
    vector<AnalysisModule*>::iterator iter;
    for( iter  = d_analysisModules.begin();
         iter != d_analysisModules.end(); iter++){
      AnalysisModule* am = *iter;
      am->scheduleRestartInitialize( sched, level);
    }
  }
}

//______________________________________________________________________
//
void MPMICE2::restartInitialize()
{
  if (cout_doing.active())
    cout_doing <<"Doing restartInitialize \t\t\t MPMICE2" << endl;

  d_mpm->restartInitialize();
  d_ice->restartInitialize();
  
  if(d_analysisModules.size() != 0){
    vector<AnalysisModule*>::iterator iter;
    for( iter  = d_analysisModules.begin();
         iter != d_analysisModules.end(); iter++){
      AnalysisModule* am = *iter;
      am->restartInitialize();
    }
  }
}

//______________________________________________________________________
//
void MPMICE2::scheduleComputeStableTimeStep(const LevelP& level,
                                      SchedulerP& sched)
{
  // Schedule computing the ICE stable timestep
  d_ice->scheduleComputeStableTimeStep(level, sched);
  // MPM stable timestep is a by product of the CM
}

//______________________________________________________________________
//
void
MPMICE2::scheduleTimeAdvance(const LevelP& inlevel, SchedulerP& sched)
{
  // Only do scheduling on level 0 for lockstep AMR
  if(inlevel->getIndex() > 0 && isLockstepAMR())
    return;

  // If we have a finer level, then assume that we are doing multilevel MPMICE2
  // Otherwise, it is plain-ole MPMICE2
  do_mlmpmice2 = false;
  if(inlevel->hasFinerLevel()){
    do_mlmpmice2 = true;
  }
  const LevelP& mpm_level = do_mlmpmice2? inlevel->getGrid()->getLevel(inlevel->getGrid()->numLevels()-1) : inlevel;

  const PatchSet* mpm_patches = mpm_level->eachPatch();
  const MaterialSet* ice_matls = m_materialManager->allMaterials( "ICE" );
  const MaterialSet* mpm_matls = m_materialManager->allMaterials( "MPM" );
  const MaterialSet* all_matls = m_materialManager->allMaterials();
  MaterialSubset* press_matl   = d_ice->d_press_matl;
  MaterialSubset* one_matl     = d_ice->d_press_matl;

  const MaterialSubset* ice_matls_sub = ice_matls->getUnion();
  const MaterialSubset* mpm_matls_sub = mpm_matls->getUnion();
  cout_doing << "---------------------------------------------------------Level ";
  if(do_mlmpmice2){
    cout_doing << inlevel->getIndex() << " (ICE) " << mpm_level->getIndex() << " (MPM)"<< endl;;
  } else {
    cout_doing << inlevel->getIndex()<< endl;
  }

 //__________________________________
 // Scheduling
  for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
    const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
    d_ice->scheduleComputeThermoTransportProperties(sched, ice_level,ice_matls);
    
    d_ice->scheduleMaxMach_on_Lodi_BC_Faces(        sched, ice_level,ice_matls);
  }
  
  // diagnostic task
  //d_mpm->scheduleTotalParticleCount(          sched, mpm_patches, mpm_matls);
   
  d_mpm->scheduleApplyExternalLoads(          sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeCurrentParticleSize(  sched, mpm_patches, mpm_matls);
  d_mpm->scheduleInterpolateParticlesToGrid(  sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeHeatExchange(         sched, mpm_patches, mpm_matls);

  if(d_mpm->flags->d_computeNormals){
    d_mpm->scheduleComputeNormals(            sched, mpm_patches, mpm_matls);
  }
  d_mpm->scheduleExMomInterpolated(           sched, mpm_patches, mpm_matls);

  // schedule the interpolation of mass and volume to the cell centers
  d_mpmice->scheduleInterpolateNCToCC_0(                sched, mpm_patches, one_matl,
                                                                  mpm_matls);
                                                                  
  // do coarsens in reverse order, and before the other tasks
  if(do_mlmpmice2){
    for (int l = inlevel->getGrid()->numLevels() - 2; l >= 0; l--) {
      const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
      const PatchSet* ice_patches = ice_level->eachPatch();

      scheduleCoarsenCC_0(                      sched, ice_patches, mpm_matls);
      scheduleCoarsenNCMass(                    sched, ice_patches, mpm_matls);
    }
  }

  for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
    const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
    const PatchSet* ice_patches = ice_level->eachPatch();
    
    d_mpmice->scheduleComputePressure(                  sched, ice_patches, ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  press_matl,
                                                                  all_matls);
    
  
    d_ice->scheduleComputeTempFC(             sched, ice_patches, ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  all_matls);
    d_ice->scheduleComputeModelSources(       sched, ice_level,   all_matls);

    d_ice->scheduleUpdateVolumeFraction(      sched, ice_level,   press_matl,
                                                                  all_matls);
  
    d_ice->scheduleComputeVel_FC(             sched, ice_patches, ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  press_matl, 
                                                                  all_matls);
                                                                  
    d_ice->d_exchModel->sched_PreExchangeTasks( sched, ice_patches, ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  all_matls);

    d_ice->d_exchModel->sched_AddExch_VelFC(  sched, ice_patches, ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  all_matls,
                                                                  d_ice->d_BC_globalVars,
                                                                  false);
 
  }
  if(d_ice->d_impICE) {        //  I M P L I C I T, won't work with AMR yet
    // we should use the AMR multi-level pressure solve
    for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
      const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
      const PatchSet* ice_patches = ice_level->eachPatch();

      d_ice->scheduleSetupRHS(                sched, ice_patches, one_matl, 
                                                                  all_matls,
                                                                  false,
                                                                  "computes");
      d_ice->scheduleCompute_maxRHS(          sched, ice_level,    one_matl,
                                                                   all_matls);
                                                                  
      d_ice->scheduleImplicitPressureSolve(   sched, ice_level,   ice_patches,
                                                                  one_matl, 
                                                                  press_matl,
                                                                  ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  all_matls);
                                                           
      d_ice->scheduleComputeDel_P(            sched, ice_level,   ice_patches, 
                                                                  one_matl, 
                                                                  press_matl,
                                                                  all_matls);
    }
  }                           //  IMPLICIT AND EXPLICIT

                                                                  
  if(!(d_ice->d_impICE)){       //  E X P L I C I T 
    for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
      const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
      const PatchSet* ice_patches = ice_level->eachPatch();

      d_ice->scheduleComputeDelPressAndUpdatePressCC(
                                              sched, ice_patches, press_matl,
                                                                  ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  all_matls);
    }
  } 
  
  for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
    const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
    const PatchSet* ice_patches = ice_level->eachPatch();

    d_ice->scheduleComputePressFC(            sched, ice_patches, press_matl,
                                                                    all_matls);
                                                                    
    d_ice->scheduleVelTau_CC(                 sched, ice_patches, ice_matls);
    
    d_ice->scheduleViscousShearStress(        sched, ice_patches, ice_matls);
   
    d_ice->scheduleAccumulateMomentumSourceSinks(
                                              sched, ice_patches, press_matl,
                                                                  ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  all_matls);
                                                                  
    d_ice->scheduleAccumulateEnergySourceSinks(sched, ice_patches,ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  press_matl,
                                                                  all_matls);
  }

  if(!d_rigidMPM){
//    if(do_mlmpmice2){
//      scheduleRefinePressCC(                  sched, mpm_patches, press_matl,
//                                                                  mpm_matls);
//    }
      
      d_mpmice->scheduleInterpolatePressCCToPressNC(      sched, mpm_patches, press_matl,
                                                                  mpm_matls);

      d_mpmice->scheduleInterpolatePAndGradP(             sched, mpm_patches, press_matl,
                                                                  one_matl,
                                                                  mpm_matls_sub,
                                                                  mpm_matls);
  }
   
  d_mpm->scheduleComputeInternalForce(        sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeInternalHeatRate(     sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeNodalHeatFlux(        sched, mpm_patches, mpm_matls);
  d_mpm->scheduleSolveHeatEquations(          sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeAndIntegrateAcceleration(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleIntegrateTemperatureRate(    sched, mpm_patches, mpm_matls);
  
  d_mpmice->scheduleComputeLagrangianValuesMPM(         sched, mpm_patches, one_matl,
                                                                  mpm_matls); 

  // do coarsens in reverse order, and before the other tasks
  if(do_mlmpmice2){
    for (int l = inlevel->getGrid()->numLevels() - 2; l >= 0; l--) {
      const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
      const PatchSet* ice_patches = ice_level->eachPatch();
      scheduleCoarsenLagrangianValuesMPM(     sched, ice_patches, mpm_matls);
    }
  }

  for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
    const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
    const PatchSet* ice_patches = ice_level->eachPatch();

    d_ice->scheduleComputeLagrangianValues(   sched, ice_patches, ice_matls);
                                                                  
    d_ice->d_exchModel->sched_AddExch_Vel_Temp_CC(   
                                              sched, ice_patches, ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  all_matls,
                                                                  d_ice->d_BC_globalVars); 

    d_ice->scheduleComputeLagrangianSpecificVolume(
                                              sched, ice_patches, ice_matls_sub,
                                                                  mpm_matls_sub,
                                                                  press_matl,
                                                                  all_matls);
                                                                  
    d_ice->scheduleComputeLagrangian_Transported_Vars(
                                              sched, ice_patches, ice_matls);

  }

  d_mpmice->scheduleComputeCCVelAndTempRates(           sched, mpm_patches, mpm_matls);

//  if(do_mlmpmice2){
//    scheduleRefineCC(                         sched, mpm_patches, mpm_matls);
//  }

  d_mpmice->scheduleInterpolateCCToNC(                  sched, mpm_patches, mpm_matls);

  d_mpm->scheduleExMomIntegrated(             sched, mpm_patches, mpm_matls);
  d_mpm->scheduleSetGridBoundaryConditions(   sched, mpm_patches, mpm_matls);
  d_mpm->scheduleInterpolateToParticlesAndUpdate(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeParticleGradients(       sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeStressTensor(         sched, mpm_patches, mpm_matls);
  d_mpm->scheduleFinalParticleUpdate(         sched, mpm_patches, mpm_matls);
  if( d_mpm->flags->d_computeScaleFactor ){
    d_mpm->scheduleComputeParticleScaleFactor(sched, mpm_patches, mpm_matls);
  }
  
  for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
    const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
    const PatchSet* ice_patches = ice_level->eachPatch();
                                   
    d_ice->scheduleAdvectAndAdvanceInTime(   sched, ice_patches,ice_matls_sub,
                                                                ice_matls);
                                                                
    d_ice->scheduleConservedtoPrimitive_Vars(sched, ice_patches,ice_matls_sub,
                                                    ice_matls,"afterAdvection");
  }
} // end scheduleTimeAdvance()

// Optional function //
//______________________________________________________________________
//
void MPMICE2::scheduleRefinePressCC(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSubset* press_matl,
                                   const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"MPMICE2::scheduleRefinePressCC");
    
  MaterialSet* press_matls = scinew MaterialSet();
  press_matls->add(0);
  press_matls->addReference();

  scheduleRefineVariableCC(sched,patches, press_matls,Ilb->press_CCLabel);
  if(press_matls->removeReference())
    delete press_matls;
}

//if (do_mlmpmice2)//

//______________________________________________________________________
//
void MPMICE2::scheduleCoarsenLagrangianValuesMPM(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* mpm_matls)
{
  printSchedule(patches,cout_doing,"MPMICE2:scheduleCoarsenLagrangianValues mpm_matls");

  scheduleCoarsenVariableCC(sched, patches, mpm_matls, Ilb->rho_CCLabel,
                            1e-12,          true, "std"); // modifies
  scheduleCoarsenVariableCC(sched, patches, mpm_matls, Ilb->mass_L_CCLabel,
                            1.9531e-15,     false, "sum");
  scheduleCoarsenVariableCC(sched, patches, mpm_matls, Ilb->mom_L_CCLabel,
                            Vector(0, 0, 0),false, "sum");
  scheduleCoarsenVariableCC(sched, patches, mpm_matls, Ilb->int_eng_L_CCLabel,
                             0.0,           false, "sum");
}

//______________________________________________________________________
//
void MPMICE2::scheduleRefineCC(SchedulerP& sched,
                              const PatchSet* patches,
                              const MaterialSet* mpm_matls)
{
  if(!d_mpm->flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing, "MPMICE2::scheduleRefineCC");
  scheduleRefineVariableCC(sched, patches, mpm_matls, Ilb->dTdt_CCLabel);
  scheduleRefineVariableCC(sched, patches, mpm_matls, Ilb->dVdt_CCLabel);
}

void MPMICE2::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  if (d_switchCriteria) {
    d_switchCriteria->scheduleSwitchTest(level,sched);
  }
}

//______________________________________________________________________
void MPMICE2::scheduleRefineInterface(const LevelP& fineLevel,
                                     SchedulerP& scheduler,
                                     bool needOld, bool needNew)
{
  d_ice->scheduleRefineInterface(fineLevel, scheduler, needOld, needNew);
  d_mpm->scheduleRefineInterface(fineLevel, scheduler, needOld, needNew);

  if(fineLevel->getIndex() > 0 && scheduler->isCopyDataTimestep() &&
     d_mpm->flags->doMPMOnLevel(fineLevel->getIndex(),
                                fineLevel->getGrid()->numLevels())) {
    cout_doing << d_myworld->myRank() 
               << " MPMICE2::scheduleRefineInterface \t\t\tL-"
               << fineLevel->getIndex() << endl;

    Task* task = scinew Task("MPMICE2::refineCoarseFineInterface",
                             this, &MPMICE2::refineCoarseFineInterface);

    const MaterialSet* all_matls   = m_materialManager->allMaterials();
    const MaterialSubset* one_matl = d_ice->d_press_matl;

    task->modifies(Mlb->NC_CCweightLabel, one_matl);

    scheduler->addTask(task, fineLevel->eachPatch(), all_matls);
  }
}


void MPMICE2::refineCoarseFineInterface(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset*,
                                       DataWarehouse* fine_old_dw,
                                       DataWarehouse* fine_new_dw)
{
  // This isn't actually refining anything, it is simply reinitializing
  // NC_CCweight after regridding on all levels finer than 0 because
  // copyData doesn't copy extra cell data.
  const Level* level = getLevel(patches);
  if(level->getIndex() > 0){
    cout_doing << d_myworld->myRank()
               << " Doing refineCoarseFineInterface"<< "\t\t\t MPMICE2 L-"
               << level->getIndex() << " Patches: " << *patches << endl;

    for(int p=0;p<patches->size();p++){
      const Patch* patch = patches->get(p);
      //__________________________________
      //NC_CCweight
      NCVariable<double> NC_CCweight;
      fine_new_dw->getModifiable(NC_CCweight, Mlb->NC_CCweightLabel, 0, patch);
      //__________________________________
      // - Initialize NC_CCweight = 0.125
      // - Find the walls with symmetry BC and double NC_CCweight
      NC_CCweight.initialize(0.125);
      vector<Patch::FaceType>::const_iterator iter;
      vector<Patch::FaceType> bf;
      patch->getBoundaryFaces(bf);
      
      for (iter  = bf.begin(); iter != bf.end(); ++iter){
        Patch::FaceType face = *iter;
        int mat_id = 0;
        if (patch->haveBC(face,mat_id,"symmetry","Symmetric")) {
             
          for(CellIterator iter = patch->getFaceIterator(face,Patch::FaceNodes);
              !iter.done(); iter++) {
            NC_CCweight[*iter] = 2.0*NC_CCweight[*iter];
          } // cell iterator
        } // if symmetry
      } // for patch faces
    } // for patches
  } // if level
}
//______________________________________________________________________
void MPMICE2::scheduleRefine(const PatchSet* patches, 
                            SchedulerP& sched)
{
  d_ice->scheduleRefine(patches, sched);
  d_mpm->scheduleRefine(patches, sched);

  printSchedule(patches,cout_doing,"MPMICE2::scheduleRefine");

  Task* task = scinew Task("MPMICE2::refine", this, &MPMICE2::refine);
  
  task->requires(Task::OldDW, Ilb->timeStepLabel);
  
  task->computes(Mlb->heatRate_CCLabel);
  task->computes(Ilb->sp_vol_CCLabel);
  task->computes(MIlb->vel_CCLabel);
  task->computes(Ilb->temp_CCLabel);

  sched->addTask(task, patches, m_materialManager->allMaterials( "MPM" ));
}

//______________________________________________________________________    
void MPMICE2::scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& sched)
{
  d_ice->scheduleCoarsen(coarseLevel, sched);
  d_mpm->scheduleCoarsen(coarseLevel, sched);
}

//______________________________________________________________________
void MPMICE2::scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                                          SchedulerP& sched)
{
  d_ice->scheduleInitialErrorEstimate(coarseLevel, sched);
  d_mpm->scheduleInitialErrorEstimate(coarseLevel, sched);
}

//______________________________________________________________________
void MPMICE2::scheduleErrorEstimate(const LevelP& coarseLevel,
                                   SchedulerP& sched)
{
  d_ice->scheduleErrorEstimate(coarseLevel, sched);
  d_mpm->scheduleErrorEstimate(coarseLevel, sched);
}

//______________________________________________________________________
 void MPMICE2::scheduleRefineVariableCC(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls,
                                       const VarLabel* variable)
 {
   ostringstream taskName;
   taskName << "MPMICE2::refineVariable(" << variable->getName() << ")";
   Task* t;

   // the sgis don't like accepting a templated function over a function call for some reason...
   void (MPMICE2::*func)(const ProcessorGroup*, const PatchSubset*, const MaterialSubset*,
                        DataWarehouse*, DataWarehouse*, const VarLabel*);
                        
   switch(variable->typeDescription()->getSubType()->getType()){
   case TypeDescription::double_type:
     func = &MPMICE2::refineVariableCC<double>;
     t=scinew Task(taskName.str().c_str(),this, func, variable);
     break;
   case TypeDescription::Vector:
     func = &MPMICE2::refineVariableCC<Vector>;
     t=scinew Task(taskName.str().c_str(),this, func, variable);
     break;
   default:
     throw InternalError("Unknown variable type for refine", __FILE__, __LINE__);
   }

   Ghost::GhostType  gac = Ghost::AroundCells;
   t->requires(Task::NewDW, variable, 0, Task::CoarseLevel, 0, Task::NormalDomain, gac, 1);
   t->computes(variable);
   sched->addTask(t, patches, matls);
 }

 //______________________________________________________________________
 template<typename T>
   void MPMICE2::scheduleCoarsenVariableCC(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls,
                                          const VarLabel* variable,
                                          T defaultValue, 
                                          bool modifies,
                                          const string& coarsenMethod)
{
  // The SGI compiler does't like accepting a templated function over
  // a function call for some reason...  We use this hack to force it
  // to figure out the correct type of the function.
  void (MPMICE2::*func)(const ProcessorGroup*, const PatchSubset*,
                       const MaterialSubset*, DataWarehouse*, DataWarehouse*, 
                       const VarLabel*, T, bool, string);
  func = &MPMICE2::coarsenVariableCC<T>;
  ostringstream taskName;

  taskName << "MPMICE2::coarsenVariableCC(" << variable->getName() 
           << (modifies?" modified":"") << ")";
  
  Task* t=scinew Task(taskName.str().c_str(),this, func, 
                       variable, defaultValue, modifies, coarsenMethod);
  
  Ghost::GhostType  gn = Ghost::None;
  Task::MaterialDomainSpec ND   = Task::NormalDomain;
  
  t->requires(Task::OldDW, Ilb->timeStepLabel);

  t->requires(Task::NewDW, variable, 0, Task::FineLevel, 0, ND,gn,0);
  
  if(coarsenMethod == "massWeighted"){
    t->requires(Task::NewDW, MIlb->cMassLabel, 0, Task::FineLevel, 0, ND,gn,0);
  }
  
  if(modifies){
    t->modifies(variable);
  }else{
    t->computes(variable);
  }
  sched->addTask(t, patches, matls);
}
 

 //______________________________________________________________________
 template<typename T>
   void MPMICE2::scheduleCoarsenVariableNC(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls,
                                          const VarLabel* variable,
                                          T defaultValue,
                                          bool modifies,
                                          string coarsenMethod)
{
  // The SGI compiler does't like accepting a templated function over
  // a function call for some reason...  We use this hack to force it
  // to figure out the correct type of the function.
  void (MPMICE2::*func)(const ProcessorGroup*, const PatchSubset*,
                       const MaterialSubset*, DataWarehouse*, DataWarehouse*,
                       const VarLabel*, T, bool, string);
  func = &MPMICE2::coarsenVariableNC<T>;
  ostringstream taskName;

  taskName << "MPMICE2::coarsenVariableNC(" << variable->getName() 
           << (modifies?" modified":"") << ")";

  Task* t=scinew Task(taskName.str().c_str(),this, func,
                       variable, defaultValue, modifies, coarsenMethod);

  //Ghost::GhostType  gn = Ghost::None;
  Ghost::GhostType  gan = Ghost::AroundNodes;
  Task::MaterialDomainSpec ND   = Task::NormalDomain;

  const LevelP fineLevel = getLevel(patches)->getFinerLevel();
  IntVector refineRatio(fineLevel->getRefinementRatio());
  int ghost = max(refineRatio.x(),refineRatio.y());
  ghost = max(ghost,refineRatio.z());

  t->requires(Task::NewDW, variable, 0, Task::FineLevel, 0, ND, gan, ghost);

  if(modifies){
    t->modifies(variable);
  }else{
    t->computes(variable);
  }
  sched->addTask(t, patches, matls);
}

//______________________________________________________________________
void
MPMICE2::refine(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* /*matls*/,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, VarLabel::find( timeStep_name) );

  bool isNotInitialTimeStep = (timeStep > 0);

  for (int p = 0; p<patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,patch,cout_doing,"Doing refine");
     
    unsigned int numMPMMatls=m_materialManager->getNumMatls( "MPM" );
    
    for(unsigned int m = 0; m < numMPMMatls; m++){
      MPMMaterial* mpm_matl = (MPMMaterial*) m_materialManager->getMaterial( "MPM",  m );
      int dwi = mpm_matl->getDWIndex();

      cout_doing << d_myworld->myRank() << " Doing refine on patch "
           << patch->getID() << " material # = " << dwi << endl;
      
      // for now, create 0 heat flux
      CCVariable<double> heatFlux;
      new_dw->allocateAndPut(heatFlux, Mlb->heatRate_CCLabel, dwi, patch);
      heatFlux.initialize(0.0);

      CCVariable<double> rho_micro, sp_vol_CC, rho_CC, Temp_CC, vol_frac_CC;
      CCVariable<Vector> vel_CC;

      new_dw->allocateTemporary(rho_micro,   patch);
      new_dw->allocateTemporary(rho_CC,      patch);
      new_dw->allocateTemporary(vol_frac_CC, patch);

      new_dw->allocateAndPut(sp_vol_CC,   Ilb->sp_vol_CCLabel,    dwi,patch);
      new_dw->allocateAndPut(Temp_CC,     MIlb->temp_CCLabel,     dwi,patch);
      new_dw->allocateAndPut(vel_CC,      MIlb->vel_CCLabel,      dwi,patch);

      mpm_matl->initializeDummyCCVariables(rho_micro,   rho_CC,
                                           Temp_CC,     vel_CC,
                                           vol_frac_CC, patch);  
      //__________________________________
      //  Set boundary conditions                                     
      setBC(rho_micro, "Density",      patch, m_materialManager, dwi, new_dw, isNotInitialTimeStep);
      setBC(Temp_CC,   "Temperature",  patch, m_materialManager, dwi, new_dw, isNotInitialTimeStep);
      setBC(vel_CC,    "Velocity",     patch, m_materialManager, dwi, new_dw, isNotInitialTimeStep);

      for (CellIterator iter = patch->getExtraCellIterator();
           !iter.done();iter++){
        sp_vol_CC[*iter] = 1.0/rho_micro[*iter];
      }

      //__________________________________
      //    B U L L E T   P R O O F I N G
      IntVector neg_cell;
      ostringstream warn;
      if( !areAllValuesPositive( rho_CC, neg_cell ) ) {
        Point pt = patch->getCellPosition(neg_cell);
        
        warn<<"ERROR MPMICE2::actuallyInitialize, mat "<<dwi<< " cell: "
            << neg_cell << ", position: " << pt << ", rho_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__ );
      }
      if( !areAllValuesPositive( Temp_CC, neg_cell ) ) {
        Point pt = patch->getCellPosition(neg_cell);
        
        warn<<"ERROR MPMICE2::actuallyInitialize, mat "<<dwi<< " cell "
            << neg_cell << ", position: " << pt << ", Temp_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__ );
      }
      if( !areAllValuesPositive( sp_vol_CC, neg_cell ) ) {
        Point pt = patch->getCellPosition(neg_cell);
        
        warn<<"ERROR MPMICE2::actuallyInitialize, mat "<<dwi<< " cell "
            << neg_cell << ", position: " << pt << ", sp_vol_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__ );
      }
    }  //mpmMatls
  }  //patches
}

//______________________________________________________________________
//
template<typename T>
void MPMICE2::refineVariableCC(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse*,
                              DataWarehouse* new_dw,
                              const VarLabel* variable)
{
  const Level* fineLevel = getLevel(patches);
  const Level* coarseLevel = fineLevel->getCoarserLevel().get_rep();
  IntVector refineRatio(fineLevel->getRefinementRatio());

  for(int p=0;p<patches->size();p++){
    const Patch* finePatch = patches->get(p);
    ostringstream message;
    message<<"Doing refineVariableCC (" << variable->getName() << ")\t\t\t";
    printTask(patches,finePatch,cout_doing,message.str());    
    
    // region of fine space that will correspond to the coarse we need to get
    IntVector cl, ch, fl, fh;
    IntVector bl(0,0,0);  // boundary layer cells
    int nGhostCells = 1;
    bool returnExclusiveRange=true;
    
    getCoarseLevelRange(finePatch, coarseLevel, cl, ch, fl, fh, bl, 
                        nGhostCells, returnExclusiveRange);

    for(int m = 0;m<matls->size();m++){
      int indx = matls->get(m);

      CCVariable<T> fine_q_CC;
      new_dw->allocateAndPut(fine_q_CC, variable, indx, finePatch);
      
      constCCVariable<T> coarse_q_CC;

      new_dw->getRegion(coarse_q_CC, variable, indx, coarseLevel, cl, ch, false);
    
      // Only interpolate over the intersection of the fine and coarse patches
      // coarse cell 
//      linearInterpolation<T>(coarse_q_CC, coarseLevel, fineLevel,
//                             refineRatio, lo, hi, fine_q_CC);

      piecewiseConstantInterpolation<T>(coarse_q_CC, fineLevel,
                                        fl, fh, fine_q_CC);
    }
  }
}



//__________________________________
//
template<typename T>
void MPMICE2::coarsenDriver_stdNC(IntVector cl,
                                 IntVector ch,
                                 IntVector refinementRatio,
                                 double ratio,
                                 const Level* coarseLevel,
                                 constNCVariable<T>& fine_q_NC,
                                 NCVariable<T>& coarse_q_NC )
{
  T zero(0.0);
  // iterate over coarse level cells
  const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();
  Vector DX = coarseLevel->dCell();
  IntVector range(refinementRatio.x()/2, refinementRatio.y()/2, refinementRatio.z()/2);

  IntVector varLow = fine_q_NC.getLowIndex();
  IntVector varHigh = fine_q_NC.getHighIndex();

  for(NodeIterator iter(cl, ch); !iter.done(); iter++){
    IntVector c = *iter;
    IntVector fineNode = coarseLevel->mapNodeToFiner(c);
    Point coarseLoc=coarseLevel->getNodePosition(c);

    IntVector start = Max(fineNode-range, varLow);
    IntVector end = Min(fineNode+range, varHigh);

    // for each coarse level cell iterate over the fine level cells
    T q_NC_tmp(zero);

    for (NodeIterator inner(start, end); !inner.done(); inner++) {
      IntVector fc(*inner);
      Point fineLoc=fineLevel->getNodePosition(fc);
      Vector C2F = fineLoc - coarseLoc;
      Vector Vweight = C2F/DX;
      double weight = (1.-fabs(Vweight.x()))*
        (1.-fabs(Vweight.y()))*
        (1.-fabs(Vweight.z()));
      q_NC_tmp += fine_q_NC[fc]*weight;
    }
    coarse_q_NC[c] =q_NC_tmp;
  }
}

//______________________________________________________________________
//
template<typename T>
void MPMICE2::coarsenVariableCC(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset* matls,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw,
                               const VarLabel* variable,
                               T defaultValue, 
                               bool modifies,
                               string coarsenMethod)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, Ilb->timeStepLabel );

  bool isNotInitialTimeStep = (timeStep > 0);

  const Level* coarseLevel = getLevel(patches);
  const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();
  
  IntVector refineRatio(fineLevel->getRefinementRatio());
  double ratio = 1./(refineRatio.x()*refineRatio.y()*refineRatio.z());

  for(int p=0;p<patches->size();p++){
    const Patch* coarsePatch = patches->get(p);
    ostringstream message;
    message<<"Doing CoarsenVariableCC (" << variable->getName() << ")\t\t\t";
    printTask(patches,coarsePatch,cout_doing,message.str());

    for(int m = 0;m<matls->size();m++){
      int indx = matls->get(m);

      CCVariable<T> coarse_q_CC;
      if(modifies){
        new_dw->getModifiable(coarse_q_CC, variable, indx, coarsePatch);
      }else{
        new_dw->allocateAndPut(coarse_q_CC, variable, indx, coarsePatch);
      }
      coarse_q_CC.initialize(defaultValue);

      Level::selectType finePatches;
      coarsePatch->getFineLevelPatches(finePatches);
      for(unsigned int i=0;i<finePatches.size();i++){
        const Patch* finePatch = finePatches[i];
  
        IntVector cl, ch, fl, fh;
        getFineLevelRange(coarsePatch, finePatch, cl, ch, fl, fh);
        if (fh.x() <= fl.x() || fh.y() <= fl.y() || fh.z() <= fl.z()) {
          continue;
        }
        
        constCCVariable<T> fine_q_CC;
        new_dw->getRegion(fine_q_CC,  variable, indx, fineLevel, fl, fh, false);
        
        //__________________________________
        //  call the coarsening function
        ASSERT((coarsenMethod=="std" || coarsenMethod=="sum" 
                                     || coarsenMethod=="massWeighted"));
        if(coarsenMethod == "std"){
          coarsenDriver_std(cl, ch, fl, fh, refineRatio,ratio, coarseLevel, 
                            fine_q_CC, coarse_q_CC);
        }
        if(coarsenMethod =="sum"){
          ratio = 1.0;
          coarsenDriver_std(cl, ch, fl, fh, refineRatio,ratio, coarseLevel, 
                            fine_q_CC, coarse_q_CC);
        }
        if(coarsenMethod == "massWeighted"){
          constCCVariable<double> cMass;
          new_dw->getRegion(cMass,  MIlb->cMassLabel, indx, fineLevel, fl, fh, false);
          
          coarsenDriver_massWeighted(cl,ch, fl, fh, refineRatio,coarseLevel,
                                     cMass, fine_q_CC, coarse_q_CC );
        }
      }  // fine patches
      // Set BCs on coarsened data.  This sucks--Steve
      if(variable->getName()=="temp_CC"){
       setBC(coarse_q_CC, "Temperature",coarsePatch,m_materialManager,indx,new_dw, isNotInitialTimeStep);
      }
      else if(variable->getName()=="rho_CC"){
       setBC(coarse_q_CC, "Density",    coarsePatch,m_materialManager,indx,new_dw, isNotInitialTimeStep);
      }
      else if(variable->getName()=="vel_CC"){
       setBC(coarse_q_CC, "Velocity",   coarsePatch,m_materialManager,indx,new_dw, isNotInitialTimeStep);
      }
      else if(variable->getName()=="c.mass"       ||
              variable->getName()=="sp_vol_CC"    ||
              variable->getName()=="mom_L_CC"     ||
              variable->getName()=="int_eng_L_CC" ){
       setBC(coarse_q_CC,"set_if_sym_BC",coarsePatch,m_materialManager,indx,new_dw, isNotInitialTimeStep);
      }
    }  // matls
  }  // coarse level
}

//______________________________________________________________________
//
template<typename T>
void MPMICE2::coarsenVariableNC(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset* matls,
                               DataWarehouse*,
                               DataWarehouse* new_dw,
                               const VarLabel* variable,
                               T defaultValue, 
                               bool modifies,
                               string coarsenMethod)
{
  const Level* coarseLevel = getLevel(patches);
  const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();
  
  IntVector refineRatio(fineLevel->getRefinementRatio());
  double ratio = 1./(refineRatio.x()*refineRatio.y()*refineRatio.z());
  
  for(int p=0;p<patches->size();p++){  
    const Patch* coarsePatch = patches->get(p);
    ostringstream message;
    message<<"Doing CoarsenVariableNC (" << variable->getName() << ")\t\t\t";
    printTask(patches,coarsePatch,cout_doing,message.str());

    for(int m = 0;m<matls->size();m++){
      int indx = matls->get(m);

      NCVariable<T> coarse_q_NC;
      if(modifies){
        new_dw->getModifiable(coarse_q_NC, variable, indx, coarsePatch);
      }else{
        new_dw->allocateAndPut(coarse_q_NC, variable, indx, coarsePatch);
      }
      coarse_q_NC.initialize(defaultValue);

      Level::selectType finePatches;
      coarsePatch->getFineLevelPatches(finePatches);
      for(unsigned int i=0;i<finePatches.size();i++){
        const Patch* finePatch = finePatches[i];
        
        IntVector cl, ch, fl, fh;

        IntVector padding(refineRatio.x()/2,refineRatio.y()/2,refineRatio.z()/2);
        getFineLevelRangeNodes(coarsePatch, finePatch, cl, ch, fl, fh,padding);
        

        if (fh.x() <= fl.x() || fh.y() <= fl.y() || fh.z() <= fl.z()) {
          continue;
        }

        constNCVariable<T> fine_q_NC;
        new_dw->getRegion(fine_q_NC,  variable, indx, fineLevel, fl, fh, false);

        //__________________________________
        //  call the coarsening function
        ASSERT(coarsenMethod=="sum"); 
        if(coarsenMethod == "sum"){
          coarsenDriver_stdNC(cl, ch, refineRatio,ratio, coarseLevel, 
                              fine_q_NC, coarse_q_NC);
        }
      }  // fine patches
    }  // matls
  }  // coarse level
}

/* _____________________________________________________________________
MPMICE2::scheduleFinalizeTimestep--
This task called at the very bottom of the timestep,
after scheduleTimeAdvance and the scheduleCoarsen.

This is scheduled on every level.
_____________________________________________________________________*/
void
MPMICE2::scheduleFinalizeTimestep(const LevelP& level, SchedulerP& sched)
{
    cout_doing << "----------------------------" << endl;
    cout_doing << d_myworld->myRank() << " MPMICE2::scheduleFinalizeTimestep\t\t\t\tL-" << level->getIndex() << endl;

    const PatchSet* ice_patches = level->eachPatch();
    const MaterialSet* ice_matls = m_materialManager->allMaterials("ICE");
    const MaterialSet* all_matls = m_materialManager->allMaterials();
    const MaterialSet* mpm_matls = m_materialManager->allMaterials("MPM");
    const MaterialSubset* ice_matls_sub = ice_matls->getUnion();

    d_ice->scheduleConservedtoPrimitive_Vars(sched, ice_patches, ice_matls_sub,
        ice_matls,
        "finalizeTimestep");

    d_ice->scheduleTestConservation(sched, ice_patches, ice_matls_sub,
        all_matls);

    if (d_analysisModules.size() != 0) {
        vector<AnalysisModule*>::iterator iter;
        for (iter = d_analysisModules.begin();
            iter != d_analysisModules.end(); iter++) {
            AnalysisModule* am = *iter;
            am->scheduleDoAnalysis_preReloc(sched, level);
        }
    }

    // only do on finest level until we get AMR MPM
    if (level->getIndex() == level->getGrid()->numLevels() - 1)
        sched->scheduleParticleRelocation(level,
            Mlb->pXLabel_preReloc,
            d_mpm->d_particleState_preReloc,
            Mlb->pXLabel,
            d_mpm->d_particleState,
            Mlb->pParticleIDLabel, mpm_matls);

    //__________________________________
    //  on the fly analysis
    if (d_analysisModules.size() != 0) {
        vector<AnalysisModule*>::iterator iter;
        for (iter = d_analysisModules.begin();
            iter != d_analysisModules.end(); iter++) {
            AnalysisModule* am = *iter;
            am->scheduleDoAnalysis(sched, level);
        }
    }

    cout_doing << "---------------------------------------------------------" << endl;
}
