/*

The MIT License

Copyright (c) 1997-2010 Center for the Simulation of Accidental Fires and 
Explosions (CSAFE), and  Scientific Computing and Imaging Institute (SCI), 
University of Utah.

License for the specific language governing rights and limitations under
Permission is hereby granted, free of charge, to any person obtaining a 
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.

*/

#include <CCA/Components/Examples/ExamplesLabel.h>
#include <CCA/Components/Examples/RMCRT_Test.h>
#include <CCA/Components/Regridder/PerPatchVars.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/BBox.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/GeometryPiece/UnionGeometryPiece.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/SimpleMaterial.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/ProcessorGroup.h>

using SCIRun::Point;
using SCIRun::Vector;
using SCIRun::DebugStream;

static DebugStream dbg("RMCRT_Test", false);

namespace Uintah
{
RMCRT_Test::RMCRT_Test ( const ProcessorGroup* myworld ): UintahParallelComponent( myworld )
{
  //d_examplesLabel = scinew ExamplesLabel();
  d_colorLabel        = VarLabel::create("color",         CCVariable<double>::getTypeDescription());
  d_sumColorDiffLabel = VarLabel::create("sumColorDiff",      CCVariable<double>::getTypeDescription());
  d_gac = Ghost::AroundCells;
  d_gn  = Ghost::None;
  d_matl = 0;
}

RMCRT_Test::~RMCRT_Test ( void )
{
  dbg << UintahParallelComponent::d_myworld->myrank() << " Doing: RMCRT destructor " << endl;

  for (int i = 0; i< (int)d_refine_geom_objs.size(); i++) {
    delete d_refine_geom_objs[i];
  }
  //delete d_examplesLabel;
}

//______________________________________________________________________
void RMCRT_Test::problemSetup(const ProblemSpecP& prob_spec, 
                              const ProblemSpecP& restart_prob_spec, 
                              GridP& grid, 
                              SimulationStateP& state )
{
  d_sharedState = state;
  d_material = scinew SimpleMaterial();
  d_sharedState->registerSimpleMaterial( d_material );


  //manually manipulate the scheduling of copy data for the shootRay task
  Scheduler* sched = dynamic_cast<Scheduler*>(getPort("scheduler"));
  sched->overrideVariableBehavior("color",false, false, true);

  ProblemSpecP spec = prob_spec->findBlock("RMCRT");

  BBox gridBoundingBox;
  grid->getSpatialRange( gridBoundingBox );
  d_gridMax = gridBoundingBox.max().asVector();
  d_gridMin = gridBoundingBox.min().asVector();
  d_centerOfDomain   = (( d_gridMax - d_gridMin ) / 2.0 ) + d_gridMin;

  //defaults
  d_radiusOfBall     = 0.10 * d_gridMax.x();
  d_radiusOfOrbit    = 0.25 * d_gridMax.x();
  d_angularVelocity  = 10;

  spec->get("ballRadius",       d_radiusOfBall);
  spec->get("orbitRadius",      d_radiusOfOrbit);
  spec->get("angularVelocity",  d_angularVelocity);

  //__________________________________
  //  Read in the AMR section
  ProblemSpecP rmcrt_ps;
  ProblemSpecP amr_ps = prob_spec->findBlock("AMR");
  if (amr_ps){
    rmcrt_ps = amr_ps->findBlock("RMCRT");
  }

  if(!rmcrt_ps){
    string warn;
    warn ="\n INPUT FILE ERROR:\n <RMCRT>  block not found inside of <AMR> block \n";
    throw ProblemSetupException(warn, __FILE__, __LINE__);
  }

  rmcrt_ps->require( "orderOfInterpolation", d_orderOfInterpolation);

  //__________________________________
  // read in the regions that user would like 
  // refined if the grid has not been setup manually
  bool manualGrid;
  rmcrt_ps->getWithDefault("manualGrid", manualGrid, false);

  if(!manualGrid){
    ProblemSpecP refine_ps = rmcrt_ps->findBlock("Refine_Regions");
    if(!refine_ps ){
      string warn;
      warn ="\n INPUT FILE ERROR:\n <Refine_Regions> "
           " block not found inside of <RMCRT> block \n";
      throw ProblemSetupException(warn, __FILE__, __LINE__);
    }

    // Read in the refined regions geometry objects
    int piece_num = 0;
    list<GeometryObject::DataItem> geom_obj_data;
    geom_obj_data.push_back(GeometryObject::DataItem("level", GeometryObject::Integer));

    for (ProblemSpecP geom_obj_ps = refine_ps->findBlock("geom_object");
          geom_obj_ps != 0;
          geom_obj_ps = geom_obj_ps->findNextBlock("geom_object") ) {

        vector<GeometryPieceP> pieces;
        GeometryPieceFactory::create(geom_obj_ps, pieces);

        GeometryPieceP mainpiece;
        if(pieces.size() == 0){
           throw ParameterNotFound("No piece specified in geom_object", __FILE__, __LINE__);
        } else if(pieces.size() > 1){
           mainpiece = scinew UnionGeometryPiece(pieces);
        } else {
           mainpiece = pieces[0];
        }
        piece_num++;
        d_refine_geom_objs.push_back(scinew GeometryObject(mainpiece,geom_obj_ps,geom_obj_data));
     }
   }

  //__________________________________
  //  bulletproofing
  if(!d_sharedState->isLockstepAMR()){
    ostringstream msg;
    msg << "\n ERROR: You must add \n"
        << " <useLockStep> true </useLockStep> \n"
        << " inside of the <AMR> section. \n"; 
    throw ProblemSetupException(msg.str(),__FILE__, __LINE__);
  }  
}
  
//______________________________________________________________________
void RMCRT_Test::scheduleInitialize ( const LevelP& level, 
                                      SchedulerP& scheduler )
{
  printSchedule(level,dbg,"RMCRT_Test::scheduleInitialize");

  Task* task = scinew Task( "RMCRT_Test::initialize", this, 
                            &RMCRT_Test::initialize );

  task->computes( d_colorLabel );
  scheduler->addTask( task, level->eachPatch(), d_sharedState->allMaterials() );
}

//______________________________________________________________________
void RMCRT_Test::scheduleComputeStableTimestep ( const LevelP& level, SchedulerP& scheduler )
{
  printSchedule(level,dbg,"RMCRT_Test::scheduleComputeStableTimestep");

  Task* task = scinew Task( "RMCRT_Test::computeStableTimestep", this, 
                            &RMCRT_Test::computeStableTimestep );

  task->computes( d_sharedState->get_delt_label(),level.get_rep() );

  scheduler->addTask( task, level->eachPatch(), d_sharedState->allMaterials() );
}
//______________________________________________________________________
void RMCRT_Test::scheduleTimeAdvance ( const LevelP& level,
                                       SchedulerP& sched)
{
  if(level->getIndex() > 0){  // only schedule once
    return;
  }

  const MaterialSet* matls = d_sharedState->allMaterials();
  GridP grid = level->getGrid();
  int maxLevels = level->getGrid()->numLevels();

  // pseudo RMCRT
  for (int l = 0; l < maxLevels-1; l++) {
    const LevelP& level = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleShootRays( sched, patches, matls );
    //scheduleComputeQ
    //scheduleRefine_divQ
  }


  // only schedule CFD on the finest level
  const LevelP& fineLevel = grid->getLevel(maxLevels-1);
  const PatchSet* patches = fineLevel->eachPatch();
  schedulePseudoCFD( sched, patches, matls );
}
//______________________________________________________________________
//
void RMCRT_Test::schedulePseudoCFD(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls)
{
  printSchedule(patches,dbg,"RMCRT_Test::schedulePseudoCFD");
  
  Task* t = scinew Task("RMCRT_Test::pseudoCFD",
                  this, &RMCRT_Test::pseudoCFD);

  t->computes( d_colorLabel );

  sched->addTask(t, patches, matls);
}
//______________________________________________________________________
void RMCRT_Test::pseudoCFD ( const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw )
{ 
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);

    printTask(patches, patch,dbg,"pseudoCFD");

    for(int m = 0;m<matls->size();m++){
      int matl = matls->get(m);

     CCVariable<double> color;
      new_dw->allocateAndPut(color,    d_colorLabel,    matl, patch);

      for ( CellIterator iter(patch->getCellIterator()); !iter.done(); iter++) {
        IntVector c(*iter);

        Vector whereThisCellIs( patch->cellPosition( c ) );
        Vector distanceToCenterOfDomain = whereThisCellIs - d_centerOfBall;

       color[c] = distanceToCenterOfDomain.length();
      }
    }
  }
}

//______________________________________________________________________
//
void RMCRT_Test::scheduleShootRays(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls)
{
  printSchedule(patches,dbg,"RMCRT_Test::scheduleShootRays");
  
  Task* t = scinew Task("RMCRT_Test::shootRays",
                  this, &RMCRT_Test::shootRays);

  Ghost::GhostType  gn  = Ghost::None;
  Task::DomainSpec  ND  = Task::NormalDomain;
  #define allPatches 0
  #define allMatls 0

  t->requires(Task::OldDW, d_colorLabel,  allPatches, ND,allMatls, ND, gn,0);
  //t->requires(Task::OldDW, d_colorLabel,   d_gn, 0);
  t->computes( d_sumColorDiffLabel );
  sched->addTask(t, patches, matls);
}
  
//______________________________________________________________________
void RMCRT_Test::shootRays ( const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw )
{ 
  const Level* level = getLevel(patches);
  IntVector L_lo,L_hi;
  bool useBoundaryCells = false;

  level->findInteriorCellIndexRange(L_lo, L_hi);
  constCCVariable<double> color;

  old_dw->getRegion(color, d_colorLabel,  d_matl, level, L_lo, L_hi);

  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);

    printTask(patches, patch,dbg,"shootRays");

    CCVariable<double> sumColorDiff;
    new_dw->allocateAndPut(sumColorDiff, d_sumColorDiffLabel, d_matl, patch);
    sumColorDiff.initialize(0.0);

    IntVector P_lo = patch->getCellLowIndex();
    IntVector P_hi = patch->getCellHighIndex();
    IntVector middle = P_lo + (P_hi - P_lo)/IntVector(2,2,2);

    // Ray in x dir                            
    int j = middle.y();                             
    int k = middle.z();                             
    for(int i=L_lo.x(); i<L_hi.x(); i++){      
      IntVector h (i,j,k);                     
      //cout << "    " << h << endl;             
      sumColorDiff[middle] += color[h] - color[middle];  
    }                                          
    // Ray in Y dir                            
    int i = middle.x();                             
    k = middle.z();                                 
    for(int j=L_lo.y(); j<L_hi.y(); j++){      
      IntVector h (i,j,k);                     
      sumColorDiff[middle] += color[h] - color[middle];  
    }                                          

    // Ray in Z dir                            
    i = middle.x();                                 
    j = middle.y();                                 
    for(int k=L_lo.z(); k<L_hi.z(); k++){      
      IntVector h (i,j,k);                     
      sumColorDiff[middle] += color[h] - color[middle];  
    } 
      // do something here
  }
}  
//______________________________________________________________________
//
void RMCRT_Test::scheduleRefine_Q(const PatchSet* patches,
                                  SchedulerP& sched,
                                  const MaterialSet* matls)
{
  const Level* fineLevel = getLevel(patches);
  int L_indx = fineLevel->getIndex();
  
  if(L_indx > 0 ){
     printSchedule(patches,dbg,"RMCRT_Test::scheduleRefine_Q");

    Task* task = scinew Task("RMCRT_Test::refine_Q",this, 
                             &RMCRT_Test::refine_Q);
    
    Task::DomainSpec  ND  = Task::NormalDomain;
    #define allPatches 0
    #define allMatls 0
    task->requires(Task::NewDW, d_sumColorDiffLabel, allPatches, Task::CoarseLevel, allMatls, ND, d_gn,0);
     
    task->computes(d_sumColorDiffLabel);
    sched->addTask(task, patches, matls);
  }
}
  
//______________________________________________________________________
//
void RMCRT_Test::refine_Q(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse*,
                          DataWarehouse* new_dw)
{
  const Level* fineLevel = getLevel(patches);
  const Level* coarseLevel = fineLevel->getCoarserLevel().get_rep();
  
  for(int p=0;p<patches->size();p++){  
    const Patch* finePatch = patches->get(p);
    printTask(patches, finePatch,dbg,"Doing refineQ");

    Level::selectType coarsePatches;
    finePatch->getCoarseLevelPatches(coarsePatches);

    CCVariable<double> sumColorDiff_fine;
    new_dw->allocateAndPut(sumColorDiff_fine, d_sumColorDiffLabel, d_matl, finePatch);
    sumColorDiff_fine.initialize(-999);
    
    IntVector refineRatio = fineLevel->getRefinementRatio();

    // region of fine space that will correspond to the coarse we need to get
    IntVector cl, ch, fl, fh;
    IntVector bl(0,0,0);  // boundary layer or padding
    int nghostCells = 0;
    bool returnExclusiveRange=true;
    
    getCoarseLevelRange(finePatch, coarseLevel, cl, ch, fl, fh, bl, 
                        nghostCells, returnExclusiveRange);

    dbg <<" refineQ: " 
        <<" finePatch  "<< finePatch->getID() << " fl " << fl << " fh " << fh
        <<" coarseRegion " << cl << " " << ch <<endl;

    constCCVariable<double> sumColorDiff_coarse;
    new_dw->getRegion(sumColorDiff_coarse, d_sumColorDiffLabel, d_matl, coarseLevel, cl, ch);

    selectInterpolator(sumColorDiff_coarse, d_orderOfInterpolation, coarseLevel, fineLevel,
                       refineRatio, fl, fh,sumColorDiff_fine);

  }  // course patch loop 
}
  
//______________________________________________________________________
void RMCRT_Test::scheduleErrorEstimate ( const LevelP& level, SchedulerP& scheduler )
{
  printSchedule(level,dbg,"RMCRT_Test::errorEstimate");

  Task* task = scinew Task( "RMCRT_Test::errorEstimate", this, 
                            &RMCRT_Test::errorEstimate, false );

  task->requires( Task::NewDW, d_colorLabel,    d_gn, 0 );

  task->modifies( d_sharedState->get_refineFlag_label(),      d_sharedState->refineFlagMaterials() );
  task->modifies( d_sharedState->get_refinePatchFlag_label(), d_sharedState->refineFlagMaterials() );

  scheduler->addTask( task, scheduler->getLoadBalancer()->getPerProcessorPatchSet(level), d_sharedState->allMaterials() );
}

//______________________________________________________________________
void RMCRT_Test::scheduleInitialErrorEstimate ( const LevelP& level, SchedulerP& scheduler )
{
  printSchedule(level,dbg,"RMCRT_Test::scheduleInitialErrorEstimate");

  Task* task = scinew Task( "RMCRT_Test::initialErrorEstimate", this, 
                            &RMCRT_Test::errorEstimate, true );

  task->requires( Task::NewDW, d_colorLabel, d_gn, 0 );
  task->modifies( d_sharedState->get_refineFlag_label(),      d_sharedState->refineFlagMaterials() );
  task->modifies( d_sharedState->get_refinePatchFlag_label(), d_sharedState->refineFlagMaterials() );

  scheduler->addTask( task, level->eachPatch(), d_sharedState->allMaterials() );
}
//______________________________________________________________________
void RMCRT_Test::scheduleCoarsen ( const LevelP& coarseLevel, SchedulerP& scheduler )
{
  printSchedule(coarseLevel,dbg,"RMCRT_Test::scheduleCoarsen");

  Task* task = scinew Task( "RMCRT_Test::coarsen", this, 
                            &RMCRT_Test::coarsen );

  task->requires(Task::NewDW, d_colorLabel, 0, Task::FineLevel, 0, Task::NormalDomain, d_gn, 0);
  task->computes(d_colorLabel);

  scheduler->addTask( task, coarseLevel->eachPatch(), d_sharedState->allMaterials() );
}
//______________________________________________________________________
void RMCRT_Test::scheduleRefine ( const PatchSet* patches, 
                                  SchedulerP& scheduler )
{
  printSchedule(patches,dbg,"RMCRT_Test::scheduleRefine");

  Task* task = scinew Task( "RMCRT_Test::refine", this, 
                            &RMCRT_Test::refine );

  task->requires(Task::NewDW, d_colorLabel, 0, Task::CoarseLevel, 0, Task::NormalDomain, d_gn, 0);
  //    task->requires(Task::NewDW, d_oldcolorLabel, 0, Task::CoarseLevel, 0,
  //             Task::NormalDomain, d_gn, 0);
  scheduler->addTask( task, patches, d_sharedState->allMaterials() );
}
//______________________________________________________________________
void RMCRT_Test::scheduleRefineInterface ( const LevelP&, 
                                           SchedulerP&, 
                                           bool, 
                                           bool)
{
}
//______________________________________________________________________
void RMCRT_Test::initialize (const ProcessorGroup*,
                             const PatchSubset* patches, 
                             const MaterialSubset* matls,
                             DataWarehouse*, 
                             DataWarehouse* new_dw)
{

  d_centerOfBall     = d_centerOfDomain;

  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);

    for(int m = 0;m<matls->size();m++){
      int matl = matls->get(m);

     CCVariable<double> color;
      new_dw->allocateAndPut(color, d_colorLabel, matl, patch);

     for ( CellIterator iter(patch->getCellIterator()); !iter.done(); iter++) {

       IntVector idx(*iter);
        Vector whereThisCellIs( patch->cellPosition( idx ) );
        Vector distanceToCenterOfDomain = whereThisCellIs - d_centerOfBall;
        color[idx] = distanceToCenterOfDomain.length();
      }
    }
  }
}


//______________________________________________________________________
void RMCRT_Test::computeStableTimestep (const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset* /*matls*/,
                                        DataWarehouse* /*old_dw*/,
                                        DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);
  double delt = level->dCell().x();
  new_dw->put(delt_vartype(delt), d_sharedState->get_delt_label(), level);
}

//______________________________________________________________________
void RMCRT_Test::errorEstimate ( const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset* matls,
                                 DataWarehouse* old_dw, 
                                 DataWarehouse* new_dw, 
                                 bool initial )
{ 
  //__________________________________
  //   initial refinement region
  if(initial){

    const Level* level = getLevel(patches);

    for(int p=0;p<patches->size();p++){
      const Patch* patch = patches->get(p);
      printTask(patches, patch,dbg,"Doing initialErrorEstimate");

      CCVariable<int> refineFlag;
      PerPatch<PatchFlagP> refinePatchFlag;
      new_dw->getModifiable(refineFlag, d_sharedState->get_refineFlag_label(), 0, patch);
      new_dw->get(refinePatchFlag, d_sharedState->get_refinePatchFlag_label(), 0, patch);

      PatchFlag* refinePatch = refinePatchFlag.get().get_rep();

      // loop over all the geometry objects
      for(int obj=0; obj<(int)d_refine_geom_objs.size(); obj++){
        GeometryPieceP piece = d_refine_geom_objs[obj]->getPiece();
        Vector dx = patch->dCell();

        int geom_level =  d_refine_geom_objs[obj]->getInitialData_int("level");

        //don't add refinement flags if the current level is greater than the geometry level specification
        if(geom_level!=-1 && level->getIndex()>=geom_level)
          continue;

        for(CellIterator iter = patch->getCellIterator(); !iter.done();iter++){
          IntVector c = *iter;
          Point  lower  = patch->nodePosition(c);
          Vector upperV = lower.asVector() + dx; 
          Point  upper  = upperV.asPoint();

          if(piece->inside(upper) && piece->inside(lower))
            refineFlag[c] = true;
            refinePatch->set();
        }
      }  // object loop
    }  // patches loop
  }
  else{

    //__________________________________
    //   flag regions inside of the ball
    double pi = 3.141592653589;
    if ( getLevel(patches)->getIndex() == getLevel(patches)->getGrid()->numLevels()-1 ) {
      d_centerOfBall     = d_centerOfDomain;
    }

    for ( int p = 0; p < patches->size(); p++ ) {
      const Patch* patch = patches->get(p);

      printTask(patches, patch,dbg,"Doing errorEstimate");

      CCVariable<int> refineFlag;
      new_dw->getModifiable(refineFlag, d_sharedState->get_refineFlag_label(), 0, patch);

      PerPatch<PatchFlagP> refinePatchFlag;
      new_dw->get(refinePatchFlag, d_sharedState->get_refinePatchFlag_label(), 0, patch);
      PatchFlag* refinePatch = refinePatchFlag.get().get_rep();

      bool foundErrorOnPatch = false;
      constCCVariable<double> color;
      new_dw->get( color,    d_colorLabel,    d_matl, patch, d_gn, 0 );

      for(CellIterator iter = patch->getCellIterator(); !iter.done();iter++){
        IntVector c = *iter;

        if ( color[c] <= d_radiusOfBall ) {
          refineFlag[c]=true;
          foundErrorOnPatch = true;
        } else {
          refineFlag[c]=false;
        }
      }

      // flag this patch
      if ( foundErrorOnPatch ) {
        refinePatch->flag = true;
      } else {
        refinePatch->flag = false;
      }

    }  // patch loop
  }  // not initial timestep
}
//______________________________________________________________________
void RMCRT_Test::coarsen ( const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse*, DataWarehouse* new_dw )
{
  const Level* coarseLevel = getLevel(patches);
  const LevelP fineLevel = coarseLevel->getFinerLevel();
  IntVector rr(fineLevel->getRefinementRatio());
  double ratio = 1./(rr.x()*rr.y()*rr.z());

  for(int p=0;p<patches->size();p++){  
    const Patch* coarsePatch = patches->get(p);

    printTask(patches, coarsePatch,dbg,"Doing coarsen");

    // Find the overlapping regions...
    Level::selectType finePatches;
    coarsePatch->getFineLevelPatches(finePatches);

    for(int m = 0;m<matls->size();m++){
      int matl = matls->get(m);

      CCVariable<double> color_coarse;
      new_dw->allocateAndPut(color_coarse, d_colorLabel, matl, coarsePatch);

      for(int i=0;i<finePatches.size();i++){
        const Patch* finePatch = finePatches[i];

        constCCVariable<double> color_fine;
        new_dw->get(color_fine, d_colorLabel, matl, finePatch, d_gn, 0);

        IntVector fl(finePatch->getCellLowIndex());
        IntVector fh(finePatch->getCellHighIndex());

        IntVector l(fineLevel->mapCellToCoarser(fl));
        IntVector h(fineLevel->mapCellToCoarser(fh));

        l = Max(l, coarsePatch->getCellLowIndex());
        h = Min(h, coarsePatch->getCellHighIndex());

        for(CellIterator iter(l, h); !iter.done(); iter++){
          IntVector c = *iter;

          double sumColorDiff=0;
          IntVector fineStart(coarseLevel->mapCellToFiner(c));

          for(CellIterator inside(IntVector(0,0,0), fineLevel->getRefinementRatio());
              !inside.done(); inside++){
            sumColorDiff += color_fine[fineStart+*inside];
          }
          color_coarse[c]=sumColorDiff*ratio;
        }  // intersection loop
      }  // fine patch loop
    }
  }  // course patch loop 
}

//______________________________________________________________________
void RMCRT_Test::refine ( const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse*, 
                          DataWarehouse* new_dw )
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);

    printTask(patches, patch,dbg,"Doing refine");

    for(int m = 0;m<matls->size();m++){
      int matl = matls->get(m);

     CCVariable<double> color;
      new_dw->allocateAndPut(color, d_colorLabel, matl, patch);

     for ( CellIterator iter(patch->getCellIterator()); !iter.done(); iter++) {
        IntVector c = *iter;
        Vector whereThisCellIs( patch->cellPosition( c ) );
        Vector distanceToCenterOfDomain = whereThisCellIs - d_centerOfBall;
        color[c] = distanceToCenterOfDomain.length();
      }
    }
  }
}
//__________________________________
//  
void RMCRT_Test::printSchedule(const PatchSet* patches,
                              DebugStream& dbg,
                              const string& where)
{
  if (dbg.active()){
    dbg << UintahParallelComponent::d_myworld->myrank() << " ";
    dbg << left;
    dbg.width(50);
    dbg  << where << "L-"
        << getLevel(patches)->getIndex()<< endl;
  }  
}
//__________________________________
//
void RMCRT_Test::printSchedule(const LevelP& level,
                              DebugStream& dbg,
                              const string& where)
{
  if (dbg.active()){
    dbg << UintahParallelComponent::d_myworld->myrank() << " ";
    dbg << left;
    dbg.width(50);
    dbg << where << "L-"
        << level->getIndex()<< endl;
  }  
}
//__________________________________
//
void RMCRT_Test::printTask(const PatchSubset* patches,
                          const Patch* patch,
                          DebugStream& dbg,
                          const string& where)
{
  if (dbg.active()){
    dbg << UintahParallelComponent::d_myworld->myrank() << " ";
    dbg << left;
    dbg.width(50);
    dbg << where << " MPM \tL-"
        << getLevel(patches)->getIndex()
        << " patch " << patch->getGridIndex()<< endl;
  }  
}

//__________________________________
//
void RMCRT_Test::printTask(const Patch* patch,
                          DebugStream& dbg,
                          const string& where)
{
  if (dbg.active()){
    dbg << UintahParallelComponent::d_myworld->myrank() << " ";
    dbg << left;
    dbg.width(50);
    dbg << where << " MPM \tL-"
        << patch->getLevel()->getIndex()
        << " patch " << patch->getGridIndex()<< endl;
  }  
}  
 
} // namespace Uintah
