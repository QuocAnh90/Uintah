//----- ScalarSolver.cc ----------------------------------------------

#include <Packages/Uintah/CCA/Components/Arches/ScalarSolver.h>
#include <Packages/Uintah/CCA/Components/Arches/Arches.h>
#include <Packages/Uintah/CCA/Components/Arches/ArchesLabel.h>
#include <Packages/Uintah/CCA/Components/Arches/ArchesMaterial.h>
#include <Packages/Uintah/CCA/Components/Arches/BoundaryCondition.h>
#include <Packages/Uintah/CCA/Components/Arches/CellInformationP.h>
#include <Packages/Uintah/CCA/Components/Arches/Discretization.h>
#include <Packages/Uintah/CCA/Components/Arches/PetscSolver.h>
#include <Packages/Uintah/CCA/Components/Arches/PhysicalConstants.h>
#include <Packages/Uintah/CCA/Components/Arches/RBGSSolver.h>
#include <Packages/Uintah/CCA/Components/Arches/Source.h>
#include <Packages/Uintah/CCA/Components/Arches/TurbulenceModel.h>
#include <Packages/Uintah/CCA/Ports/DataWarehouse.h>
#include <Packages/Uintah/CCA/Ports/Scheduler.h>
#include <Packages/Uintah/Core/Exceptions/InvalidValue.h>
#include <Packages/Uintah/Core/Grid/CCVariable.h>
#include <Packages/Uintah/Core/Grid/Level.h>
#include <Packages/Uintah/Core/Grid/Patch.h>
#include <Packages/Uintah/Core/Grid/PerPatch.h>
#include <Packages/Uintah/Core/Grid/SFCXVariable.h>
#include <Packages/Uintah/Core/Grid/SFCYVariable.h>
#include <Packages/Uintah/Core/Grid/SFCZVariable.h>
#include <Packages/Uintah/Core/Grid/SimulationState.h>
#include <Packages/Uintah/Core/Grid/Task.h>
#include <Packages/Uintah/Core/Grid/VarTypes.h>
#include <Packages/Uintah/Core/ProblemSpec/ProblemSpec.h>


using namespace Uintah;
using namespace std;

//****************************************************************************
// Default constructor for PressureSolver
//****************************************************************************
ScalarSolver::ScalarSolver(const ArchesLabel* label,
			   const MPMArchesLabel* MAlb,
			   TurbulenceModel* turb_model,
			   BoundaryCondition* bndry_cond,
			   PhysicalConstants* physConst) :
                                 d_lab(label), d_MAlab(MAlb),
                                 d_turbModel(turb_model), 
                                 d_boundaryCondition(bndry_cond),
				 d_physicalConsts(physConst)
{
  d_discretize = 0;
  d_source = 0;
  d_linearSolver = 0;
}

//****************************************************************************
// Destructor
//****************************************************************************
ScalarSolver::~ScalarSolver()
{
  delete d_discretize;
  delete d_source;
  delete d_linearSolver;
}

//****************************************************************************
// Problem Setup
//****************************************************************************
void 
ScalarSolver::problemSetup(const ProblemSpecP& params)
{
  ProblemSpecP db = params->findBlock("MixtureFractionSolver");
  string finite_diff;
  db->require("finite_difference", finite_diff);
  if (finite_diff == "second") 
    d_discretize = scinew Discretization();
  else {
    throw InvalidValue("Finite Differencing scheme "
		       "not supported: " + finite_diff);
    //throw InvalidValue("Finite Differencing scheme "
	//	       "not supported: " + finite_diff, db);
  }
  // make source and boundary_condition objects
  d_source = scinew Source(d_turbModel, d_physicalConsts);
  string linear_sol;
  db->require("linear_solver", linear_sol);
  if (linear_sol == "linegs")
    d_linearSolver = scinew RBGSSolver();
  else if (linear_sol == "petsc")
     d_linearSolver = scinew PetscSolver(0); // CHEAT - steve d_myworld);
  else {
    throw InvalidValue("linear solver option"
		       " not supported" + linear_sol);
    //throw InvalidValue("linear solver option"
	//	       " not supported" + linear_sol, db);
  }
  d_linearSolver->problemSetup(db);
}

//****************************************************************************
// Schedule solve of linearized scalar equation
//****************************************************************************
void 
ScalarSolver::solve(SchedulerP& sched,
		    const PatchSet* patches,
		    const MaterialSet* matls,
		    int index)
{
  //create a new data warehouse to store matrix coeff
  // and source terms. It gets reinitialized after every 
  // pressure solve.
   //DataWarehouseP matrix_dw = sched->createDataWarehouse(new_dw);

  //computes stencil coefficients and source terms
  // requires : scalarIN, [u,v,w]VelocitySPBC, densityIN, viscosityIN
  // computes : scalCoefSBLM, scalLinSrcSBLM, scalNonLinSrcSBLM
  sched_buildLinearMatrix(sched, patches, matls, index);
  
  // Schedule the scalar solve
  // require : scalarIN, scalCoefSBLM, scalNonLinSrcSBLM
  // compute : scalResidualSS, scalCoefSS, scalNonLinSrcSS, scalarSP
  //d_linearSolver->sched_scalarSolve(level, sched, new_dw, matrix_dw, index);
  sched_scalarLinearSolve(sched, patches, matls, index);
}

//****************************************************************************
// Schedule build of linear matrix
//****************************************************************************
void 
ScalarSolver::sched_buildLinearMatrix(SchedulerP& sched, const PatchSet* patches,
				      const MaterialSet* matls,
				      int index)
{
  Task* tsk = scinew Task("ScalarSolver::BuildCoeff",
			  this,
			  &ScalarSolver::buildLinearMatrix,
			  index);

  // This task requires scalar and density from old time step for transient
  // calculation
  //DataWarehouseP old_dw = new_dw->getTop();
  tsk->requires(Task::OldDW, d_lab->d_sharedState->get_delt_label());
  
  tsk->requires(Task::NewDW, d_lab->d_cellTypeLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_scalarINLabel,
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalarOUTBCLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_densityINLabel, 
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_viscosityINLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_uVelocityOUTBCLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_vVelocityOUTBCLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_wVelocityOUTBCLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);

      // added one more argument of index to specify scalar component
  tsk->computes(d_lab->d_scalCoefSBLMLabel, d_lab->d_stencilMatl,
		Task::OutOfDomain);
  tsk->computes(d_lab->d_scalNonLinSrcSBLMLabel);

  sched->addTask(tsk, patches, matls);
}

//****************************************************************************
// Schedule linear solve of scalar
//****************************************************************************
void
ScalarSolver::sched_scalarLinearSolve(SchedulerP& sched, const PatchSet* patches,
				      const MaterialSet* matls,
				      int index)
{
  Task* tsk = scinew Task("ScalarSolver::scalarLinearSolve",
			  this,
			  &ScalarSolver::scalarLinearSolve, index);
  
  tsk->requires(Task::OldDW, d_lab->d_sharedState->get_delt_label());
  
  // coefficient for the variable for which solve is invoked
  tsk->requires(Task::NewDW, d_lab->d_densityINLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalarOUTBCLabel, 
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_scalCoefSBLMLabel, 
		d_lab->d_stencilMatl, Task::OutOfDomain,
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalNonLinSrcSBLMLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->computes(d_lab->d_scalarSPLabel);
  
  sched->addTask(tsk, patches, matls);
}
      
//****************************************************************************
// Actually build linear matrix
//****************************************************************************
void ScalarSolver::buildLinearMatrix(const ProcessorGroup* pc,
				     const PatchSubset* patches,
				     const MaterialSubset*,
				     DataWarehouse* old_dw,
				     DataWarehouse* new_dw,
				     int index)
{
  delt_vartype delT;
  old_dw->get(delT, d_lab->d_sharedState->get_delt_label() );
  double delta_t = delT;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int archIndex = 0; // only one arches material
    int matlIndex = d_lab->d_sharedState->getArchesMaterial(archIndex)->getDWIndex(); 
    ArchesVariables scalarVars;
    
    // Get the PerPatch CellInformation data
    PerPatch<CellInformationP> cellInfoP;
    // checkpointing
    //  old_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    // new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    if (new_dw->exists(d_lab->d_cellInfoLabel, matlIndex, patch)) 
      new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    else {
      cellInfoP.setData(scinew CellInformation(patch));
      new_dw->put(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    }
    CellInformation* cellinfo = cellInfoP.get().get_rep();

    // from old_dw get PCELL, DENO, FO(index)
    new_dw->getCopy(scalarVars.cellType, d_lab->d_cellTypeLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.old_scalar, d_lab->d_scalarINLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);

    // from new_dw get DEN, VIS, F(index), U, V, W
    new_dw->getCopy(scalarVars.density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.viscosity, d_lab->d_viscosityINLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.scalar, d_lab->d_scalarOUTBCLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    // for explicit get old values
    new_dw->getCopy(scalarVars.uVelocity, d_lab->d_uVelocityOUTBCLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.vVelocity, d_lab->d_vVelocityOUTBCLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.wVelocity, d_lab->d_wVelocityOUTBCLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);

  // allocate matrix coeffs
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++) {
      new_dw->allocate(scalarVars.scalarCoeff[ii], 
		       d_lab->d_scalCoefSBLMLabel, ii, patch);
      new_dw->allocate(scalarVars.scalarConvectCoeff[ii],
		       d_lab->d_scalConvCoefSBLMLabel, ii, patch);
    }
    new_dw->allocate(scalarVars.scalarLinearSrc, 
		     d_lab->d_scalLinSrcSBLMLabel, matlIndex, patch);
    new_dw->allocate(scalarVars.scalarNonlinearSrc, 
		     d_lab->d_scalNonLinSrcSBLMLabel, matlIndex, patch);
 
  // compute ith component of scalar stencil coefficients
  // inputs : scalarSP, [u,v,w]VelocityMS, densityCP, viscosityCTS
  // outputs: scalCoefSBLM
    d_discretize->calculateScalarCoeff(pc, patch,
				       delta_t, index, cellinfo, 
				       &scalarVars);

    // Calculate scalar source terms
    // inputs : [u,v,w]VelocityMS, scalarSP, densityCP, viscosityCTS
    // outputs: scalLinSrcSBLM, scalNonLinSrcSBLM
    d_source->calculateScalarSource(pc, patch,
				    delta_t, index, cellinfo, 
				    &scalarVars );

    // Calculate the scalar boundary conditions
    // inputs : scalarSP, scalCoefSBLM
    // outputs: scalCoefSBLM
    d_boundaryCondition->scalarBC(pc, patch,  index, cellinfo, 
				  &scalarVars);
  // apply multimaterial intrusion wallbc
    if (d_MAlab)
      d_boundaryCondition->mmscalarWallBC(pc, patch, cellinfo,
					  &scalarVars);

    // similar to mascal
    // inputs :
    // outputs:
    d_source->modifyScalarMassSource(pc, patch, delta_t, index, &scalarVars);
    
    // Calculate the scalar diagonal terms
    // inputs : scalCoefSBLM, scalLinSrcSBLM
    // outputs: scalCoefSBLM
    d_discretize->calculateScalarDiagonal(pc, patch, index, &scalarVars);
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++) {
      new_dw->put(scalarVars.scalarCoeff[ii], 
		  d_lab->d_scalCoefSBLMLabel, ii, patch);
    }
    new_dw->put(scalarVars.scalarNonlinearSrc, 
		d_lab->d_scalNonLinSrcSBLMLabel, matlIndex, patch);

  }
}
//****************************************************************************
// Actual scalar solve .. may be changed after recursive tasks are added
//****************************************************************************
void 
ScalarSolver::scalarLinearSolve(const ProcessorGroup* pc,
				const PatchSubset* patches,
				const MaterialSubset*,
				DataWarehouse* old_dw,
				DataWarehouse* new_dw,
				int index)
{
  delt_vartype delT;
  old_dw->get(delT, d_lab->d_sharedState->get_delt_label() );
  double delta_t = delT;
  
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int archIndex = 0; // only one arches material
    int matlIndex = d_lab->d_sharedState->getArchesMaterial(archIndex)->getDWIndex(); 
    ArchesVariables scalarVars;
    //DataWarehouseP old_dw = new_dw->getTop();
    // Get the PerPatch CellInformation data
    PerPatch<CellInformationP> cellInfoP;
    // get old_dw from getTop function
    // checkpointing
    //  old_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    //  new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    //  old_dw->get(cellInfoP, d_cellInfoLabel, matlIndex, patch);
    if (new_dw->exists(d_lab->d_cellInfoLabel, matlIndex, patch)) 
      new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    else {
      cellInfoP.setData(scinew CellInformation(patch));
      new_dw->put(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    }
    CellInformation* cellinfo = cellInfoP.get().get_rep();
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    // for explicit calculation
    {
    new_dw->allocate(scalarVars.scalar, d_lab->d_scalarSPLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->copyOut(scalarVars.scalar, d_lab->d_scalarOUTBCLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    }
    scalarVars.old_scalar.allocate(scalarVars.scalar.getLowIndex(),
				   scalarVars.scalar.getHighIndex());
    scalarVars.old_scalar.copy(scalarVars.scalar);
    
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++)
      new_dw->getCopy(scalarVars.scalarCoeff[ii], d_lab->d_scalCoefSBLMLabel, 
		  ii, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.scalarNonlinearSrc, d_lab->d_scalNonLinSrcSBLMLabel,
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->allocate(scalarVars.residualScalar, d_lab->d_scalarRes,
		     matlIndex, patch);

#if 0  
  // compute eqn residual
    d_linearSolver->computeScalarResidual(pc, patch, new_dw, new_dw, index, 
					  &scalarVars);
    new_dw->put(sum_vartype(scalarVars.residScalar), d_lab->d_scalarResidLabel);
    new_dw->put(sum_vartype(scalarVars.truncScalar), d_lab->d_scalarTruncLabel);
#endif
  // apply underelax to eqn
    d_linearSolver->computeScalarUnderrelax(pc, patch, index, 
					    &scalarVars);
    // make it a separate task later
    d_linearSolver->scalarLisolve(pc, patch, index, delta_t, 
				  &scalarVars, cellinfo, d_lab);
  // put back the results
    new_dw->put(scalarVars.scalar, d_lab->d_scalarSPLabel, 
		matlIndex, patch);
  }
}

//****************************************************************************
// Schedule solve of linearized scalar equation
//****************************************************************************
void 
ScalarSolver::solvePred(SchedulerP& sched,
			const PatchSet* patches,
			const MaterialSet* matls,
			int index)
{
  //computes stencil coefficients and source terms
  // requires : scalarIN, [u,v,w]VelocitySPBC, densityIN, viscosityIN
  // computes : scalCoefSBLM, scalLinSrcSBLM, scalNonLinSrcSBLM
  sched_buildLinearMatrixPred(sched, patches, matls, index);
  
  // Schedule the scalar solve
  // require : scalarIN, scalCoefSBLM, scalNonLinSrcSBLM
  // compute : scalResidualSS, scalCoefSS, scalNonLinSrcSS, scalarSP
  //d_linearSolver->sched_scalarSolve(level, sched, new_dw, matrix_dw, index);
  sched_scalarLinearSolvePred(sched, patches, matls, index);
}

//****************************************************************************
// Schedule build of linear matrix
//****************************************************************************
void 
ScalarSolver::sched_buildLinearMatrixPred(SchedulerP& sched,
					  const PatchSet* patches,
					  const MaterialSet* matls,
					  int index)
{
  Task* tsk = scinew Task("ScalarSolver::BuildCoeffPred",
			  this,
			  &ScalarSolver::buildLinearMatrixPred,
			  index);


  tsk->requires(Task::OldDW, d_lab->d_sharedState->get_delt_label());
  
  // This task requires scalar and density from old time step for transient
  // calculation
  //DataWarehouseP old_dw = new_dw->getTop();  
  tsk->requires(Task::NewDW, d_lab->d_cellTypeLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_scalarOUTBCLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_densityINLabel, 
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_viscosityINLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_uVelocityOUTBCLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_vVelocityOUTBCLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_wVelocityOUTBCLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);

      // added one more argument of index to specify scalar component
  tsk->computes(d_lab->d_scalCoefPredLabel, d_lab->d_stencilMatl,
		Task::OutOfDomain);
  tsk->computes(d_lab->d_scalDiffCoefPredLabel, d_lab->d_stencilMatl,
		Task::OutOfDomain);
  tsk->computes(d_lab->d_scalNonLinSrcPredLabel);

  sched->addTask(tsk, patches, matls);
}

      
//****************************************************************************
// Actually build linear matrix
//****************************************************************************
void ScalarSolver::buildLinearMatrixPred(const ProcessorGroup* pc,
					 const PatchSubset* patches,
					 const MaterialSubset*,
					 DataWarehouse* old_dw,
					 DataWarehouse* new_dw,
					 int index)
{
  delt_vartype delT;
  old_dw->get(delT, d_lab->d_sharedState->get_delt_label() );
  double delta_t = delT;
#ifdef correctorstep
  delta_t /= 2.0;
#endif
#ifdef Runge_Kutta_2nd
  delta_t *= 2.0; 
#endif
#ifdef Runge_Kutta_3d
  double gamma_1 = 8.0/15.0;
  delta_t *= 2.0; // since correctorstep is also defined for Runge-Kutta
  delta_t *= gamma_1; 
#ifdef Runge_Kutta_3d_ssp
  delta_t /= gamma_1; 
#endif
#endif

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int archIndex = 0; // only one arches material
    int matlIndex = d_lab->d_sharedState->
                    getArchesMaterial(archIndex)->getDWIndex(); 
    ArchesVariables scalarVars;
    
    // Get the PerPatch CellInformation data
    PerPatch<CellInformationP> cellInfoP;
    // checkpointing
    //  old_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    // new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    if (new_dw->exists(d_lab->d_cellInfoLabel, matlIndex, patch)) 
      new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    else {
      cellInfoP.setData(scinew CellInformation(patch));
      new_dw->put(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    }
    CellInformation* cellinfo = cellInfoP.get().get_rep();

    // from old_dw get PCELL, DENO, FO(index)
    new_dw->getCopy(scalarVars.cellType, d_lab->d_cellTypeLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    //    new_dw->get(scalarVars.old_scalar, d_lab->d_scalarINLabel, 
    //		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.old_scalar, d_lab->d_scalarOUTBCLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);

    // from new_dw get DEN, VIS, F(index), U, V, W
    new_dw->getCopy(scalarVars.density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.viscosity, d_lab->d_viscosityINLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.scalar, d_lab->d_scalarOUTBCLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    // for explicit get old values
    new_dw->getCopy(scalarVars.uVelocity, d_lab->d_uVelocityOUTBCLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.vVelocity, d_lab->d_vVelocityOUTBCLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.wVelocity, d_lab->d_wVelocityOUTBCLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);

  // allocate matrix coeffs
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++) {
      new_dw->allocate(scalarVars.scalarCoeff[ii], 
		       d_lab->d_scalCoefPredLabel, ii, patch);
      new_dw->allocate(scalarVars.scalarConvectCoeff[ii],
		       d_lab->d_scalConvCoefPredLabel, ii, patch);
      new_dw->allocate(scalarVars.scalarDiffusionCoeff[ii],
		       d_lab->d_scalDiffCoefPredLabel, ii, patch);
    }
    new_dw->allocate(scalarVars.scalarLinearSrc, 
		     d_lab->d_scalLinSrcPredLabel, matlIndex, patch);
    new_dw->allocate(scalarVars.scalarNonlinearSrc, 
		     d_lab->d_scalNonLinSrcPredLabel, matlIndex, patch);
 
  // compute ith component of scalar stencil coefficients
  // inputs : scalarSP, [u,v,w]VelocityMS, densityCP, viscosityCTS
  // outputs: scalCoefSBLM
    d_discretize->calculateScalarCoeff(pc, patch,
				       delta_t, index, cellinfo, 
				       &scalarVars);

    // Calculate scalar source terms
    // inputs : [u,v,w]VelocityMS, scalarSP, densityCP, viscosityCTS
    // outputs: scalLinSrcSBLM, scalNonLinSrcSBLM
    d_source->calculateScalarSource(pc, patch,
				    delta_t, index, cellinfo, 
				    &scalarVars );

    // Calculate the scalar boundary conditions
    // inputs : scalarSP, scalCoefSBLM
    // outputs: scalCoefSBLM
    d_boundaryCondition->scalarBC(pc, patch,  index, cellinfo, 
				  &scalarVars);
  // apply multimaterial intrusion wallbc
    if (d_MAlab)
      d_boundaryCondition->mmscalarWallBC(pc, patch, cellinfo,
					  &scalarVars);

    // similar to mascal
    // inputs :
    // outputs:
    d_source->modifyScalarMassSource(pc, patch, delta_t, index, &scalarVars);
    
    // Calculate the scalar diagonal terms
    // inputs : scalCoefSBLM, scalLinSrcSBLM
    // outputs: scalCoefSBLM
    d_discretize->calculateScalarDiagonal(pc, patch, index, &scalarVars);
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++) {
      new_dw->put(scalarVars.scalarCoeff[ii], 
		  d_lab->d_scalCoefPredLabel, ii, patch);
      new_dw->put(scalarVars.scalarDiffusionCoeff[ii],
		  d_lab->d_scalDiffCoefPredLabel, ii, patch);
    }
    new_dw->put(scalarVars.scalarNonlinearSrc, 
		d_lab->d_scalNonLinSrcPredLabel, matlIndex, patch);

  }
}


//****************************************************************************
// Schedule linear solve of scalar
//****************************************************************************
void
ScalarSolver::sched_scalarLinearSolvePred(SchedulerP& sched,
					  const PatchSet* patches,
					  const MaterialSet* matls,
					  int index)
{
  Task* tsk = scinew Task("ScalarSolver::scalarLinearSolvePred",
			  this,
			  &ScalarSolver::scalarLinearSolvePred,
			  index);
  
  tsk->requires(Task::OldDW, d_lab->d_sharedState->get_delt_label());

  // coefficient for the variable for which solve is invoked
  tsk->requires(Task::NewDW, d_lab->d_densityINLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalarOUTBCLabel, 
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_scalCoefPredLabel, 
		d_lab->d_stencilMatl, Task::OutOfDomain,
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalNonLinSrcPredLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
#ifdef correctorstep
  tsk->computes(d_lab->d_scalarPredLabel);
#else
  tsk->computes(d_lab->d_scalarSPLabel);
#endif
  
#ifdef Runge_Kutta_3d
#ifndef Runge_Kutta_3d_ssp
  tsk->computes(d_lab->d_scalarTempLabel);
#endif
#endif
  
  sched->addTask(tsk, patches, matls);
}
//****************************************************************************
// Actual scalar solve .. may be changed after recursive tasks are added
//****************************************************************************
void 
ScalarSolver::scalarLinearSolvePred(const ProcessorGroup* pc,
                                const PatchSubset* patches,
				const MaterialSubset*,
				DataWarehouse* old_dw,
				DataWarehouse* new_dw,
				int index)
{
  delt_vartype delT;
  old_dw->get(delT, d_lab->d_sharedState->get_delt_label() );
  double delta_t = delT;
#ifdef correctorstep
  delta_t /= 2.0;
#endif
#ifdef Runge_Kutta_2nd
  delta_t *= 2.0; 
#endif
#ifdef Runge_Kutta_3d
  double gamma_1 = 8.0/15.0;
  delta_t *= 2.0; // since correctorstep is also defined for Runge-Kutta
  delta_t *= gamma_1; 
#ifdef Runge_Kutta_3d_ssp
  delta_t /= gamma_1; 
#endif
#endif
  
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int archIndex = 0; // only one arches material
    int matlIndex = d_lab->d_sharedState->
                    getArchesMaterial(archIndex)->getDWIndex(); 
    ArchesVariables scalarVars;
    //DataWarehouseP old_dw = new_dw->getTop();
    // Get the PerPatch CellInformation data
    PerPatch<CellInformationP> cellInfoP;
    // get old_dw from getTop function
    // checkpointing
    //  old_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    //  new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    //  old_dw->get(cellInfoP, d_cellInfoLabel, matlIndex, patch);
    if (new_dw->exists(d_lab->d_cellInfoLabel, matlIndex, patch)) 
      new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    else {
      cellInfoP.setData(scinew CellInformation(patch));
      new_dw->put(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    }
    CellInformation* cellinfo = cellInfoP.get().get_rep();
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    // for explicit calculation
    {
#ifdef correctorstep
    new_dw->allocate(scalarVars.scalar, d_lab->d_scalarPredLabel, 
                matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
#else
    new_dw->allocate(scalarVars.scalar, d_lab->d_scalarSPLabel, 
                matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
#endif
    new_dw->copyOut(scalarVars.scalar, d_lab->d_scalarOUTBCLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    }
    scalarVars.old_scalar.allocate(scalarVars.scalar.getLowIndex(),
				   scalarVars.scalar.getHighIndex());
    scalarVars.old_scalar.copy(scalarVars.scalar);
    
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++)
      new_dw->getCopy(scalarVars.scalarCoeff[ii], d_lab->d_scalCoefPredLabel, 
		  ii, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.scalarNonlinearSrc, d_lab->d_scalNonLinSrcPredLabel,
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->allocate(scalarVars.residualScalar, d_lab->d_scalarRes,
		     matlIndex, patch);

  // apply underelax to eqn
    d_linearSolver->computeScalarUnderrelax(pc, patch, index, 
					    &scalarVars);
    // make it a separate task later
    d_linearSolver->scalarLisolve(pc, patch, index, delta_t, 
    &scalarVars, cellinfo, d_lab);
				  // put back the results
#if 0
    cerr << "print scalar solve after predict" << endl;
    scalarVars.scalar.print(cerr);
#endif

#ifdef Runge_Kutta_3d
#ifndef Runge_Kutta_3d_ssp
    CCVariable<double> temp_scalar;
    constCCVariable<double> old_density;

    new_dw->get(old_density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->allocate(temp_scalar, d_lab->d_scalarTempLabel, matlIndex, patch);
    temp_scalar.initialize(0.0);
    
    IntVector indexLow = patch->getCellFORTLowIndex();
    IntVector indexHigh = patch->getCellFORTHighIndex();
    
    for (int colZ = indexLow.z(); colZ <= indexHigh.z(); colZ ++) {
      for (int colY = indexLow.y(); colY <= indexHigh.y(); colY ++) {
        for (int colX = indexLow.x(); colX <= indexHigh.x(); colX ++) {

            IntVector currCell(colX, colY, colZ);

            temp_scalar[currCell] = old_density[currCell]*
	    (scalarVars.scalar[currCell]-
            scalarVars.old_scalar[currCell])/gamma_1;
        }
      }
    }
    new_dw->put(temp_scalar, d_lab->d_scalarTempLabel, matlIndex, patch);
#endif
#endif

#ifdef correctorstep
    new_dw->put(scalarVars.scalar, d_lab->d_scalarPredLabel, 
                matlIndex, patch);
#else
    new_dw->put(scalarVars.scalar, d_lab->d_scalarSPLabel, 
                matlIndex, patch);
#endif

  }
}

//****************************************************************************
// Schedule solve of linearized scalar equation, corrector step
//****************************************************************************
void 
ScalarSolver::solveCorr(SchedulerP& sched,
			const PatchSet* patches,
			const MaterialSet* matls,
			int index)
{
  //computes stencil coefficients and source terms
  // requires : scalarIN, [u,v,w]VelocitySPBC, densityIN, viscosityIN
  // computes : scalCoefSBLM, scalLinSrcSBLM, scalNonLinSrcSBLM
  sched_buildLinearMatrixCorr(sched, patches, matls, index);
  
  // Schedule the scalar solve
  // require : scalarIN, scalCoefSBLM, scalNonLinSrcSBLM
  // compute : scalResidualSS, scalCoefSS, scalNonLinSrcSS, scalarSP
  //d_linearSolver->sched_scalarSolve(level, sched, new_dw, matrix_dw, index);
  sched_scalarLinearSolveCorr(sched, patches, matls, index);
}

//****************************************************************************
// Schedule build of linear matrix
//****************************************************************************
void 
ScalarSolver::sched_buildLinearMatrixCorr(SchedulerP& sched,
					  const PatchSet* patches,
					  const MaterialSet* matls,
					  int index)
{
  Task* tsk = scinew Task("ScalarSolver::BuildCoeffCorr",
			  this,
			  &ScalarSolver::buildLinearMatrixCorr,
			  index);


  tsk->requires(Task::OldDW, d_lab->d_sharedState->get_delt_label());
  
  // This task requires scalar and density from old time step for transient
  // calculation
  //DataWarehouseP old_dw = new_dw->getTop();
  tsk->requires(Task::NewDW, d_lab->d_cellTypeLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  #ifdef Runge_Kutta_3d
//  tsk->requires(Task::NewDW, d_lab->d_scalarPredLabel,
//		Ghost::AroundCells, Arches::ONEGHOSTCELL);
//  tsk->requires(Task::NewDW, d_lab->d_densityPredLabel, 
//		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalarIntermLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_densityIntermLabel, 
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_viscosityIntermLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_uVelocityIntermLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_vVelocityIntermLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_wVelocityIntermLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  #else
  #ifndef Runge_Kutta_2nd
  //tsk->requires(Task::NewDW, d_lab->d_scalarINLabel,
  //		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalarOUTBCLabel,
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_densityINLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  #endif
  tsk->requires(Task::NewDW, d_lab->d_scalarPredLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_densityPredLabel, 
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_viscosityPredLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_uVelocityPredLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_vVelocityPredLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_wVelocityPredLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  #endif

      // added one more argument of index to specify scalar component
  tsk->computes(d_lab->d_scalCoefCorrLabel, d_lab->d_stencilMatl,
		Task::OutOfDomain);
  tsk->computes(d_lab->d_scalDiffCoefCorrLabel, d_lab->d_stencilMatl,
		Task::OutOfDomain);
 
  tsk->computes(d_lab->d_scalNonLinSrcCorrLabel);

  sched->addTask(tsk, patches, matls);
}

      
//****************************************************************************
// Actually build linear matrix
//****************************************************************************
void ScalarSolver::buildLinearMatrixCorr(const ProcessorGroup* pc,
					 const PatchSubset* patches,
					 const MaterialSubset*,
					 DataWarehouse* old_dw,
					 DataWarehouse* new_dw,
					 int index)
{
  delt_vartype delT;
  old_dw->get(delT, d_lab->d_sharedState->get_delt_label() );
  double delta_t = delT;
#ifdef Runge_Kutta_3d
  double gamma_3 = 3.0/4.0;
  delta_t *= gamma_3;
#ifdef Runge_Kutta_3d_ssp
  delta_t /= gamma_3;
#endif
#endif

  
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int archIndex = 0; // only one arches material
    int matlIndex = d_lab->d_sharedState->
                    getArchesMaterial(archIndex)->getDWIndex(); 
    ArchesVariables scalarVars;
    
    // Get the PerPatch CellInformation data
    PerPatch<CellInformationP> cellInfoP;
    // checkpointing
    //  old_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    // new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    if (new_dw->exists(d_lab->d_cellInfoLabel, matlIndex, patch)) 
      new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    else {
      cellInfoP.setData(scinew CellInformation(patch));
      new_dw->put(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    }
    CellInformation* cellinfo = cellInfoP.get().get_rep();

    // from old_dw get PCELL, DENO, FO(index)
    new_dw->getCopy(scalarVars.cellType, d_lab->d_cellTypeLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);

    // from new_dw get DEN, VIS, F(index), U, V, W
  #ifdef Runge_Kutta_3d
    // old_density and old_scalar for Runge-Kutta are NOT from initial timestep
    // but from previous (Interm) time step
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityIntermLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.old_scalar, d_lab->d_scalarIntermLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.density, d_lab->d_densityIntermLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.viscosity, d_lab->d_viscosityIntermLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.scalar, d_lab->d_scalarIntermLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
  #else
  #ifdef Runge_Kutta_2nd
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.old_scalar, d_lab->d_scalarPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
  #else
    // ***warning* 21st July changed from IN to Pred
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.old_scalar, d_lab->d_scalarOUTBCLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
  #endif
    new_dw->getCopy(scalarVars.density, d_lab->d_densityPredLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.viscosity, d_lab->d_viscosityPredLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.scalar, d_lab->d_scalarPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
  #endif
    // for explicit get old values
  #ifdef Runge_Kutta_3d
    new_dw->getCopy(scalarVars.uVelocity, d_lab->d_uVelocityIntermLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.vVelocity, d_lab->d_vVelocityIntermLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.wVelocity, d_lab->d_wVelocityIntermLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  #else
    new_dw->getCopy(scalarVars.uVelocity, d_lab->d_uVelocityPredLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.vVelocity, d_lab->d_vVelocityPredLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.wVelocity, d_lab->d_wVelocityPredLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  #endif

  // allocate matrix coeffs
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++) {
      new_dw->allocate(scalarVars.scalarCoeff[ii], 
		       d_lab->d_scalCoefCorrLabel, ii, patch);
      new_dw->allocate(scalarVars.scalarConvectCoeff[ii],
		       d_lab->d_scalConvCoefCorrLabel, ii, patch);
      new_dw->allocate(scalarVars.scalarDiffusionCoeff[ii],
		       d_lab->d_scalDiffCoefCorrLabel, ii, patch);

    }
    new_dw->allocate(scalarVars.scalarLinearSrc, 
		     d_lab->d_scalLinSrcCorrLabel, matlIndex, patch);
    new_dw->allocate(scalarVars.scalarNonlinearSrc, 
		     d_lab->d_scalNonLinSrcCorrLabel, matlIndex, patch);
 
  // compute ith component of scalar stencil coefficients
  // inputs : scalarSP, [u,v,w]VelocityMS, densityCP, viscosityCTS
  // outputs: scalCoefSBLM
    d_discretize->calculateScalarCoeff(pc, patch,
				       delta_t, index, cellinfo, 
				       &scalarVars);

    // Calculate scalar source terms
    // inputs : [u,v,w]VelocityMS, scalarSP, densityCP, viscosityCTS
    // outputs: scalLinSrcSBLM, scalNonLinSrcSBLM
    d_source->calculateScalarSource(pc, patch,
				    delta_t, index, cellinfo, 
				    &scalarVars );

    // Calculate the scalar boundary conditions
    // inputs : scalarSP, scalCoefSBLM
    // outputs: scalCoefSBLM
    d_boundaryCondition->scalarBC(pc, patch,  index, cellinfo, 
				  &scalarVars);
  // apply multimaterial intrusion wallbc
    if (d_MAlab)
      d_boundaryCondition->mmscalarWallBC(pc, patch, cellinfo,
					  &scalarVars);

    // similar to mascal
    // inputs :
    // outputs:
    d_source->modifyScalarMassSource(pc, patch, delta_t, index, &scalarVars);
    
    // Calculate the scalar diagonal terms
    // inputs : scalCoefSBLM, scalLinSrcSBLM
    // outputs: scalCoefSBLM
    d_discretize->calculateScalarDiagonal(pc, patch, index, &scalarVars);
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++) {
      new_dw->put(scalarVars.scalarCoeff[ii], 
		  d_lab->d_scalCoefCorrLabel, ii, patch);
      new_dw->put(scalarVars.scalarDiffusionCoeff[ii],
		  d_lab->d_scalDiffCoefCorrLabel, ii, patch);

    }
    new_dw->put(scalarVars.scalarNonlinearSrc, 
		d_lab->d_scalNonLinSrcCorrLabel, matlIndex, patch);

  }
}


//****************************************************************************
// Schedule linear solve of scalar
//****************************************************************************
void
ScalarSolver::sched_scalarLinearSolveCorr(SchedulerP& sched,
					  const PatchSet* patches,
					  const MaterialSet* matls,
					  int index)
{
  Task* tsk = scinew Task("ScalarSolver::scalarLinearSolveCorr",
			  this,
			  &ScalarSolver::scalarLinearSolveCorr,
			  index);
  
  tsk->requires(Task::OldDW, d_lab->d_sharedState->get_delt_label());

  // coefficient for the variable for which solve is invoked
  //***warning changed in to pred  
  #ifdef Runge_Kutta_3d
  tsk->requires(Task::NewDW, d_lab->d_densityIntermLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalarIntermLabel, 
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  #ifndef Runge_Kutta_3d_ssp
  tsk->requires(Task::NewDW, d_lab->d_scalarTempLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  #endif
  #else
  tsk->requires(Task::NewDW, d_lab->d_densityPredLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalarPredLabel, 
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  #endif
  tsk->requires(Task::NewDW, d_lab->d_scalCoefCorrLabel, 
		d_lab->d_stencilMatl, Task::OutOfDomain,
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalNonLinSrcCorrLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  #ifdef Runge_Kutta_2nd
  tsk->requires(Task::NewDW, d_lab->d_scalarOUTBCLabel,
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_densityINLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  #endif
  #ifdef Runge_Kutta_3d_ssp
  tsk->requires(Task::NewDW, d_lab->d_scalarOUTBCLabel,
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_densityINLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  #endif
  tsk->computes(d_lab->d_scalarSPLabel);
  
  sched->addTask(tsk, patches, matls);
}
//****************************************************************************
// Actual scalar solve .. may be changed after recursive tasks are added
//****************************************************************************
void 
ScalarSolver::scalarLinearSolveCorr(const ProcessorGroup* pc,
				const PatchSubset* patches,
				const MaterialSubset*,
				DataWarehouse* old_dw,
				DataWarehouse* new_dw,
				int index)
{
  delt_vartype delT;
  old_dw->get(delT, d_lab->d_sharedState->get_delt_label() );
  double delta_t = delT;
#ifdef Runge_Kutta_3d
  double gamma_3 = 3.0/4.0;
  double zeta_2 = -5.0/12.0;
  delta_t *= gamma_3;
#ifdef Runge_Kutta_3d_ssp
  delta_t /= gamma_3;
#endif
#endif
  
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int archIndex = 0; // only one arches material
    int matlIndex = d_lab->d_sharedState->
                    getArchesMaterial(archIndex)->getDWIndex(); 
    ArchesVariables scalarVars;
    // Get the PerPatch CellInformation data
    PerPatch<CellInformationP> cellInfoP;
    // get old_dw from getTop function
    // checkpointing
    //  old_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    //  new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    //  old_dw->get(cellInfoP, d_cellInfoLabel, matlIndex, patch);
    if (new_dw->exists(d_lab->d_cellInfoLabel, matlIndex, patch)) 
      new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    else {
      cellInfoP.setData(scinew CellInformation(patch));
      new_dw->put(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    }
    CellInformation* cellinfo = cellInfoP.get().get_rep();
    // ***warning* 21st July changed from IN to Pred
  #ifdef Runge_Kutta_3d
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityIntermLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
  #else
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
  #endif
    // for explicit calculation
    {
    new_dw->allocate(scalarVars.scalar, d_lab->d_scalarSPLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
  #ifdef Runge_Kutta_3d
    new_dw->copyOut(scalarVars.scalar, d_lab->d_scalarIntermLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
  #else
    new_dw->copyOut(scalarVars.scalar, d_lab->d_scalarPredLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
  #endif
    }
    scalarVars.old_scalar.allocate(scalarVars.scalar.getLowIndex(),
				   scalarVars.scalar.getHighIndex());
    scalarVars.old_scalar.copy(scalarVars.scalar);
    
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++)
      new_dw->getCopy(scalarVars.scalarCoeff[ii], d_lab->d_scalCoefCorrLabel, 
		  ii, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.scalarNonlinearSrc, d_lab->d_scalNonLinSrcCorrLabel,
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->allocate(scalarVars.residualScalar, d_lab->d_scalarRes,
		     matlIndex, patch);
  // apply underelax to eqn
    d_linearSolver->computeScalarUnderrelax(pc, patch, index, 
					    &scalarVars);
    // make it a separate task later
    d_linearSolver->scalarLisolve(pc, patch, index, delta_t, 
				  &scalarVars, cellinfo, d_lab);
  #ifdef Runge_Kutta_3d
  #ifndef Runge_Kutta_3d_ssp
    constCCVariable<double> temp_scalar;
    constCCVariable<double> old_density;

    new_dw->get(old_density, d_lab->d_densityIntermLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->get(temp_scalar, d_lab->d_scalarTempLabel,
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    
    IntVector indexLow = patch->getCellFORTLowIndex();
    IntVector indexHigh = patch->getCellFORTHighIndex();
    
    for (int colZ = indexLow.z(); colZ <= indexHigh.z(); colZ ++) {
      for (int colY = indexLow.y(); colY <= indexHigh.y(); colY ++) {
        for (int colX = indexLow.x(); colX <= indexHigh.x(); colX ++) {

            IntVector currCell(colX, colY, colZ);

            scalarVars.scalar[currCell] += zeta_2*temp_scalar[currCell]/
            old_density[currCell];
            if (scalarVars.scalar[currCell] > 1.0) 
		scalarVars.scalar[currCell] = 1.0;
            else if (scalarVars.scalar[currCell] < 0.0)
            	scalarVars.scalar[currCell] = 0.0;
        }
      }
    }
  #endif
  #endif
  #ifdef Runge_Kutta_2nd
    constCCVariable<double> old_scalar;
    constCCVariable<double> old_density;
    constCCVariable<double> new_density;

    new_dw->get(old_scalar, d_lab->d_scalarOUTBCLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->get(old_density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->get(new_density, d_lab->d_densityPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    
    IntVector indexLow = patch->getCellFORTLowIndex();
    IntVector indexHigh = patch->getCellFORTHighIndex();
    
    for (int colZ = indexLow.z(); colZ <= indexHigh.z(); colZ ++) {
      for (int colY = indexLow.y(); colY <= indexHigh.y(); colY ++) {
        for (int colX = indexLow.x(); colX <= indexHigh.x(); colX ++) {

            IntVector currCell(colX, colY, colZ);

            scalarVars.scalar[currCell] = (scalarVars.scalar[currCell]+
            old_density[currCell]/new_density[currCell]*
	    old_scalar[currCell])/2.0;
            if (scalarVars.scalar[currCell] > 1.0) 
		scalarVars.scalar[currCell] = 1.0;
            else if (scalarVars.scalar[currCell] < 0.0)
            	scalarVars.scalar[currCell] = 0.0;
        }
      }
    }
  #endif
  #ifdef Runge_Kutta_3d_ssp
    constCCVariable<double> old_scalar;
    constCCVariable<double> old_density;
    constCCVariable<double> new_density;

    new_dw->get(old_scalar, d_lab->d_scalarOUTBCLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->get(old_density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->get(new_density, d_lab->d_densityIntermLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    
    IntVector indexLow = patch->getCellFORTLowIndex();
    IntVector indexHigh = patch->getCellFORTHighIndex();
    
    for (int colZ = indexLow.z(); colZ <= indexHigh.z(); colZ ++) {
      for (int colY = indexLow.y(); colY <= indexHigh.y(); colY ++) {
        for (int colX = indexLow.x(); colX <= indexHigh.x(); colX ++) {

            IntVector currCell(colX, colY, colZ);

            scalarVars.scalar[currCell] = (2.0*scalarVars.scalar[currCell]+
            old_density[currCell]/new_density[currCell]*
	    old_scalar[currCell])/3.0;
            if (scalarVars.scalar[currCell] > 1.0) 
		scalarVars.scalar[currCell] = 1.0;
            else if (scalarVars.scalar[currCell] < 0.0)
            	scalarVars.scalar[currCell] = 0.0;
        }
      }
    }
  #endif

  // put back the results
    new_dw->put(scalarVars.scalar, d_lab->d_scalarSPLabel, 
		matlIndex, patch);
  }
}

//****************************************************************************
// Schedule solve of linearized scalar equation, intermediate step
//****************************************************************************
void 
ScalarSolver::solveInterm(SchedulerP& sched,
			const PatchSet* patches,
			const MaterialSet* matls,
			int index)
{
  //computes stencil coefficients and source terms
  // requires : scalarIN, [u,v,w]VelocitySPBC, densityIN, viscosityIN
  // computes : scalCoefSBLM, scalLinSrcSBLM, scalNonLinSrcSBLM
  sched_buildLinearMatrixInterm(sched, patches, matls, index);
  
  // Schedule the scalar solve
  // require : scalarIN, scalCoefSBLM, scalNonLinSrcSBLM
  // compute : scalResidualSS, scalCoefSS, scalNonLinSrcSS, scalarSP
  //d_linearSolver->sched_scalarSolve(level, sched, new_dw, matrix_dw, index);
  sched_scalarLinearSolveInterm(sched, patches, matls, index);
}

//****************************************************************************
// Schedule build of linear matrix
//****************************************************************************
void 
ScalarSolver::sched_buildLinearMatrixInterm(SchedulerP& sched,
					  const PatchSet* patches,
					  const MaterialSet* matls,
					  int index)
{
  Task* tsk = scinew Task("ScalarSolver::BuildCoeffInterm",
			  this,
			  &ScalarSolver::buildLinearMatrixInterm,
			  index);


  tsk->requires(Task::OldDW, d_lab->d_sharedState->get_delt_label());
  
  // This task requires scalar and density from old time step for transient
  // calculation
  //DataWarehouseP old_dw = new_dw->getTop();
  tsk->requires(Task::NewDW, d_lab->d_cellTypeLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
//  tsk->requires(Task::NewDW, d_lab->d_scalarOUTBCLabel,
//		Ghost::AroundCells, Arches::ONEGHOSTCELL);
//  tsk->requires(Task::NewDW, d_lab->d_densityINLabel, 
//		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalarPredLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_densityPredLabel, 
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_viscosityPredLabel,
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_uVelocityPredLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_vVelocityPredLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_wVelocityPredLabel,
		Ghost::AroundFaces, Arches::ONEGHOSTCELL);

      // added one more argument of index to specify scalar component
  tsk->computes(d_lab->d_scalCoefIntermLabel, d_lab->d_stencilMatl,
		Task::OutOfDomain);
  tsk->computes(d_lab->d_scalDiffCoefIntermLabel, d_lab->d_stencilMatl,
		Task::OutOfDomain);
  tsk->computes(d_lab->d_scalNonLinSrcIntermLabel);

  sched->addTask(tsk, patches, matls);
}

      
//****************************************************************************
// Actually build linear matrix
//****************************************************************************
void ScalarSolver::buildLinearMatrixInterm(const ProcessorGroup* pc,
					 const PatchSubset* patches,
					 const MaterialSubset*,
					 DataWarehouse* old_dw,
					 DataWarehouse* new_dw,
					 int index)
{
  delt_vartype delT;
  old_dw->get(delT, d_lab->d_sharedState->get_delt_label() );
  double delta_t = delT;
  double gamma_2 = 5.0/12.0;
  delta_t *= gamma_2; 
  #ifdef Runge_Kutta_3d_ssp
  delta_t /= gamma_2; 
  #endif
  
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int archIndex = 0; // only one arches material
    int matlIndex = d_lab->d_sharedState->
                    getArchesMaterial(archIndex)->getDWIndex(); 
    ArchesVariables scalarVars;
    
    // Get the PerPatch CellInformation data
    PerPatch<CellInformationP> cellInfoP;
    // checkpointing
    //  old_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    // new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    if (new_dw->exists(d_lab->d_cellInfoLabel, matlIndex, patch)) 
      new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    else {
      cellInfoP.setData(scinew CellInformation(patch));
      new_dw->put(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    }
    CellInformation* cellinfo = cellInfoP.get().get_rep();

    // from old_dw get PCELL, DENO, FO(index)
    new_dw->getCopy(scalarVars.cellType, d_lab->d_cellTypeLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    // ***warning* 21st July changed from IN to Pred
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.old_scalar, d_lab->d_scalarPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);

    // from new_dw get DEN, VIS, F(index), U, V, W
    new_dw->getCopy(scalarVars.density, d_lab->d_densityPredLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.viscosity, d_lab->d_viscosityPredLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.scalar, d_lab->d_scalarPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    // for explicit get old values
    new_dw->getCopy(scalarVars.uVelocity, d_lab->d_uVelocityPredLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.vVelocity, d_lab->d_vVelocityPredLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);
    new_dw->getCopy(scalarVars.wVelocity, d_lab->d_wVelocityPredLabel, 
		matlIndex, patch, Ghost::AroundFaces, Arches::ONEGHOSTCELL);

  // allocate matrix coeffs
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++) {
      new_dw->allocate(scalarVars.scalarCoeff[ii], 
		       d_lab->d_scalCoefIntermLabel, ii, patch);
      new_dw->allocate(scalarVars.scalarConvectCoeff[ii],
		       d_lab->d_scalConvCoefIntermLabel, ii, patch);
      new_dw->allocate(scalarVars.scalarDiffusionCoeff[ii],
		       d_lab->d_scalDiffCoefIntermLabel, ii, patch);
    }
    new_dw->allocate(scalarVars.scalarLinearSrc, 
		     d_lab->d_scalLinSrcIntermLabel, matlIndex, patch);
    new_dw->allocate(scalarVars.scalarNonlinearSrc, 
		     d_lab->d_scalNonLinSrcIntermLabel, matlIndex, patch);
 
  // compute ith component of scalar stencil coefficients
  // inputs : scalarSP, [u,v,w]VelocityMS, densityCP, viscosityCTS
  // outputs: scalCoefSBLM
    d_discretize->calculateScalarCoeff(pc, patch,
				       delta_t, index, cellinfo, 
				       &scalarVars);

    // Calculate scalar source terms
    // inputs : [u,v,w]VelocityMS, scalarSP, densityCP, viscosityCTS
    // outputs: scalLinSrcSBLM, scalNonLinSrcSBLM
    d_source->calculateScalarSource(pc, patch,
				    delta_t, index, cellinfo, 
				    &scalarVars );

    // Calculate the scalar boundary conditions
    // inputs : scalarSP, scalCoefSBLM
    // outputs: scalCoefSBLM
    d_boundaryCondition->scalarBC(pc, patch,  index, cellinfo, 
				  &scalarVars);
  // apply multimaterial intrusion wallbc
    if (d_MAlab)
      d_boundaryCondition->mmscalarWallBC(pc, patch, cellinfo,
					  &scalarVars);

    // similar to mascal
    // inputs :
    // outputs:
    d_source->modifyScalarMassSource(pc, patch, delta_t, index, &scalarVars);
    
    // Calculate the scalar diagonal terms
    // inputs : scalCoefSBLM, scalLinSrcSBLM
    // outputs: scalCoefSBLM
    d_discretize->calculateScalarDiagonal(pc, patch, index, &scalarVars);
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++) {
      new_dw->put(scalarVars.scalarCoeff[ii], 
		  d_lab->d_scalCoefIntermLabel, ii, patch);
      new_dw->put(scalarVars.scalarDiffusionCoeff[ii],
		  d_lab->d_scalDiffCoefIntermLabel, ii, patch);
    }
    new_dw->put(scalarVars.scalarNonlinearSrc, 
		d_lab->d_scalNonLinSrcIntermLabel, matlIndex, patch);

  }
}


//****************************************************************************
// Schedule linear solve of scalar
//****************************************************************************
void
ScalarSolver::sched_scalarLinearSolveInterm(SchedulerP& sched,
					  const PatchSet* patches,
					  const MaterialSet* matls,
					  int index)
{
  Task* tsk = scinew Task("ScalarSolver::scalarLinearSolveInterm",
			  this,
			  &ScalarSolver::scalarLinearSolveInterm,
			  index);
  
  tsk->requires(Task::OldDW, d_lab->d_sharedState->get_delt_label());

  // coefficient for the variable for which solve is invoked
  //***warning changed in to pred  
  tsk->requires(Task::NewDW, d_lab->d_densityPredLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalarPredLabel, 
		Ghost::AroundCells, Arches::ONEGHOSTCELL);
  tsk->requires(Task::NewDW, d_lab->d_scalCoefIntermLabel, 
		d_lab->d_stencilMatl, Task::OutOfDomain,
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_scalNonLinSrcIntermLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  #ifndef Runge_Kutta_3d_ssp
  tsk->modifies(d_lab->d_scalarTempLabel);
  #endif
  #ifdef Runge_Kutta_3d_ssp
  tsk->requires(Task::NewDW, d_lab->d_scalarOUTBCLabel,
		Ghost::None, Arches::ZEROGHOSTCELLS);
  tsk->requires(Task::NewDW, d_lab->d_densityINLabel, 
		Ghost::None, Arches::ZEROGHOSTCELLS);
  #endif
  tsk->computes(d_lab->d_scalarIntermLabel);
  
  sched->addTask(tsk, patches, matls);
}
//****************************************************************************
// Actual scalar solve .. may be changed after recursive tasks are added
//****************************************************************************
void 
ScalarSolver::scalarLinearSolveInterm(const ProcessorGroup* pc,
				const PatchSubset* patches,
				const MaterialSubset*,
				DataWarehouse* old_dw,
				DataWarehouse* new_dw,
				int index)
{
  delt_vartype delT;
  old_dw->get(delT, d_lab->d_sharedState->get_delt_label() );
  double delta_t = delT;
  double gamma_2 = 5.0/12.0;
  double zeta_1 = -17.0/60.0;
  delta_t *= gamma_2; 
  #ifdef Runge_Kutta_3d_ssp
  delta_t /= gamma_2; 
  #endif
  
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int archIndex = 0; // only one arches material
    int matlIndex = d_lab->d_sharedState->
                    getArchesMaterial(archIndex)->getDWIndex(); 
    ArchesVariables scalarVars;
    // Get the PerPatch CellInformation data
    PerPatch<CellInformationP> cellInfoP;
    // get old_dw from getTop function
    // checkpointing
    //  old_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    //  new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    //  old_dw->get(cellInfoP, d_cellInfoLabel, matlIndex, patch);
    if (new_dw->exists(d_lab->d_cellInfoLabel, matlIndex, patch)) 
      new_dw->get(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    else {
      cellInfoP.setData(scinew CellInformation(patch));
      new_dw->put(cellInfoP, d_lab->d_cellInfoLabel, matlIndex, patch);
    }
    CellInformation* cellinfo = cellInfoP.get().get_rep();
    // ***warning* 21st July changed from IN to Pred
    new_dw->getCopy(scalarVars.old_density, d_lab->d_densityPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    // for explicit calculation
    {
    new_dw->allocate(scalarVars.scalar, d_lab->d_scalarIntermLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    new_dw->copyOut(scalarVars.scalar, d_lab->d_scalarPredLabel, 
		matlIndex, patch, Ghost::AroundCells, Arches::ONEGHOSTCELL);
    }
    scalarVars.old_scalar.allocate(scalarVars.scalar.getLowIndex(),
				   scalarVars.scalar.getHighIndex());
    scalarVars.old_scalar.copy(scalarVars.scalar);
    
    for (int ii = 0; ii < d_lab->d_stencilMatl->size(); ii++)
      new_dw->getCopy(scalarVars.scalarCoeff[ii], d_lab->d_scalCoefIntermLabel, 
		  ii, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getCopy(scalarVars.scalarNonlinearSrc, d_lab->d_scalNonLinSrcIntermLabel,
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->allocate(scalarVars.residualScalar, d_lab->d_scalarRes,
		     matlIndex, patch);

  // apply underelax to eqn
    d_linearSolver->computeScalarUnderrelax(pc, patch, index, 
					    &scalarVars);
    // make it a separate task later
    d_linearSolver->scalarLisolve(pc, patch, index, delta_t, 
				  &scalarVars, cellinfo, d_lab);

  #ifndef Runge_Kutta_3d_ssp
    CCVariable<double> temp_scalar;
    constCCVariable<double> old_density;

    new_dw->get(old_density, d_lab->d_densityPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->getModifiable(temp_scalar, d_lab->d_scalarTempLabel,
                matlIndex, patch);
    
    IntVector indexLow = patch->getCellFORTLowIndex();
    IntVector indexHigh = patch->getCellFORTHighIndex();
    
    for (int colZ = indexLow.z(); colZ <= indexHigh.z(); colZ ++) {
      for (int colY = indexLow.y(); colY <= indexHigh.y(); colY ++) {
        for (int colX = indexLow.x(); colX <= indexHigh.x(); colX ++) {

            IntVector currCell(colX, colY, colZ);

            scalarVars.scalar[currCell] += zeta_1*temp_scalar[currCell]/
            old_density[currCell];
            temp_scalar[currCell] = old_density[currCell]*
	    (scalarVars.scalar[currCell]-
            scalarVars.old_scalar[currCell])/
            gamma_2-zeta_1*temp_scalar[currCell]/gamma_2;
            if (scalarVars.scalar[currCell] > 1.0) 
		scalarVars.scalar[currCell] = 1.0;
            else if (scalarVars.scalar[currCell] < 0.0)
            	scalarVars.scalar[currCell] = 0.0;
        }
      }
    } 
//    new_dw->put(temp_scalar, d_lab->d_scalarTempLabel, matlIndex, patch);
  #endif
  #ifdef Runge_Kutta_3d_ssp
    constCCVariable<double> old_scalar;
    constCCVariable<double> old_density;
    constCCVariable<double> new_density;

    new_dw->get(old_scalar, d_lab->d_scalarOUTBCLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->get(old_density, d_lab->d_densityINLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    new_dw->get(new_density, d_lab->d_densityPredLabel, 
		matlIndex, patch, Ghost::None, Arches::ZEROGHOSTCELLS);
    
    IntVector indexLow = patch->getCellFORTLowIndex();
    IntVector indexHigh = patch->getCellFORTHighIndex();
    
    for (int colZ = indexLow.z(); colZ <= indexHigh.z(); colZ ++) {
      for (int colY = indexLow.y(); colY <= indexHigh.y(); colY ++) {
        for (int colX = indexLow.x(); colX <= indexHigh.x(); colX ++) {

            IntVector currCell(colX, colY, colZ);

            scalarVars.scalar[currCell] = (scalarVars.scalar[currCell]+
            old_density[currCell]/new_density[currCell]*
	    3.0*old_scalar[currCell])/4.0;
            if (scalarVars.scalar[currCell] > 1.0) 
		scalarVars.scalar[currCell] = 1.0;
            else if (scalarVars.scalar[currCell] < 0.0)
            	scalarVars.scalar[currCell] = 0.0;
        }
      }
    }
  #endif
  
    // put back the results
    new_dw->put(scalarVars.scalar, d_lab->d_scalarIntermLabel, 
		matlIndex, patch);
  }
}
