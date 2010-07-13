#include <CCA/Components/Arches/CoalModels/ParticleVelocity.h>
#include <CCA/Components/Arches/TransportEqns/EqnFactory.h>
#include <CCA/Components/Arches/TransportEqns/EqnBase.h>
#include <CCA/Components/Arches/TransportEqns/DQMOMEqn.h>
#include <CCA/Components/Arches/ArchesLabel.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Parallel/Parallel.h>

//===========================================================================

using namespace std;
using namespace Uintah; 

ParticleVelocity::ParticleVelocity( std::string modelName, 
                                    SimulationStateP& sharedState,
                                    const ArchesLabel* fieldLabels,
                                    vector<std::string> icLabelNames, 
                                    vector<std::string> scalarLabelNames,
                                    int qn ) 
: ModelBase(modelName, sharedState, fieldLabels, icLabelNames, scalarLabelNames, qn)
{
  d_quadNode = qn;

  // Create a label for this model
  d_modelLabel = VarLabel::create( modelName, CCVariable<Vector>::getTypeDescription() );

  // Create the gas phase source term associated with this model
  std::string gasSourceName = modelName + "_gasSource";
  d_gasLabel = VarLabel::create( gasSourceName, CCVariable<Vector>::getTypeDescription() );

  // Create velocity vector label (to store velocity components into a vector)
  std::string qnode;
  std::stringstream out;
  out << d_quadNode;
  qnode = out.str();
  std::string velname = "vel_qn";
  d_velocity_label = VarLabel::create( velname+qnode, CCVariable<Vector>::getTypeDescription() );
}

ParticleVelocity::~ParticleVelocity()
{
  delete d_boundaryCond;
}

//---------------------------------------------------------------------------
// Method: Problem Setup
//---------------------------------------------------------------------------
  void 
ParticleVelocity::problemSetup(const ProblemSpecP& params)
{
  // This method is called by child class problemSetup()'s

  ProblemSpecP db = params; 

  const ProblemSpecP db_root = db->getRootNode();
  ProblemSpecP db_physicalConstants = db_root->findBlock("PhysicalConstants");
  db_physicalConstants->require("viscosity", d_visc );
  db_physicalConstants->require("gravity", d_gravity );

  // set model clipping (not used yet...)
  //db->getWithDefault( "low_clip",  d_lowModelClip,  1.0e-6 );
  //db->getWithDefault( "high_clip", d_highModelClip, 999999 );

  db->getWithDefault( "partvelBC_eq_gasvelBC", d_gasBC, false ); 
  if(d_gasBC) {
    proc0cout << endl << " WARNING: Arches: ParticleVelocity: Setting particle velocities equal to gas velocities at the boundary using the <partvelBC_eq_gasvelBC> tag will cause errors!" << endl << endl;
  }

  // grab weight scaling factor and small value
  DQMOMEqnFactory& dqmom_eqn_factory = DQMOMEqnFactory::self();

  std::string temp_weight_name = "w_qn";
  std::string node;
  std::stringstream out;
  out << d_quadNode;
  node = out.str();
  temp_weight_name += node;
  EqnBase& t_weight_eqn = dqmom_eqn_factory.retrieve_scalar_eqn( temp_weight_name );
  DQMOMEqn& weight_eqn = dynamic_cast<DQMOMEqn&>(t_weight_eqn);

  d_w_small = weight_eqn.getSmallClip();
  d_w_scaling_factor = weight_eqn.getScalingConstant();
  d_weight_label = weight_eqn.getTransportEqnLabel();

  if( db->findBlock("scaling_const") ) {
    // several classes use particle velocity without dealing with scaling... This should be changed, but for now don't allow the user to scale particle velocity to prevent inadvertent errors
    string err = "ERROR: Arches: ParticleVelocity: You specified a scaling constant for the particle velocity internal coordinate.  This is not allowed.";
    throw ProblemSetupException(err,__FILE__,__LINE__);
  }

  d_boundaryCond = scinew BoundaryCondition_new( d_fieldLabels );

}


//-------------------------------------------------------------------------
// Method: Actually do the dummy initialization
//-------------------------------------------------------------------------
/** @details
 This is called from ExplicitSolver::noSolve(), which skips the first timestep
 so that the initial conditions are correct.

This method was originally in ModelBase, but it requires creating CCVariables
 for the model and gas source terms, and the CCVariable type (double, Vector, &c.)
 is model-dependent.  Putting the method here eliminates if statements in 
 ModelBase and keeps the ModelBase class as generic as possible.

@see ExplicitSolver::noSolve()
 */
void
ParticleVelocity::dummyInit( const ProcessorGroup* pc,
                             const PatchSubset* patches, 
                             const MaterialSubset* matls, 
                             DataWarehouse* old_dw, 
                             DataWarehouse* new_dw )
{
  for( int p=0; p < patches->size(); ++p ) {

    Ghost::GhostType  gn = Ghost::None;

    const Patch* patch = patches->get(p);
    int archIndex = 0;
    int matlIndex = d_fieldLabels->d_sharedState->getArchesMaterial(archIndex)->getDWIndex(); 

    CCVariable<Vector> ModelTerm;
    CCVariable<Vector> GasModelTerm;
    
    constCCVariable<Vector> oldModelTerm;
    constCCVariable<Vector> oldGasModelTerm;

    new_dw->allocateAndPut( ModelTerm,    d_modelLabel, matlIndex, patch );
    new_dw->allocateAndPut( GasModelTerm, d_gasLabel,   matlIndex, patch ); 

    old_dw->get( oldModelTerm,    d_modelLabel, matlIndex, patch, gn, 0 );
    old_dw->get( oldGasModelTerm, d_gasLabel,   matlIndex, patch, gn, 0 );
    
    ModelTerm.copyData(oldModelTerm);
    GasModelTerm.copyData(oldGasModelTerm);
  }
}

//---------------------------------------------------------------------------
// Method: Schedule the initialization of special variables unique to model
//---------------------------------------------------------------------------
void 
ParticleVelocity::sched_initVars( const LevelP& level, SchedulerP& sched )
{
  std::string taskname = "ParticleVelocity::initVars";
  Task* tsk = scinew Task(taskname, this, &ParticleVelocity::initVars);

  tsk->computes( d_modelLabel );
  tsk->computes( d_gasLabel   );
  tsk->computes( d_velocity_label );

  sched->addTask(tsk, level->eachPatch(), d_sharedState->allArchesMaterials()); 
}

//-------------------------------------------------------------------------
// Method: Initialize special variables unique to the model
//-------------------------------------------------------------------------
void
ParticleVelocity::initVars( const ProcessorGroup * pc, 
                            const PatchSubset    * patches, 
                            const MaterialSubset * matls, 
                            DataWarehouse        * old_dw, 
                            DataWarehouse        * new_dw )
{
  for( int p=0; p < patches->size(); p++ ) {

    const Patch* patch = patches->get(p);
    int archIndex = 0;
    int matlIndex = d_fieldLabels->d_sharedState->getArchesMaterial(archIndex)->getDWIndex(); 

    CCVariable<Vector> model_value; 
    new_dw->allocateAndPut( model_value, d_modelLabel, matlIndex, patch ); 
    model_value.initialize( Vector(0.0,0.0,0.0) );

    CCVariable<Vector> gas_value; 
    new_dw->allocateAndPut( gas_value, d_gasLabel, matlIndex, patch ); 
    gas_value.initialize( Vector(0.0,0.0,0.0) );

    CCVariable<Vector> particle_velocity; 
    new_dw->allocateAndPut( particle_velocity, d_velocity_label, matlIndex, patch ); 
    particle_velocity.initialize( Vector(0.0,0.0,0.0) );

  }
}


