#include <CCA/Components/Arches/Utility/BoundaryInfo.h>

using namespace Uintah;

BoundaryInfo::BoundaryInfo( std::string task_name, int matl_index ) :
TaskInterface( task_name, matl_index ) {
}

BoundaryInfo::~BoundaryInfo(){
}

//--------------------------------------------------------------------------------------------------
TaskAssignedExecutionSpace BoundaryInfo::loadTaskComputeBCsFunctionPointers()
{
  return create_portable_arches_tasks( this
                                     , TaskInterface::BC
                                     , &BoundaryInfo::compute_bcs<UINTAH_CPU_TAG>     // Task supports non-Kokkos builds
                                     //, &BoundaryInfo::compute_bcs<KOKKOS_OPENMP_TAG>  // Task supports Kokkos::OpenMP builds
                                     //, &BoundaryInfo::compute_bcs<KOKKOS_CUDA_TAG>    // Task supports Kokkos::Cuda builds
                                     );
}

//--------------------------------------------------------------------------------------------------
TaskAssignedExecutionSpace BoundaryInfo::loadTaskInitializeFunctionPointers()
{
  return create_portable_arches_tasks( this
                                     , TaskInterface::INITIALIZE
                                     , &BoundaryInfo::initialize<UINTAH_CPU_TAG>     // Task supports non-Kokkos builds
                                     //, &BoundaryInfo::initialize<KOKKOS_OPENMP_TAG>  // Task supports Kokkos::OpenMP builds
                                     //, &BoundaryInfo::initialize<KOKKOS_CUDA_TAG>    // Task supports Kokkos::Cuda builds
                                     );
}

//--------------------------------------------------------------------------------------------------
TaskAssignedExecutionSpace BoundaryInfo::loadTaskEvalFunctionPointers()
{
  return create_portable_arches_tasks( this
                                     , TaskInterface::TIMESTEP_EVAL
                                     , &BoundaryInfo::eval<UINTAH_CPU_TAG>     // Task supports non-Kokkos builds
                                     //, &BoundaryInfo::eval<KOKKOS_OPENMP_TAG>  // Task supports Kokkos::OpenMP builds
                                     //, &BoundaryInfo::eval<KOKKOS_CUDA_TAG>    // Task supports Kokkos::Cuda builds
                                     );
}


void
BoundaryInfo::problemSetup( ProblemSpecP& db ){
}

void
BoundaryInfo::create_local_labels(){

  register_new_variable<SFCXVariable<double> >( "area_fraction_x" );
  register_new_variable<SFCYVariable<double> >( "area_fraction_y" );
  register_new_variable<SFCZVariable<double> >( "area_fraction_z" );

}

//
//------------------------------------------------
//-------------- INITIALIZATION ------------------
//------------------------------------------------
//

typedef std::vector<ArchesFieldContainer::VariableInformation> VarInfoVecT;

void
BoundaryInfo::register_initialize( VarInfoVecT& variable_registry , const bool packed_tasks){

  //FUNCITON CALL     STRING NAME(VL)     TYPE       DEPENDENCY    GHOST DW     VR
  register_variable( "area_fraction_x" , ArchesFieldContainer::COMPUTES , variable_registry );
  register_variable( "area_fraction_y" , ArchesFieldContainer::COMPUTES , variable_registry );
  register_variable( "area_fraction_z" , ArchesFieldContainer::COMPUTES , variable_registry );

}

template<typename ExecutionSpace, typename MemorySpace>
void BoundaryInfo::initialize( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace, MemorySpace>& executionObject ){}

//
//------------------------------------------------
//------------- TIMESTEP INIT --------------------
//------------------------------------------------
//
void
BoundaryInfo::register_timestep_init( VarInfoVecT& variable_registry , const bool packed_tasks){

  register_variable( "area_fraction_x", ArchesFieldContainer::COMPUTES, variable_registry );
  register_variable( "area_fraction_y", ArchesFieldContainer::COMPUTES, variable_registry );
  register_variable( "area_fraction_z", ArchesFieldContainer::COMPUTES, variable_registry );

  register_variable( "area_fraction_x", ArchesFieldContainer::REQUIRES, 0, ArchesFieldContainer::OLDDW, variable_registry );
  register_variable( "area_fraction_y", ArchesFieldContainer::REQUIRES, 0, ArchesFieldContainer::OLDDW, variable_registry );
  register_variable( "area_fraction_z", ArchesFieldContainer::REQUIRES, 0, ArchesFieldContainer::OLDDW, variable_registry );

}

void
BoundaryInfo::timestep_init( const Patch* patch, ArchesTaskInfoManager* tsk_info ){}
//
//------------------------------------------------
//------------- TIMESTEP WORK --------------------
//------------------------------------------------
//

void
BoundaryInfo::register_timestep_eval( VarInfoVecT& variable_registry,
                                      const int time_substep, const bool packed_tasks ){

}

template<typename ExecutionSpace, typename MemorySpace>
void BoundaryInfo::eval( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace, MemorySpace>& executionObject ){}
