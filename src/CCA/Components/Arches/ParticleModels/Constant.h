#ifndef Uintah_Component_Arches_Constant_h
#define Uintah_Component_Arches_Constant_h

#include <CCA/Components/Arches/Task/TaskInterface.h>
#include <CCA/Components/Arches/ParticleModels/ParticleTools.h>

//-------------------------------------------------------

/**
 * @class    Constant
 * @author   Alex Abboud
 * @date     October 2014
 *
 * @brief    This class sets a constant source term for particles
 *
 * @details  A constant source term for easier debugging
 *
 */

//-------------------------------------------------------

namespace Uintah{

  //T is the dependent variable type
  template <typename T>
  class Constant : public TaskInterface {

  public:

    Constant<T>( std::string task_name, int matl_index, const std::string var_name, const int N );
    ~Constant<T>(){};

    TaskAssignedExecutionSpace loadTaskComputeBCsFunctionPointers();

    TaskAssignedExecutionSpace loadTaskInitializeFunctionPointers();

    TaskAssignedExecutionSpace loadTaskEvalFunctionPointers();

    TaskAssignedExecutionSpace loadTaskTimestepInitFunctionPointers();

    TaskAssignedExecutionSpace loadTaskRestartInitFunctionPointers();

    void problemSetup( ProblemSpecP& db );

    void create_local_labels();

    class Builder : public TaskInterface::TaskBuilder {

    public:

      Builder( std::string task_name, int matl_index, std::string base_var_name, const int N ) :
      m_task_name(task_name), m_matl_index(matl_index), _base_var_name(base_var_name), _N(N){}
      ~Builder(){}

      Constant* build()
      { return scinew Constant<T>( m_task_name, m_matl_index, _base_var_name, _N ); }

    private:

      std::string m_task_name;
      int m_matl_index;
      std::string _base_var_name;
      const int _N;

    };

  protected:

    void register_initialize( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry , const bool packed_tasks);

    void register_timestep_init( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry , const bool packed_tasks);

    void register_timestep_eval( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry, const int time_substep , const bool packed_tasks){}

    void register_compute_bcs( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry, const int time_substep , const bool packed_tasks){};

    template <typename ExecutionSpace, typename MemSpace>
    void compute_bcs( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace, MemSpace>& executionObject ){}

    template <typename ExecutionSpace, typename MemSpace>
    void initialize( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace, MemSpace>& executionObject );

    template<typename ExecutionSpace, typename MemSpace> void timestep_init( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace,MemSpace>& exObj);

    template <typename ExecutionSpace, typename MemSpace>
    void eval( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace, MemSpace>& executionObject ){}

  private:

    const std::string _base_var_name;

    const int _N;                      //<<< The number of "environments"
    std::vector<double> _const;        //<<< constant source value/environment

    /** @brief Set the actual value of the constant to the grid variable **/ 
    void set_value( const Patch* patch, ArchesTaskInfoManager* tsk_info );

  };

  //Function definitions:
  //------------------------------------------------------------------------------------------------
  template <typename T>
  void Constant<T>::create_local_labels(){
    for ( int i = 0; i < _N; i++ ){

      std::string name = ArchesCore::append_env( _base_var_name, i ); 
      register_new_variable<T>(name);

      name = ArchesCore::append_qn_env( _base_var_name, i ); 
      register_new_variable<T>(name); 

    }
  }

  //------------------------------------------------------------------------------------------------
  template <typename T>
  Constant<T>::Constant( std::string task_name, int matl_index,
                         const std::string base_var_name, const int N ) :
  TaskInterface( task_name, matl_index ), _base_var_name(base_var_name), _N(N){
  }

  //--------------------------------------------------------------------------------------------------
  template <typename T>
  TaskAssignedExecutionSpace Constant<T>::loadTaskComputeBCsFunctionPointers()
  {
    return create_portable_arches_tasks<TaskInterface::BC>( this
                                       , &Constant<T>::compute_bcs<UINTAH_CPU_TAG>     // Task supports non-Kokkos builds
                                       //, &Constant<T>::compute_bcs<KOKKOS_OPENMP_TAG>  // Task supports Kokkos::OpenMP builds
                                       //, &Constant<T>::compute_bcs<KOKKOS_CUDA_TAG>    // Task supports Kokkos::Cuda builds
                                       );
  }

  //--------------------------------------------------------------------------------------------------
  template <typename T>
  TaskAssignedExecutionSpace Constant<T>::loadTaskInitializeFunctionPointers()
  {
    return create_portable_arches_tasks<TaskInterface::INITIALIZE>( this
                                       , &Constant<T>::initialize<UINTAH_CPU_TAG>     // Task supports non-Kokkos builds
                                       //, &Constant<T>::initialize<KOKKOS_OPENMP_TAG>  // Task supports Kokkos::OpenMP builds
                                       //, &Constant<T>::initialize<KOKKOS_CUDA_TAG>    // Task supports Kokkos::Cuda builds
                                       );
  }

  //--------------------------------------------------------------------------------------------------
  template <typename T>
  TaskAssignedExecutionSpace Constant<T>::loadTaskEvalFunctionPointers()
  {
    return create_portable_arches_tasks<TaskInterface::TIMESTEP_EVAL>( this
                                       , &Constant<T>::eval<UINTAH_CPU_TAG>     // Task supports non-Kokkos builds
                                       //, &Constant<T>::eval<KOKKOS_OPENMP_TAG>  // Task supports Kokkos::OpenMP builds
                                       //, &Constant<T>::eval<KOKKOS_CUDA_TAG>    // Task supports Kokkos::Cuda builds
                                       );
  }

  //--------------------------------------------------------------------------------------------------
  template <typename T>
  TaskAssignedExecutionSpace Constant<T>::loadTaskTimestepInitFunctionPointers()
  {
    return create_portable_arches_tasks<TaskInterface::TIMESTEP_INITIALIZE>( this
                                       , &Constant<T>::timestep_init<UINTAH_CPU_TAG>     // Task supports non-Kokkos builds
                                       , &Constant<T>::timestep_init<KOKKOS_OPENMP_TAG>  // Task supports Kokkos::OpenMP builds
                                       );
  }

  //--------------------------------------------------------------------------------------------------
  template <typename T>
  TaskAssignedExecutionSpace Constant<T>::loadTaskRestartInitFunctionPointers()
  {
    return  TaskAssignedExecutionSpace::NONE_EXECUTION_SPACE;
  }

  //------------------------------------------------------------------------------------------------
  template <typename T>
  void Constant<T>::problemSetup( ProblemSpecP& db ){

    db->require("constant",_const);

  }

  //------------------------------------------------------------------------------------------------
  template <typename T>
  void Constant<T>::register_initialize( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry , const bool packed_tasks){

    for ( int ei = 0; ei < _N; ei++ ){

      std::string name = ArchesCore::append_env(_base_var_name, ei);
      register_variable( name, ArchesFieldContainer::COMPUTES, 0, ArchesFieldContainer::NEWDW, variable_registry );
      name = ArchesCore::append_qn_env(_base_var_name, ei);
      register_variable( name, ArchesFieldContainer::COMPUTES, 0, ArchesFieldContainer::NEWDW, variable_registry );

    }
  }

  //------------------------------------------------------------------------------------------------
  template <typename T>
  template<typename ExecutionSpace, typename MemSpace>
  void Constant<T>::initialize( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace, MemSpace>& executionObject ){

    set_value( patch, tsk_info ); 

  }

  //------------------------------------------------------------------------------------------------
  template <typename T>
  void Constant<T>::register_timestep_init( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry , const bool packed_tasks){

    for ( int ei = 0; ei < _N; ei++ ){

      std::string name = ArchesCore::append_env(_base_var_name, ei);
      register_variable( name, ArchesFieldContainer::COMPUTES, 0, ArchesFieldContainer::NEWDW, variable_registry );
      name = ArchesCore::append_qn_env(_base_var_name, ei);
      register_variable( name, ArchesFieldContainer::COMPUTES, 0, ArchesFieldContainer::NEWDW, variable_registry );

    }
  }

  template <typename T>
  template<typename ExecutionSpace, typename MemSpace> void
  Constant<T>::timestep_init( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace, MemSpace>& executionObject ){

    set_value( patch, tsk_info ); 

  }

  template <typename T> 
  void Constant<T>::set_value( const Patch* patch, ArchesTaskInfoManager* tsk_info ){ 

    for ( int ei = 0; ei < _N; ei++ ){

      std::string name = ArchesCore::append_env( _base_var_name, ei );
      T& model_value = tsk_info->get_uintah_field_add<T>(name);
      model_value.initialize(_const[ei]); 

      name = ArchesCore::append_qn_env( _base_var_name, ei );
      T& model_qn_value = tsk_info->get_uintah_field_add<T>(name);
      model_qn_value.initialize(_const[ei]); 

    }

  }

}
#endif
