#ifndef Uintah_Component_Arches_VelRhoHatBC_h
#define Uintah_Component_Arches_VelRhoHatBC_h

#include <CCA/Components/Arches/Task/AtomicTaskInterface.h>

//===============================================================

/**
* @class  VelRhoHatBC
* @author Jeremy Thornock
* @date   2016
*
* @brief Computes the BC on the hatted vel*rho
*
**/

//===============================================================

namespace Uintah{

  class VelRhoHatBC : public AtomicTaskInterface {

public:

    /** @brief Default constructor **/
    VelRhoHatBC( std::string task_name, int matl_index );

    /** @brief Default destructor **/
    ~VelRhoHatBC();

    TaskAssignedExecutionSpace loadTaskComputeBCsFunctionPointers();

    TaskAssignedExecutionSpace loadTaskInitializeFunctionPointers();

    TaskAssignedExecutionSpace loadTaskEvalFunctionPointers();

    TaskAssignedExecutionSpace loadTaskRestartInitFunctionPointers();
  
    TaskAssignedExecutionSpace loadTaskTimestepInitFunctionPointers();

    /** @brief Input file interface **/
    void problemSetup( ProblemSpecP& db );

    /** @brief Create local labels for the task **/
    void create_local_labels();

    /** @brief Registers all variables with pertinent information for the
     *         uintah dw interface **/
    void register_timestep_eval( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry,
                                const int time_substep, const bool pack_tasks );

    template <typename ExecutionSpace, typename MemSpace>
    void compute_bcs( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace, MemSpace>& executionObject ){}

    template <typename ExecutionSpace, typename MemSpace>
    void initialize( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace, MemSpace>& executionObject ){}

    template <typename ExecutionSpace, typename MemSpace>
    void eval( const Patch* patch, ArchesTaskInfoManager* tsk_info, ExecutionObject<ExecutionSpace, MemSpace>& executionObject );

    /** @brief Builder class containing instructions on how to build the task **/
    class Builder : public AtomicTaskInterface::AtomicTaskBuilder {

      public:

      Builder( std::string name, int matl_index ) : m_task_name(name), m_matl_index(matl_index){}
      ~Builder(){}

      VelRhoHatBC* build()
      { return scinew VelRhoHatBC( m_task_name, m_matl_index ); }

      private:

        std::string m_task_name;
        int m_matl_index;

    };

private:

    std::string m_xmom;
    std::string m_ymom;
    std::string m_zmom;


  };
} // namespace Uintah

#endif
