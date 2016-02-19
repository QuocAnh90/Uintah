/*
 * The MIT License
 *
 * Copyright (c) 1997-2016 The University of Utah
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

#ifndef CCA_COMPONENTS_SCHEDULERS_THREADFUNNELEDSCHEDULER_H
#define CCA_COMPONENTS_SCHEDULERS_THREADFUNNELEDSCHEDULER_H

#include <CCA/Components/Schedulers/MPIScheduler.h>
#include <Core/Lockfree/Lockfree_Timer.hpp>
#include <Core/Malloc/AllocatorTags.hpp>

#include <mutex>
#include <queue>
#include <thread>

namespace Uintah {

class Task;
class DetailedTask;
class TaskRunner;

using clock_type = std::chrono::high_resolution_clock;
using nanoseconds = std::chrono::nanoseconds;

using TaskPool = Lockfree::CircularPool<   DetailedTask*
                                         , Lockfree::ENABLE_SIZE         // size model
                                         , Lockfree::EXCLUSIVE_INSTANCE // usage model
                                         , Uintah::MallocAllocator      // allocator
                                         , Uintah::MallocAllocator      // size_type allocator
                                       >;


/**************************************

 CLASS
 ThreadFunneledScheduler


 GENERAL INFORMATION
 ThreadFunneledScheduler.h

 Alan Humphrey
 Scientific Computing and Imaging Institute
 University of Utah

 KEYWORDS
 Task Scheduler, Multi-threaded MPI, CPU

 DESCRIPTION
 A multi-threaded scheduler that uses a combination of MPI + Pthreads

 Uses MPI_THREAD_FUNNELED, where ONLY the thread that called MPI_Init_thread will make MPI calls.

 Dynamic scheduling with non-deterministic, out-of-order execution of tasks at runtime.

 One MPI rank per numa region.

 Pthreads are pinned to individual CPU cores where these tasks are executed.

 Uses a decentralized model wherein all threads can access task queues,
 processes there own MPI send and recvs, with shared access to the DataWarehouse.


 WARNING
 This scheduler is still EXPERIMENTAL and undergoing extensive
 development, not all tasks/components are thread-safe yet.

 Requires MPI_THREAD_Funnled (level-1) support.

 MPI_THREAD_FUNNELED: Only the thread that called MPI_Init_thread makes MPI calls.

 ****************************************/

class ThreadFunneledScheduler : public MPIScheduler {

  public:

    ThreadFunneledScheduler( const ProcessorGroup* myworld, const Output* oport, ThreadFunneledScheduler* parentScheduler = 0 )
      : MPIScheduler(myworld, oport, parentScheduler)
      {}

    virtual ~ThreadFunneledScheduler();

    virtual void problemSetup( const ProblemSpecP& prob_spec, SimulationStateP& state );

    virtual SchedulerP createSubScheduler();

    virtual void execute( int tgnum = 0, int iteration = 0);

    virtual bool useInternalDeps() { return !d_sharedState->isCopyDataTimestep(); }

    friend class TaskRunner;


  private:

    // eliminate copy, assignment and move
    ThreadFunneledScheduler( const ThreadFunneledScheduler & )            = delete;
    ThreadFunneledScheduler& operator=( const ThreadFunneledScheduler & ) = delete;
    ThreadFunneledScheduler( ThreadFunneledScheduler && )                 = delete;
    ThreadFunneledScheduler& operator=( ThreadFunneledScheduler && )      = delete;

    DetailedTask* get_reduction_task();

    void process_mpi( int iteration );

    static void init_threads( ThreadFunneledScheduler * scheduler, int num_threads );

    void run_task( DetailedTask* task, int iteration );

    void select_tasks();

    void thread_fence();

    static void set_runner( TaskRunner*, int tid );


    TaskPool   m_task_pool{};
    TaskPool   m_mpi_pending_pool{};
    TaskPool   m_mpi_test_pool{};

    std::queue<DetailedTask*>                 m_reduction_tasks{};
    std::vector<std::vector<DetailedTask*> >  m_phase_task_list{};

    Lockfree::Timer  m_mpi_test_time{};

    // thread shared data, needs lock protection when accessed
    std::vector<int>           m_phase_tasks{};
    std::vector<int>           m_phase_tasks_done{};
//    std::vector<DetailedTask*> m_phase_sync_tasks{};
    DetailedTasks*             m_detailed_tasks{};

    bool     m_abort{ false };
    int      m_current_iteration{ 0 };
    int      m_num_tasks_done{ 0 };
    int      m_num_tasks{ 0 };
    int      m_current_phase{ 0 };
    int      m_num_phases{ 0 };
    int      m_abort_point{ 0 };
    int      m_num_threads{ 0 };

};


class TaskRunner {

  public:

    TaskRunner() = default;

    TaskRunner( ThreadFunneledScheduler* scheduler)
      : m_scheduler{ scheduler }
      , m_task_wait_time{}
      , m_task_exec_time{}
    {}

    ~TaskRunner() {};

    void run() const;


  private:

    ThreadFunneledScheduler*  m_scheduler{ nullptr };
    Lockfree::Timer           m_task_wait_time{};
    Lockfree::Timer           m_task_exec_time{};


    TaskRunner( const TaskRunner & )            = delete;
    TaskRunner& operator=( const TaskRunner & ) = delete;
    TaskRunner( TaskRunner &&)                  = delete;
    TaskRunner& operator=( TaskRunner && )      = delete;

    friend class ThreadFunneledScheduler;

}; //class ThreadFunneledScheduler

} // namespace Uintah

#endif  // CCA_COMPONENTS_SCHEDULERS_THREADFUNNELEDSCHEDULER_H
