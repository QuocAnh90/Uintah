/*
 * This file was automatically generated by SCC - do NOT edit!
 * You should edit ConditionVariable.scc instead 
 */

#ifndef SCI_THREAD_CONDITIONVARIABLE_H
#define SCI_THREAD_CONDITIONVARIABLE_H 1

/**
 * Condition variable primitive.  When a thread calls the
 * <i>wait</i> method,which will block until another thread calls
 * the <i>cond_signal</i> or <i>cond_broadcast</i> methods.  When
 * there are multiple threads waiting, <i>cond_broadcast</i> will unblock
 * all of them, while <i>cond_signal</i> will unblock only one (an
 * arbitrary one) of them.  This primitive is used to allow a thread
 * to wait for some condition to exist, such as an available resource.
 * The thread waits for that condition, until it is unblocked by another
 * thread that caused the condition to exist (<i>i.e.</i> freed the
 * resource).
 */

#include "Semaphore.h"
#include "Mutex.h"



/**************************************
 
CLASS
   ConditionVariable
   
KEYWORDS
   ConditionVariable
   
DESCRIPTION
   Condition variable primitive.  When a thread calls the
   <i>wait</i> method,which will block until another thread calls
   the <i>cond_signal</i> or <i>cond_broadcast</i> methods.  When
   there are multiple threads waiting, <i>cond_broadcast</i> will unblock
   all of them, while <i>cond_signal</i> will unblock only one (an
   arbitrary one) of them.  This primitive is used to allow a thread
   to wait for some condition to exist, such as an available resource.
   The thread waits for that condition, until it is unblocked by another
   thread that caused the condition to exist (<i>i.e.</i> freed the
   resource).
PATTERNS


WARNING
   
****************************************/

class SCICORESHARE ConditionVariable {
    const char* name;
    int nwaiters;
    Mutex mutex;
    Semaphore semaphore;
public:
    //////////
    //Create a condition variable <i>name</i> should be a static string which describes the
    //primitive for debugging purposes.
    ConditionVariable(const char* name) ;
    
    //////////
    //Destroy the condition variable
    ~ConditionVariable() ;
    
    //////////
    //Wait for a condition.  This method atomically unlicks <b>mutex</i>, and blocks.  The
    //<b>mutex</b> is typically used to guard access to the resource that the thread is
    //waiting for.
    void wait(Mutex& m) ;
    
    //////////
    //Signal a condition.  This will unblock one of the waiting threads.  No guarantee is made
    //as to which thread will be unblocked, but thread implementations should give preference
    //to the thread that has waited the longest.
    void cond_signal()  ;

    //////////
    //Signal a condition.  This will unblock all of the waiting threads.  Note that only the
    //number of waiting threads will be unblocked.  No guarantee is made that these are the
    //N threads that were blocked at the time of the broadcast.
    void cond_broadcast() ;
};

#endif


