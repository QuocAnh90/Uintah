
/*
 *  Task_pthreads.cc:  Task implementation for pthreads
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1996
 *
 *  Copyright (C) 1996 SCI Group
 */

#ifdef __linux
#define __KERNEL__
#include <linux/signal.h>
typedef __sighandler_t SIG_PF;
typedef struct sigcontext_struct sigcontext_t;
#endif

#include <Multitask/Task.h>
#include <Multitask/ITC.h>
#include <Malloc/Allocator.h>
#include <iostream.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <pthread.h>
#include <errno.h>
#include <string.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#ifdef __sgi
#include <sys/sysmp.h>
#endif

#ifdef __linux
extern "C" {
    size_t getpagesize();
};
#endif

#define DEFAULT_STACK_LENGTH 64*1024
#define INITIAL_STACK_LENGTH 16*1024
#define DEFAULT_SIGNAL_STACK_LENGTH 16*1024

#define NOT_FINISHED(what) cerr << what << ": Not Finished " << __FILE__ << " (line " << __LINE__ << ") " << endl

#if 0

extern "C" int Task_try_lock(unsigned long*);

static int aborting=0;
static int exit_code;
static int exiting=0;
static long tick;
#endif
static int pagesize;
static int devzero_fd;
static char* progname;
static pthread_key_t selfkey;

struct TaskPrivate {
    pthread_t threadid;
    caddr_t sp;
    caddr_t stackbot;
    size_t stacklen;
    size_t redlen;
};

struct TaskArgs {
    int arg;
    Task* t;
};

#define MAXTASKS 1000
static int ntasks=0;
static Task* tasks[MAXTASKS];
Mutex* sched_lock;

static void makestack(TaskPrivate* priv)
{
    priv->stacklen=DEFAULT_STACK_LENGTH+pagesize;
    priv->stackbot=(caddr_t)mmap(0, priv->stacklen, PROT_READ|PROT_WRITE,
				 MAP_SHARED, devzero_fd, 0);
    if((long)priv->stackbot == -1){
	perror("mmap");
	Task::exit_all(-1);
    }
    // Now unmap the bottom part of it...
    priv->redlen=DEFAULT_STACK_LENGTH-INITIAL_STACK_LENGTH;
    if(getenv("SCI_CATCH_SIGNALS")){
	priv->redlen=DEFAULT_STACK_LENGTH-INITIAL_STACK_LENGTH;
    } else {
	priv->redlen=pagesize;
    }
    priv->sp=(char*)priv->stackbot+priv->stacklen-1;
    if(mprotect(priv->stackbot, priv->redlen, PROT_NONE) == -1){
	perror("mprotect");
	Task::exit_all(-1);
    }
}

void* runbody(void* vargs)
{
    TaskArgs* args=(TaskArgs*)vargs;
    int arg=args->arg;
    Task* t=args->t;
    delete args;

    t->startup(arg);
    return 0;
}

int Task::startup(int task_arg)
{
    if(pthread_setspecific(selfkey, (void*)this) != 0){
	perror("pthread_setspecific");
	exit(1);
    }
    int retval=body(task_arg);
    Task::taskexit(this, retval);
    return 0; // Never reached.
}

// We are done..
void Task::taskexit(Task*, int retval)
{
    pthread_exit((void*)retval);
} 

// Fork a thread
void Task::activate(int task_arg)
{
    if(activated){
	cerr << "Error: task is being activated twice!" << endl;
	Task::exit_all(-1);
    }
    activated=1;
    priv=scinew TaskPrivate;
    TaskArgs* args=scinew TaskArgs;
    args->arg=task_arg;
    args->t=this;

    makestack(priv);

    sched_lock->lock();
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr, priv->stacklen);
    pthread_attr_setstackaddr(&attr, priv->sp);
    if(pthread_create(&priv->threadid, NULL, runbody, (void*)args) != 0) {
	perror("pthread_create");
	exit(1);
    }

    tasks[ntasks++]=this;
    sched_lock->unlock();
}

Task* Task::self()
{
    Task* t=(Task*)pthread_getspecific(selfkey);
    if(!t){
	perror("pthread_getspecific");
	exit(1);
    }
    return t;
}

static char* signal_name(int sig, int code, caddr_t addr)
{
    static char buf[1000];
    switch(sig){
    case SIGHUP:
	sprintf(buf, "SIGHUP (hangup)");
	break;
    case SIGINT:
	sprintf(buf, "SIGINT (interrupt)");
	break;
    case SIGQUIT:
	sprintf(buf, "SIGQUIT (quit)");
	break;
    case SIGILL:
	sprintf(buf, "SIGILL (illegal instruction)");
	break;
    case SIGTRAP:
	sprintf(buf, "SIGTRAP (trace trap)");
	break;
    case SIGABRT:
	sprintf(buf, "SIBABRT (Abort)");
	break;
    case SIGFPE:
	sprintf(buf, "SIGFPE (Floating Point Exception)");
	break;
    case SIGKILL:
	sprintf(buf, "SIGKILL (kill)");
	break;
    case SIGBUS:
	sprintf(buf, "SIGBUS (bus error)");
	break;
    case SIGSEGV:
	{
	    char* why;
	    switch(code){
	    case EFAULT:
		why="Invalid virtual address";
		break;
	    case EACCES:
		why="Invalid permissions for address";
		break;
	    default:
		why="Unknown code!";
		break;
	    }
	    sprintf(buf, "SIGSEGV at address %p (segmentation violation - %s)",
		    addr, why);
	}
	break;
    case SIGPIPE:
	sprintf(buf, "SIGPIPE (broken pipe)");
	break;
    case SIGALRM:
	sprintf(buf, "SIGALRM (alarm clock)");
	break;
    case SIGTERM:
	sprintf(buf, "SIGTERM (killed)");
	break;
    case SIGUSR1:
	sprintf(buf, "SIGUSR1 (user defined signal 1)");
	break;
    case SIGUSR2:
	sprintf(buf, "SIGUSR2 (user defined signal 2)");
	break;
    case SIGCLD:
	sprintf(buf, "SIGCLD (death of a child)");
	break;
    case SIGPWR:
	sprintf(buf, "SIGPWR (power fail restart)");
	break;
    case SIGWINCH:
	sprintf(buf, "SIGWINCH (window size changes)");
	break;
    case SIGURG:
	sprintf(buf, "SIGURG (urgent condition on IO channel)");
	break;
    case SIGIO:
	sprintf(buf, "SIGIO (IO possible)");
	break;
    case SIGSTOP:
	sprintf(buf, "SIGSTOP (sendable stop signal)");
	break;
    case SIGTSTP:
	sprintf(buf, "SIGTSTP (TTY stop)");
	break;
    case SIGCONT:
	sprintf(buf, "SIGCONT (continue)");
	break;
    case SIGTTIN:
	sprintf(buf, "SIGTTIN");
	break;
    case SIGTTOU:
	sprintf(buf, "SIGTTOU");
	break;
    case SIGVTALRM:
	sprintf(buf, "SIGVTALRM (virtual time alarm)");
	break;
    case SIGPROF:
	sprintf(buf, "SIGPROF (profiling alarm)");
	break;
    case SIGXCPU:
	sprintf(buf, "SIGXCPU (CPU time limit exceeded)");
	break;
    case SIGXFSZ:
	sprintf(buf, "SIGXFSZ (Filesize limit exceeded)");
	break;
    default:
	sprintf(buf, "unknown signal(%d)", sig);
	break;
    }
    return buf;
}

static void handle_halt_signals(int sig, int code, sigcontext_t* context)
{
    Task* self=Task::self();
    char* tname=self?self->get_name():"main";
    // Kill all of the threads...
#ifdef __sgi
#if defined(_LONGLONG)
    caddr_t addr=(caddr_t)context->sc_badvaddr;
#else
    caddr_t addr=(caddr_t)context->sc_badvaddr.lo32;
#endif
#endif
#ifdef __linux
    caddr_t addr=(caddr_t)context->eip;
#endif
    char* signam=signal_name(sig, code, addr);
    fprintf(stderr, "Thread \"%s\"(pid %d) caught signal %s.  Going down...\n", tname, getpid(), signam);
    Task::exit_all(-1);
}

static void handle_abort_signals(int sig, int code, sigcontext_t* context)
{
    Task* self=Task::self();
    char* tname=self?self->get_name():"main";
#ifdef __sgi
#if defined(_LONGLONG)
    caddr_t addr=(caddr_t)context->sc_badvaddr;
#else
    caddr_t addr=(caddr_t)context->sc_badvaddr.lo32;
#endif
#endif
#ifdef __linux
    caddr_t addr=(caddr_t)context->eip;
#endif
    char* signam=signal_name(sig, code, addr);

    // See if it is a segv on the stack - if so, grow it...
    if(sig==SIGSEGV && code==EACCES
       && addr >= self->priv->stackbot
       && addr < self->priv->stackbot+self->priv->stacklen){
	self->priv->redlen -= pagesize;
	if(self->priv->redlen <= 0){
	    fprintf(stderr, "%c%c%cThread \"%s\"(pid %d) ran off end of stack! \n",
		    7,7,7,tname, getpid());
	    fprintf(stderr, "Stack size was %d bytes\n", self->priv->stacklen-pagesize);
	} else {
	    if(mprotect(self->priv->stackbot+self->priv->redlen, pagesize,
			PROT_READ|PROT_WRITE) == -1){
		fprintf(stderr, "Error extending stack for thread \"%s\"", tname);
		Task::exit_all(-1);
	    }
	    fprintf(stderr, "extended stack for thread %s\n", tname);
	    fprintf(stderr, "stacksize is now %d bytes\n",
		    self->priv->stacklen-self->priv->redlen);
	    return;
	}
    }

    // Ask if we should abort...
    fprintf(stderr, "%c%c%cThread \"%s\"(pid %d) caught signal %s\ndump core? ", 7,7,7,tname, getpid(), signam);
    char buf[100];
    buf[0]='n';
    if(!fgets(buf, 100, stdin)){
	// Exit without an abort...
	Task::exit_all(-1);
    }
    if(buf[0] == 'n' || buf[0] == 'N'){
	// Exit without an abort...
	Task::exit_all(-1);
    } else if(buf[0] == 'd' || buf[0] == 'D') {
	// Start the debugger...
	Task::debug(Task::self());
	while(1){
	    sigsuspend(0);
	}
    } else {
	// Abort...
	fprintf(stderr, "Dumping core...\n");
	struct rlimit rlim;
	getrlimit(RLIMIT_CORE, &rlim);
	rlim.rlim_cur=RLIM_INFINITY;
	setrlimit(RLIMIT_CORE, &rlim);
//	signal(SIGABRT, SIG_DFL); // We will dump core, but not other threads
	kill(0, SIGABRT);
	sigsuspend(0);
    }
}

void Task::exit_all(int)
{
    fprintf(stderr, "Task::exit_all not done!\n");
    exit(1);
}

void Task::initialize(char* pn)
{
    progname=strdup(pn);
    pagesize=getpagesize();
    fprintf(stderr, "Open...\n");
    devzero_fd=open("/dev/zero", O_RDWR);
    if(devzero_fd == -1){
	perror("open");
	exit(-1);
    }
    fprintf(stderr, "Done...\n");
    if(pthread_key_create(&selfkey, NULL) != 0){
	perror("pthread_key_create\n");
	exit(1);
    }
    sched_lock=new Mutex;

#ifdef __sgi
    // Set up the signal stack so that we will be able to 
    // Catch the SEGV's that need to grow the stacks...
    int stacklen=DEFAULT_SIGNAL_STACK_LENGTH;
    caddr_t stackbot=(caddr_t)mmap(0, stacklen, PROT_READ|PROT_WRITE,
				   MAP_SHARED, devzero_fd, 0);
    if((long)stackbot == -1){
	perror("mmap");
	exit(-1);
    }
    stack_t ss;
    ss.ss_sp=stackbot;
    ss.ss_sp+=stacklen-1;
    ss.ss_size=stacklen;
    ss.ss_flags=0;
    if(sigaltstack(&ss, NULL) == -1){
	perror("sigstack");
	exit(-1);
    }
#endif
    
    // Setup the seg fault handler...
    // For SIGQUIT
    // halt all threads
    // signal(SIGINT, (SIG_PF)handle_halt_signals);
    struct sigaction action;
#ifdef __sgi
    action.sa_flags=SA_ONSTACK;
#else
    action.sa_flags=0;
#endif
    sigemptyset(&action.sa_mask);
    
    action.sa_handler=(SIG_PF)handle_halt_signals;
    if(sigaction(SIGQUIT, &action, NULL) == -1){
	perror("sigaction");
	exit(-1);
    }

    // For SIGILL, SIGABRT, SIGBUS, SIGSEGV, 
    // prompt the user for a core dump...
    if(getenv("SCI_CATCH_SIGNALS")){
	action.sa_handler=(SIG_PF)handle_abort_signals;
	if(sigaction(SIGILL, &action, NULL) == -1){
	    perror("sigaction");
	    exit(-1);
	}
	if(sigaction(SIGABRT, &action, NULL) == -1){
	    perror("sigaction");
	    exit(-1);
	}
	if(sigaction(SIGBUS, &action, NULL) == -1){
	    perror("sigaction");
	    exit(-1);
	}
	if(sigaction(SIGSEGV, &action, NULL) == -1){
	    perror("sigaction");
	    exit(-1);
	}
    }
}

int Task::nprocessors()
{
#ifdef __sgi
    static int nproc=-1;
    if(nproc==-1){
	nproc = sysmp(MP_NAPROCS);
    }
    return nproc;
#else
    return 1;
#endif
}

void Task::main_exit()
{
    fprintf(stderr, "Task::main_exit not done!\n");
    exit(1);
}


#if 0
void Task::yield()
{
    sginap(0);
}


int Task::wait_for_task(Task* task)
{
    if(uspsema(task->priv->running) == -1){
	perror("usvsema");
	Task::exit_all(-1);
    }
    if(usvsema(task->priv->running) == -1){
	perror("usvsema");
	Task::exit_all(-1);
    }
    return task->priv->retval;
}
#endif

//
// Semaphore implementation
//
struct Semaphore_private {
    int count;
    int nwaiters;
    pthread_mutex_t lock;
    pthread_cond_t cond;
};

Semaphore::Semaphore(int count)
{
    priv=scinew Semaphore_private;
    pthread_mutex_init(&priv->lock, NULL);
    pthread_cond_init(&priv->cond, NULL);
    priv->count=count;
    priv->nwaiters=0;
}

Semaphore::~Semaphore()
{
    pthread_mutex_destroy(&priv->lock);
    pthread_cond_destroy(&priv->cond);
}

void Semaphore::down()
{
    if(pthread_mutex_lock(&priv->lock) != 0){
	perror("pthread_mutex_lock");
	exit(1);
    }
    while(priv->count==0){
	priv->nwaiters++;
	if(pthread_cond_wait(&priv->cond, &priv->lock) != 0){
	    perror("pthread_cond_wait");
	    exit(1);
	}
	priv->nwaiters--;
    }
    priv->count--;
    if(pthread_mutex_unlock(&priv->lock) != 0){
	perror("pthread_mutex_unlock");
	exit(1);
    }
}

int Semaphore::try_down()
{
    if(pthread_mutex_lock(&priv->lock) != 0){
	perror("pthread_mutex_lock");
	exit(1);
    }
    int result;
    if(priv->count==0){
	result=0;
    } else {
	priv->count--;
	result=1;
    }
    if(pthread_mutex_unlock(&priv->lock) != 0){
	perror("pthread_mutex_unlock");
	exit(1);
    }
    return result;
}

void Semaphore::up()
{
    if(pthread_mutex_lock(&priv->lock) != 0){
	perror("pthread_mutex_lock");
	exit(1);
    }
    priv->count++;
    if(priv->nwaiters){
	if(pthread_cond_signal(&priv->cond) != 0){
	    perror("pthread_cond_signal");
	    exit(1);
	}
    }
    if(pthread_mutex_unlock(&priv->lock) != 0){
	perror("pthread_mutex_unlock");
	exit(1);
    }
}

struct Mutex_private {
    pthread_mutex_t lock;
};

Mutex::Mutex()
{
    priv=new Mutex_private;
    pthread_mutex_init(&priv->lock, NULL);
}

Mutex::~Mutex()
{
    pthread_mutex_destroy(&priv->lock);
    delete priv;
}

void Mutex::lock()
{
    if(pthread_mutex_lock(&priv->lock) != 0)
	perror("pthread_mutex_lock");
}

void Mutex::unlock()
{
    if(pthread_mutex_unlock(&priv->lock) != 0)
	perror("pthread_mutex_unlock");
}

int Mutex::try_lock()
{
    int err=pthread_mutex_trylock(&priv->lock);
    if(err == 0)
	return 1;
    else if(err == EBUSY)
	return 0;
    else
	perror("pthread_mutex_trylock");
    return 0;
}

//
// Condition variable implementation
//
struct ConditionVariable_private {
    pthread_cond_t cond;
};

ConditionVariable::ConditionVariable()
{
    priv=scinew ConditionVariable_private;
    pthread_cond_init(&priv->cond, NULL);
}

ConditionVariable::~ConditionVariable()
{
    if(priv){
	pthread_cond_destroy(&priv->cond);
	delete priv;
    }
}

void ConditionVariable::wait(Mutex& mutex)
{
    if(pthread_cond_wait(&priv->cond, &mutex.priv->lock) != 0)
	perror("pthread_cond_wait");
}

void ConditionVariable::cond_signal()
{
    if(pthread_cond_signal(&priv->cond) != 0)
	perror("pthread_cond_signal");
}

void ConditionVariable::broadcast()
{
    if(pthread_cond_broadcast(&priv->cond) != 0)
	perror("pthread_cond_broadcast");
}

#if 0
void Task::sleep(const TaskTime& time)
{
    sginap((time.secs*1000+time.usecs)/tick);
}

TaskInfo* Task::get_taskinfo()
{
    if(uspsema(sched_lock) == -1){
	perror("uspsema");
	Task::exit_all(-1);
    }
    TaskInfo* ti=scinew TaskInfo(ntasks);
    for(int i=0;i<ntasks;i++){
	ti->tinfo[i].name=tasks[i]->name;
	ti->tinfo[i].stacksize=tasks[i]->priv->stacklen-pagesize;
	ti->tinfo[i].stackused=tasks[i]->priv->stacklen-tasks[i]->priv->redlen;
	ti->tinfo[i].pid=tasks[i]->priv->tid;
	ti->tinfo[i].taskid=tasks[i];
    }
    if(usvsema(sched_lock) == -1){
	perror("usvsema");
	Task::exit_all(-1);
    }
    return ti;
}

void Task::coredump(Task* task)
{
    kill(task->priv->tid, SIGABRT);
}
#endif

void Task::debug(Task*)
{
    char buf[1000];
    char* dbx=getenv("SCI_DEBUGGER");
#ifdef __sgi
    if(!dbx)
	dbx="winterm -c dbx -p %d &";
    sprintf(buf, dbx, getpid(), progname);
#endif
#ifdef __linux
    if(!dbx)
	dbx="xterm -e %s %d &";
    sprintf(buf, dbx, progname, getpid());
#endif
    if(system(buf) == -1)
	perror("system");
}

#if 0

// Interface to select...
int Task::mtselect(int nfds, fd_set* readfds, fd_set* writefds,
		   fd_set* exceptfds, struct timeval* timeout)
{
    // Irix has kernel threads, so select works ok...
    return select(nfds, readfds, writefds, exceptfds, timeout);
}

static void handle_alrm(int, int, sigcontext_t*)
{
    Task* task=Task::self();
    if(task->priv->handle_alrm)
       (*task->priv->handle_alrm)(task->priv->alrm_data);
}

int Task::start_itimer(const TaskTime& start, const TaskTime& interval,
		       void (*handler)(void*), void* cbdata)
{
    if(ntimers > 0){
	NOT_FINISHED("Multiple timers in a single thread");
	return 0;
    }
    ITimer** new_timers=scinew ITimer*[ntimers+1];
    if(timers){
	for(int i=0;i<ntimers;i++){
	    new_timers[i]=timers[i];
	}
	delete[] timers;
    }
    timers=new_timers;
    ITimer* t=scinew ITimer;
    timers[ntimers]=t;
    ntimers++;
    t->start=start;
    t->interval=interval;
    t->handler=handler;
    t->cbdata=cbdata;
    t->id=timer_id++;
    struct itimerval it;
    it.it_interval.tv_sec=interval.secs;
    it.it_interval.tv_usec=interval.usecs;
    it.it_value.tv_sec=start.secs;
    it.it_value.tv_usec=start.usecs;
    struct sigaction action;
    action.sa_flags=SA_ONSTACK;
    sigemptyset(&action.sa_mask);

    action.sa_handler=(SIG_PF)handle_alrm;
    if(sigaction(SIGALRM, &action, NULL)== -1){
	perror("sigaction");
	return 0;
    }
    if(setitimer(ITIMER_REAL, &it, 0) == -1){
	perror("setitimer");
	return 0;
    }
    priv->handle_alrm=handler;
    priv->alrm_data=cbdata;
    return t->id;
}

void Task::cancel_itimer(int which_timer)
{
    priv->handle_alrm=0;
    for(int idx=0;idx<ntimers;idx++){
	if(timers[idx]->id == which_timer)
	    break;
    }
    if(idx==ntimers){
	cerr << "Cancelling bad timer!\n" << endl;
	return;
    }
    struct itimerval it;
    timerclear(&it.it_interval);
    timerclear(&it.it_value);
    if(setitimer(ITIMER_REAL, &it, 0) == -1){
	perror("setitimer");
	return;
    }
    for(int i=idx;i<ntimers-1;i++)
	timers[i]=timers[i+1];
    ntimers--;
}

#endif

int Task::start_itimer(const TaskTime& start, const TaskTime& interval,
		       void (*handler)(void*), void* cbdata)
{
    fprintf(stderr, "Task::start_itimer not done!\n");
    exit(1);
    return 0;
}
