
/*
 *  SlaveController.h: Daemon - Central Slave interface to the Master's
 *   scheduler.
 *
 *  Keeps a socket to the scheduler on the Master for control information.
 *  Dispatches control directives (messages) sent from the Master to direct
 *  remote SCIRun networks (running on a Slave).
 *
 *  Written by:
 *   Michelle Miller
 *   Department of Computer Science
 *   University of Utah
 *   Feb. 1998
 *
 *  Copyright (C) 1998 SCI Group
 */

#ifndef SCI_project_SlaveController_h
#define SCI_project_SlaveController_h

#include <PSECore/share/share.h>

#include <PSECore/Comm/MessageBase.h>
#include <SCICore/Multitask/Task.h>
#include <SCICore/Multitask/Mailbox.h>
#include <SCICore/TclInterface/Remote.h>

namespace PSECore {
  namespace Dataflow {
    class Module;
    class OPort; 
    class Connection;
    class Network;
  }
  namespace Comm {
    class MessageBase;
  }
}

namespace PSECore {
namespace Distributed {

using PSECore::Dataflow::Module;
using PSECore::Dataflow::OPort;
using PSECore::Dataflow::Connection;
using PSECore::Dataflow::Network;
using PSECore::Comm::MessageBase;

using SCICore::Multitask::Task;
using SCICore::Multitask::Mailbox;

class SCICORESHARE SlaveController : public Task {
    private:
	Network* 	net;
	char 		masterHost[HOSTNAME];	
	int 		masterPort;
	int 		master_socket;
    public:
    	Mailbox<MessageBase*> mailbox;
	SlaveController (Network*, char*, char*);
	~SlaveController();
    private:
	int createModule (char *name, char* id, int handle);
	int createLocConnect (int, int, int, int, char*, int);
	int createRemConnect (bool, int, int, int, int, char*, int, int);
	int executeModule (int handle);
	int triggerPort (int modHandle, int connHandle);
	int deleteModule (int modHandle);
 	int deleteRemConnect (int connHandle);
	int deleteLocConnect (int connHandle);
	int main_loop();
	virtual int body(int);
};

class SCICORESHARE R_Scheduler_Module_Message : public MessageBase {
public:
    Connection* conn;
    R_Scheduler_Module_Message();			// execute msg
    R_Scheduler_Module_Message(Connection* conn);	// trigger port
    virtual ~R_Scheduler_Module_Message();
};

class SCICORESHARE R_Module_Scheduler_Message : public MessageBase {
public:
    OPort* p1;
    OPort* p2;
    R_Module_Scheduler_Message();			// reschedule
    R_Module_Scheduler_Message(OPort*, OPort*); 	// multisend
    virtual ~R_Module_Scheduler_Message();
};

} // End namespace Distributed
} // End namespace PSECore

//
// $Log$
// Revision 1.2  1999/08/17 06:38:26  sparker
// Merged in modifications from PSECore to make this the new "blessed"
// version of SCIRun/Uintah.
//
// Revision 1.1  1999/07/27 16:56:01  mcq
// Initial commit
//
// Revision 1.3  1999/05/06 22:05:50  dav
// added back .h file
//
// Revision 1.1.1.1  1999/04/24 23:12:52  dav
// Import sources
//
//

#endif
