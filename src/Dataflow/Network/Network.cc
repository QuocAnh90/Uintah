/*
  The contents of this file are subject to the University of Utah Public
  License (the "License"); you may not use this file except in compliance
  with the License.
  
  Software distributed under the License is distributed on an "AS IS"
  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
  License for the specific language governing rights and limitations under
  the License.
  
  The Original Source Code is SCIRun, released March 12, 2001.
  
  The Original Source Code was developed by the University of Utah.
  Portions created by UNIVERSITY are Copyright (C) 2001, 1994 
  University of Utah. All Rights Reserved.
*/


/*
 *  Network.cc: The core of dataflow...
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 *  Distributed modifications by:
 *   Michelle Miller
 *   Dec. 1997
 *  Copyright (C) 1994 SCI Group
 */

#include <Dataflow/Network/Network.h>
#include <Core/Util/NotFinished.h>
#include <Dataflow/Network/Connection.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/PackageDB.h>
#include <Core/Malloc/Allocator.h>
#include <Core/GuiInterface/Remote.h>

#ifdef _WIN32
#include <io.h>
#endif
#include <stdio.h>
#include <iostream>
using std::cerr;
#include <string.h>

//#define DEBUG 1

namespace SCIRun {


Network::Network(int first)
  : the_lock("Network lock"),
    netedit(0), first(first), nextHandle(0), slave_socket(0)
{
}

Network::~Network()
{
    if (slave_socket != 0)
	close (slave_socket);
}

int Network::read_file(const string&)
{
    NOT_FINISHED("Network::read_file");
    return 1;
}

// For now, we just use a simple mutex for both reading and writing
void Network::read_lock()
{
    the_lock.lock();
}

void Network::read_unlock()
{
    the_lock.unlock();
}

void Network::write_lock()
{
    the_lock.lock();
}

void Network::write_unlock()
{
    the_lock.unlock();
}

int Network::nmodules()
{
    return modules.size();
}

Module* Network::module(int i)
{
    return modules[i];
}

int Network::nconnections()
{
    return connections.size();
}

Connection* Network::connection(int i)
{
    return connections[i];
}

string Network::connect(Module* m1, int p1, Module* m2, int p2)
{

    Connection* conn=scinew Connection(m1, p1, m2, p2);
    string id(m1->id+"_p"+to_string(p1)+"_to_"+m2->id+"_p"+to_string(p2));
    conn->id=id;
    conn->connect();
    connections.add(conn);
    // Reschedule next time we can...
    reschedule=1;

    // socket already set up
    Message msg;

    if (!m1->isSkeleton() && !m2->isSkeleton()) {
    }  	// do nothing - normal case
    else if (m1->isSkeleton() && m2->isSkeleton()) {
	// format CreateLocalConnectionMsg control message
        msg.type        = CREATE_LOC_CONN;
	msg.u.clc.outModHandle = m1->handle;
	msg.u.clc.oport = p1;
	msg.u.clc.inModHandle = m2->handle;
	msg.u.clc.iport = p2;
        strcpy (msg.u.clc.connID, id.c_str());
        msg.u.clc.connHandle = conn->handle = getNextHandle();
        conn_handles[conn->handle] = conn;
	

	// send control message
        char buf[BUFSIZE];
        bzero (buf, sizeof (buf));
        bcopy ((char *) &msg, buf, sizeof (msg));
        write (slave_socket, buf, sizeof(buf));
    } else {
	conn->setRemote();

	// format CreateRemoteConnectionMsg control message
        msg.type        = CREATE_REM_CONN;
	if (m1->isSkeleton()) {
	    msg.u.crc.fromRemote = true; 	// connection from slave
	    m2->handle = getNextHandle();	// assign handle to local mod
	    mod_handles[m2->handle] = m2;       // take this out - don't need
	} else {
	    msg.u.crc.fromRemote = false; 	// connection from master
	    m1->handle = getNextHandle();	// assign handle to local mod
	    mod_handles[m1->handle] = m1;       // take this out - don't need
	}
	    
        msg.u.crc.outModHandle = m1->handle;
        msg.u.crc.oport = p1;
        msg.u.crc.inModHandle = m2->handle;
        msg.u.crc.iport = p2;
        strcpy (msg.u.crc.connID, id.c_str());
        msg.u.crc.connHandle = conn->handle = getNextHandle();
        conn_handles[conn->handle] = conn;
	msg.u.crc.socketPort = conn->socketPort = BASE_PORT + conn->handle;

	// send control message
        char buf[BUFSIZE];
        bzero (buf, sizeof (buf));
        bcopy ((char *) &msg, buf, sizeof (msg));
        write (slave_socket, buf, sizeof(buf));

	// setup data connection
#ifndef _WIN32
	int listen_socket = setupConnect (conn->socketPort);
	conn->remSocket = acceptConnect (listen_socket);
	close (listen_socket);
#endif
    }
    return id;
}

int Network::disconnect(const string& connId)
{
    Message msg;

    int i;
    for (i = 0; i < connections.size(); i++)
        if (connections[i]->id == connId)
		break;
    if (i == connections.size()) {
        return 0;
    }
 
    // endpoints on two different machines - format DeleteRemConnMsg
    if (connections[i]->isRemote()) {
        msg.type        = DELETE_REM_CONN;
        msg.u.drc.connHandle = connections[i]->handle;

    	// send message to slave
	char buf[BUFSIZE];
    	bzero (buf, sizeof (buf));
    	bcopy ((char *) &msg, buf, sizeof (msg));
    	write (slave_socket, buf, sizeof(buf));
	
    	// remove from handle hash table
    	conn_handles.erase(connections[i]->handle);

    // only distrib connections will have a handle - format DeleteLocConnMsg
    } else if (connections[i]->handle) {
	msg.type 	= DELETE_LOC_CONN;
	msg.u.dlc.connHandle = connections[i]->handle;

    	// send message to slave
    	char buf[BUFSIZE];
    	bzero (buf, sizeof (buf));
    	bcopy ((char *) &msg, buf, sizeof (msg));
    	write (slave_socket, buf, sizeof(buf));

    	// remove from handle hash table
    	conn_handles.erase(connections[i]->handle);
    } 
    // remove connection ref from iport and oport

    delete connections[i];	//connection destructor,tears down data channel
    connections.remove(i);
    return 1;
}

void Network::initialize(NetworkEditor* _netedit)
{
    netedit=_netedit;
    //NOT_FINISHED("Network::initialize"); // Should read a file???
}


static string
remove_spaces(const string& str)
{
  string result;
  for (string::const_iterator i = str.begin(); i != str.end(); i++)
  {
    if (*i != ' ') { result.push_back(*i); }
  }
  return result;
}


Module* Network::add_module(const string& packageName,
                            const string& categoryName,
                            const string& moduleName)
{ 

  // Find a unique id in the Network for the new instance of this module and
  // form an instance name from it

  string instanceName;
  {
    const string name = remove_spaces(packageName + "_" +
				      categoryName + "_" +
				      moduleName + "_");
    for (int i=0; get_module_by_id(instanceName = name + to_string(i)); i++);
  }

  // Instantiate the module

  Module* mod = packageDB.instantiateModule(packageName, categoryName,
					    moduleName, instanceName);
  if(!mod) {
    cerr << "Error: can't create instance " << instanceName << "\n";
    return 0;
  }
  modules.add(mod);

//-----------------------------------------------------------------------------
// XXX: McQ, 7/19/99.  We need to re-examine things like "if(name(0)=='r')"
//      with an eye toward robust module naming.
#if 0
    Message msg;

#ifdef DEBUG
    cerr << "Network::add_module created ID\n";
#endif
    if (name(0) == 'r') {

#ifdef DEBUG
    	cerr << "Network::add_module remote module\n";
#endif

#ifndef _WIN32
	// open listen socket and startup slave if not done yet
	if (slave_socket == 0) {
	    int listen_socket = setupConnect (BASE_PORT);

/* Thread "TCLTask"(pid 8291) caught signal SIGSEGV at address 6146200 (segmentation violation - Unknown code!)
	    // rsh to startup sr, passing master port number 
	    system ("rsh burn /a/home/sci/data12/mmiller/o_Dataflow/sr -slave burn 8888");
 */ 	    slave_socket = acceptConnect (listen_socket);
            close (listen_socket);
	}
#endif

	// send message to addModule (pass ID to link 2 instances);
	msg.type 	= CREATE_MOD;
	strcpy (msg.u.cm.name, name());
	strcpy (msg.u.cm.id, id());
	msg.u.cm.handle = getNextHandle();

	char buf[BUFSIZE];
	bzero (buf, sizeof (buf));
	bcopy ((char *) &msg, buf, sizeof (msg));
   	write (slave_socket, buf, sizeof(buf));

	// pass through to normal code - skeleton module? yes, need to do
	// this so the NetworkEditor calls for remote will execute.  if I
	// don't bind a NetworkEditor, how will the calls happen?  the other
	// side (daemon) should repeat code below also.
    }

// XXX: McQ: Module instantiated moved from here to above

    if (name(0) == 'r') {		// is this a remote module?
	mod->skeleton = true;		// tag as skeleton
    	mod->handle = msg.u.cm.handle;	// skeleton & remote same handle
    }
#endif
//-----------------------------------------------------------------------------

    // Binds NetworkEditor and Network instances to module instance.  
    // Instantiates ModuleHelper and starts event loop.
    mod->set_context(netedit, this);   

    // add Module id and ptr to Module to hash table of modules in network
    module_ids[mod->id] = mod;

    // add to hash table of handles and module ptrs
    if (mod->handle > 0) {
      mod_handles[mod->handle] = mod;
    }
    
    return mod;
}

Module* Network::get_module_by_id(const string& id)
{
    MapStringModule::iterator mod;
    mod = module_ids.find(id);
    if (mod != module_ids.end()) {
	return (*mod).second;
    } else {
	return 0;
    }
}

Module* Network::get_module_by_handle (int handle)
{
    MapIntModule::iterator mod;
    mod = mod_handles.find(handle);
    if (mod != mod_handles.end()) {
	return (*mod).second;
    } else {
	return 0;
    }
}

Connection* Network::get_connect_by_handle (int handle)
{
    MapIntConnection::iterator conn;
    conn = conn_handles.find(handle);
    if (conn != conn_handles.end()) {
  	return (*conn).second;
    } else {
	return 0;
    }
}

int Network::delete_module(const string& id)
{
    Module* mod = get_module_by_id(id);
    if (!mod)
	return 0;
    
    if (mod->isSkeleton()) {

	// format deletemodule message
    	Message msg;
        msg.type        = DELETE_MOD;
        msg.u.dm.modHandle = mod->handle;

        // send msg to slave
        char buf[BUFSIZE];
        bzero (buf, sizeof (buf));
        bcopy ((char *) &msg, buf, sizeof (msg));
        write (slave_socket, buf, sizeof(buf));
	
	// remove from handle hash table
  	mod_handles.erase(mod->handle);
    }

    // traverse array of ptrs to Modules in Network to find this module
    int i;
    for (i = 0; i < modules.size(); i++)
        if (modules[i] == mod)
	    break;
    if (i == modules.size())
	return 0;

    // remove array element corresponding to module, remove from hash table
    modules.remove(i);
    module_ids.erase(id);
    delete mod;			
    return 1;
}

} // End namespace SCIRun

