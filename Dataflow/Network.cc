
/*
 *  Network.cc: The core of dataflow...
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <Dataflow/Network.h>

#include <Classlib/NotFinished.h>
#include <Dataflow/Connection.h>
#include <Dataflow/Module.h>
#include <Dataflow/ModuleList.h>
#include <Malloc/Allocator.h>

#include <iostream.h>
#include <stdlib.h>

Network::Network(int first)
: netedit(0), first(first)
{
}

Network::~Network()
{
}

int Network::read_file(const clString&)
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

clString Network::connect(Module* m1, int p1, Module* m2, int p2)
{
    Connection* conn=scinew Connection(m1, p1, m2, p2);
    clString id(m1->id+"_p"+to_string(p1)+"_to_"+m2->id+"_p"+to_string(p2));
    conn->id=id;
    conn->connect();
    connections.add(conn);
    // Reschedule next time we can...
    reschedule=1;
    return id;
}

void Network::initialize(NetworkEditor* _netedit)
{
    netedit=_netedit;
    NOT_FINISHED("Network::initialize"); // Should read a file???
}

Module* Network::add_module(const clString& name)
{
    makeModule maker=ModuleList::lookup(name);
    if(!maker){
	cerr << "Module: " << name << " not found!\n";
	return 0;
    }

    // Create the ID...
    int i=0;
    clString id;
    while(1){
	id=clString(name+"_"+to_string(i));
	if(get_module_by_id(id)){
	    i++;
	} else {
	    break;
	}
    }
    Module* mod=(*maker)(id);
    modules.add(mod);
    mod->set_context(netedit, this);
    module_ids.insert(mod->id, mod);
    return mod;
}

Module* Network::get_module_by_id(const clString& id)
{
    Module* mod;
    if(module_ids.lookup(id, mod)){
	return mod;
    } else {
	return 0;
    }
}

int Network::delete_module(const clString&)
{	
    NOT_FINISHED("Network::delete_module");
    return 0;
}
