
/*
 *  LockingHandle.h: Smart Pointers..
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#ifndef SCI_Containers_LockingHandle_h
#define SCI_Containers_LockingHandle_h 1

#include <config.h>

#include <SCICore/Util/Assert.h>

namespace SCICore {

namespace PersistentSpace {
  class Piostream;
}

namespace Containers {

using SCICore::PersistentSpace::Piostream;

template<class T>
class LockingHandle {
    T* rep;
public:
    LockingHandle();
    LockingHandle(T*);
    LockingHandle(const LockingHandle<T>&);
    LockingHandle<T>& operator=(const LockingHandle<T>&);
    LockingHandle<T>& operator=(T*);
    ~LockingHandle();

    void detach();

    inline T* operator->() const {
	ASSERT(rep != 0);
	return rep;
    }
    inline T* get_rep() const {
	return rep;
    }

    friend void TEMPLATE_TAG Pio TEMPLATE_BOX (Piostream& stream, LockingHandle<T>& data);
};

} // End namespace Containers
} // End namespace SCICore

////////////////////////////////////////////////////////////
//
// Start of included LockingHandle.cc
//

#include <SCICore/Persistent/Persistent.h>
#include <iostream.h>

namespace SCICore {
namespace Containers {

template<class T>
LockingHandle<T>::LockingHandle()
: rep(0)
{
}

template<class T>
LockingHandle<T>::LockingHandle(T* rep)
: rep(rep)
{
    if(rep){
	rep->lock.lock();
	rep->ref_cnt++;
	rep->lock.unlock();
    }
}

template<class T>
LockingHandle<T>::LockingHandle(const LockingHandle<T>& copy)
: rep(copy.rep)
{
    if(rep){
	rep->lock.lock();
	rep->ref_cnt++;
	rep->lock.unlock();
    }
}

template<class T>
LockingHandle<T>& LockingHandle<T>::operator=(const LockingHandle<T>& copy)
{
    if(rep != copy.rep){
	if(rep){
	    rep->lock.lock();
	    if(--rep->ref_cnt==0){
		rep->lock.unlock();
		delete rep;
	    } else {
		rep->lock.unlock();
	    }
	}
	if(copy.rep){
	    copy.rep->lock.lock();
	    rep=copy.rep;
	    rep->ref_cnt++;
	    copy.rep->lock.unlock();
	} else {
	    rep=copy.rep;
	}
    }
    return *this;
}

template<class T>
LockingHandle<T>& LockingHandle<T>::operator=(T* crep)
{
    if(rep){
	rep->lock.lock();
	if(--rep->ref_cnt==0){
	    rep->lock.unlock();
	    delete rep;
	} else {
	    rep->lock.unlock();
	}
    }
    if(crep){
	crep->lock.lock();
	rep=crep;
	rep->ref_cnt++;
	crep->lock.unlock();
    } else {
	rep=crep;
    }
    return *this;
}

template<class T>
LockingHandle<T>::~LockingHandle()
{
    if(rep){
	rep->lock.lock();
	if(--rep->ref_cnt==0){
	    rep->lock.unlock();
	    delete rep;
	} else {
	    rep->lock.unlock();
	}
    }
}

template<class T>
void LockingHandle<T>::detach()
{
    ASSERT(rep != 0);
    rep->lock.lock();
    if(rep->ref_cnt==1){
	rep->lock.unlock();
	return; // We have the only copy
    }
    T* oldrep=rep;
    rep=oldrep->clone();
    oldrep->ref_cnt--;
    oldrep->lock.unlock();
    rep->ref_cnt++;
}

template<class T>
void Pio(Piostream& stream, LockingHandle<T>& data)
{
    stream.begin_cheap_delim();
    PersistentSpace::Persistent* trep=data.rep;
    stream.io(trep, T::type_id);
    if(stream.reading()){
	data.rep=(T*)trep;
	if(data.rep)
	    data.rep->ref_cnt++;
    }
    stream.end_cheap_delim();
}

} // End namespace Containers
} // End namespace SCICore

//
// $Log$
// Revision 1.2  1999/08/17 06:38:36  sparker
// Merged in modifications from PSECore to make this the new "blessed"
// version of SCIRun/Uintah.
//
// Revision 1.1  1999/07/27 16:56:12  mcq
// Initial commit
//
// Revision 1.4  1999/07/07 21:10:35  dav
// added beginnings of support for g++ compilation
//
// Revision 1.3  1999/05/06 19:55:43  dav
// added back .h files
//
// Revision 1.1  1999/05/05 21:04:31  dav
// added SCICore .h files to /include directories
//
// Revision 1.1.1.1  1999/04/24 23:12:26  dav
// Import sources
//
//

#endif
