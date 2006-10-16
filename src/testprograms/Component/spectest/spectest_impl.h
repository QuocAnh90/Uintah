/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2004 Scientific Computing and Imaging Institute,
   University of Utah.

   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/

//    File   : spectest.h
//    Author : Martin Cole
//    Date   : Tue Aug 21 11:17:14 2001

#include <Core/CCA/PIDL/PIDL.h>
#include <testprograms/Component/spectest/spectest_sidl.h>
#include <map>
#include <queue>

using namespace std;

//using spectest::ref;
//using spectest::Server;

using SSIDL::array1;
//using std::istringstream;
//using std::ostringstream;

using namespace SCIRun;
using namespace CCA;
using namespace CCAPORTS;

namespace spectest {

class Port_impl : public CCA::Port {
public:
  Port_impl();
  virtual ~Port_impl();
};

class PortInfo_impl : public CCA::PortInfo {
public:
  PortInfo_impl();
  virtual ~PortInfo_impl();

  virtual ::SSIDL::string getType();
  virtual ::SSIDL::string getName();
  virtual ::SSIDL::string getProperty(const ::SSIDL::string& name);


  // Non standard interface
  void setType(const string &type) { type_ = type; }
  void setName(const string &name) { name_ = name; }
  void setProperties(const ::SSIDL::array1< ::SSIDL::string>& p) 
  { properties_ = p; }
private:
  string              type_;
  string              name_;
  ::SSIDL::array1< ::SSIDL::string> properties_;
};

class Services_impl : public CCA::Services {
public:
  struct PortData {
    PortData() {}
    PortData(const Port::pointer &p, const PortInfo::pointer &pi) :
      port_info_(pi),
      port_(p)
    {}
    PortInfo::pointer port_info_;
    Port::pointer     port_;
  };

  typedef map<string, PortData> port_map_t;

  Services_impl();
  virtual ~Services_impl();
  virtual Port::pointer getPort(const ::SSIDL::string& name);
  virtual Port::pointer getPortNonblocking(const ::SSIDL::string& name);
  virtual PortInfo::pointer createPortInfo(const ::SSIDL::string& name, 
				  const ::SSIDL::string& type, 
				  const ::SSIDL::array1< ::SSIDL::string> &properties );
  virtual void registerUsesPort(const PortInfo::pointer &name_and_type);
  virtual void unregisterUsesPort(const ::SSIDL::string&name);
  virtual void addProvidesPort(const Port::pointer &inPort, const PortInfo::pointer &name);
  virtual void removeProvidesPort(const ::SSIDL::string&name);
  virtual void releasePort(const ::SSIDL::string&name);
  virtual ComponentID::pointer  getComponentID();
private:
  port_map_t provides_;
  port_map_t uses_;
};

class Component_impl : public CCA::Component {
public:
  Component_impl();
  virtual ~Component_impl();
  virtual void setServices(const Services::pointer &svc);
protected:
  Services::pointer services_;
};

class ComponentID_impl : public CCA::ComponentID {
public:
  ComponentID_impl();
  virtual ~ComponentID_impl();
  virtual ::SSIDL::string toString();
};

class GoPort_impl : public CCAPORTS::GoPort {
public:
  GoPort_impl();
  virtual ~GoPort_impl();
};

class ConnectionEventService_impl : 
    public CCAPORTS::ConnectionEventService {
public:
  ConnectionEventService_impl();
  virtual ~ConnectionEventService_impl();

  virtual void addConnectionEventListener(int connectionEventType, 
					  const ConnectionEventListener::pointer &l);
  virtual void removeConnectionEventListener(int connectionEventType, 
					     const ConnectionEventListener::pointer &l);
  };

class ConnectionEventListener_impl : 
    public CCAPORTS::ConnectionEventListener {
public:
  ConnectionEventListener_impl();
  virtual ~ConnectionEventListener_impl();

  virtual void connectionActivity(const ConnectionEvent::pointer &evt);
};

class ConnectionEvent_impl : public CCAPORTS::ConnectionEvent {
public:
  ConnectionEvent_impl();
  virtual ~ConnectionEvent_impl();

  virtual int getEventType();
  virtual CCA::PortInfo::pointer getPortInfo();
             
  // The following enum belongs in the interface, but the sidl compiler
  // cannot currently handle it.
  enum Type { 
    Error = -1,            /* Someone got their hands on a bogus event 
			      object somehow. */
    ALL = 0,               /* Component wants to receive all event notices. 
			      ALL itself never received. */
    ConnectPending = 1,    // A connection is about to be attempted.
    Connected = 2,         // A connection has been made.
    DisconnectPending = 3, // A disconnection is about to be attempted.
    Disconnected = 4       // A disconnection has been made.
  };
  
};

/* ---------------------------------------------------------------------
 * Below is all of the non specification test code.
 * --------------------------------------------------------------------*/

class IntegerStream_impl : public IntegerStream {
  queue<int> stream_;
public:
  IntegerStream_impl();
  virtual ~IntegerStream_impl();
  virtual int pop();
  virtual void push(int i);
  virtual bool is_full();
  virtual bool is_empty();
};

class RandomInt : public Component_impl
{
public:
  RandomInt();
  ~RandomInt();
  void go(); //keep the stream full
  virtual void setServices(const CCA::Services::pointer &svc);
private:
  IntegerStream::pointer istr_;
  // IntegerStream_impl istr_;
};

class ConsumerInt : public Component_impl
{
public:
  ConsumerInt();
  ~ConsumerInt();
  void go(); //kpull from the stream
  virtual void setServices(const CCA::Services::pointer &svc);
private:
  IntegerStream::pointer istr_;
  // IntegerStream_impl istr_;
};

class Framework_impl : public Framework {
  CCA::Services::pointer services_;
public:
  Framework_impl();
  virtual ~Framework_impl();
  
  virtual CCA::Services::pointer get_services();
};

} // end spectest namespace
