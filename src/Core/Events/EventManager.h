//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2004 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : EventManager.h
//    Author : Martin Cole, McKay Davis
//    Date   : Wed May 24 07:58:40 2006

#include <Core/Thread/Runnable.h>
#include <Core/Thread/Mailbox.h>
#include <Core/Events/BaseEvent.h>
#include <Core/Events/Tools/ToolManager.h>
#include <string>
#include <map>

#if !defined(EventManager_h)
#define EventManager_h

namespace SCIRun {

using namespace std;

class EventManager : public Runnable {
public:
  typedef Mailbox<event_handle_t> event_mailbox_t;

  EventManager();
  ~EventManager();

  //! The calling object registers with a unique string id, and recieves 
  //! the mailbox that event messages come through.
  static event_mailbox_t* register_event_messages(string);
  //! Trigger the shared mailbox for the unique string id to be destroyed.
  static void             unregister_event_messages(string);

  static void add_event(event_handle_t e) 
  {
    mailbox_.send(e);
  }
  
  virtual void run();
private:
  typedef map<string, event_mailbox_t*>    id_tm_map_t;

  //! all of the threads who need to know about events.
  static id_tm_map_t                   mboxes_;

  //! the mailbox for adding events to the stream.
  static event_mailbox_t               mailbox_;

  //! for tools that process or modify events before being dispatched.
  ToolManager                          tm_;
};

} // namespace SCIRun

#endif //EventManager_h
