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

#ifndef Core_CCA_Component_Comm_NexusHandlerThread_h
#define Core_CCA_Component_Comm_NexusHandlerThread_h

#include <Core/CCA/Component/Comm/EpChannel.h>
#include <Core/Thread/Runnable.h>

namespace SCIRun {
  class NexusEpChannel;
  class NexusHandlerThread : public Runnable {
    EpChannel::HPF hfunc;
    Message* msg;
  public:
    NexusHandlerThread(NexusEpChannel *chan, int h_id);
    virtual void run();
  };
}

#endif
