
#include <Packages/Uintah/CCA/Components/Schedulers/MessageLog.h>
#include <Packages/Uintah/Core/Grid/Patch.h>
#include <Packages/Uintah/CCA/Ports/Output.h>
#include <Packages/Uintah/Core/Parallel/ProcessorGroup.h>
#include <Core/Thread/Time.h>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace Uintah;
using namespace SCIRun;
using namespace std;

MessageLog::MessageLog(const ProcessorGroup* myworld, Output* oport)
   : d_enabled(false), d_myworld(myworld), d_oport(oport)
{
}

MessageLog::~MessageLog()
{
}

void MessageLog::problemSetup(const ProblemSpecP& prob_spec)
{
   ProblemSpecP mess = prob_spec->findBlock("MessageLog");
   if(mess){
      d_enabled = true;
      ostringstream outname;
      outname << d_oport->getOutputLocation() << "/message_" << setw(5) << setfill('0') << d_myworld->myrank() << ".log";
      out.open(outname.str().c_str());
      if(!out){
	 cerr << "Message log disabled, cannot open file\n";
	 d_enabled=false;
      }
      out << "\t\t\t\t  patch\t\t processor\n";
      out << "event\time\t\tsize\tfrom\tto\tfrom\tto\tvariable\totherinfo\n";
   } else {
      d_enabled=false;
   }
}

void MessageLog::logSend(const Task::Dependency* dep, int bytes,
			 const char* msg)
{
   if(!d_enabled)
      return;
   out << "send\t";
   out << setprecision(8) << Time::currentSeconds() << "\t";
   out << bytes << "\t";
   if(dep){
      if(dep->d_task->getPatch())
	 out << dep->d_task->getPatch()->getID();
      else
	 out << "-";
      out << "\t";
      if(dep->d_patch)
	 out << dep->d_patch->getID();
      else
	 out << "-";
      out << "\t" << d_myworld->myrank() << "\t" << dep->d_task->getAssignedResourceIndex() << "\t";
      if(dep->d_var)
	 out << dep->d_var->getName();
      else
	 out << "-";
   } else {
      out << "-\t-\t-\t-\t-";
   }
   if(msg)
      out << "\t" << msg;
   out << '\n';
}

void MessageLog::logRecv(const Task::Dependency* /*dep*/, int /*bytes*/,
			 const char* /*msg*/)
{
   if(!d_enabled)
      return;
}

void MessageLog::finishTimestep()
{
   if(!d_enabled)
      return;
   out.flush();
}

