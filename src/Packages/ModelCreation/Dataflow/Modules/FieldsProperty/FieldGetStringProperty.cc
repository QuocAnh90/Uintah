/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2004 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
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

/*
 *  FieldGetFieldPropertyProperty.cc:
 *
 *  Written by:
 *  Jeroen Stinstra
 *
 */

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/StringPort.h>
#include <Core/Datatypes/String.h>
#include <Dataflow/Network/Module.h>
#include <Core/Malloc/Allocator.h>

namespace ModelCreation {

using namespace SCIRun;
using namespace std;

class FieldGetStringProperty : public Module {
public:
  FieldGetStringProperty(GuiContext*);

  virtual ~FieldGetStringProperty();

  virtual void execute();
  
private:
  GuiString             guistring1name_;
  GuiString             guistring2name_;
  GuiString             guistring3name_;
  GuiString             guistrings_;
};


DECLARE_MAKER(FieldGetStringProperty)
  FieldGetStringProperty::FieldGetStringProperty(GuiContext* ctx)
    : Module("FieldGetStringProperty", ctx, Source, "FieldsProperty", "ModelCreation"),
      guistring1name_(get_ctx()->subVar("string1-name")),
      guistring2name_(get_ctx()->subVar("string2-name")),
      guistring3name_(get_ctx()->subVar("string3-name")),
      guistrings_(get_ctx()->subVar("string-selection"))
{
}

FieldGetStringProperty::~FieldGetStringProperty(){
}


void
FieldGetStringProperty::execute()
{
  std::string string1name = guistring1name_.get();
  std::string string2name = guistring2name_.get();
  std::string string3name = guistring3name_.get();
  std::string stringlist;
        
  FieldHandle handle;
  FieldIPort  *iport;
  StringHandle fhandle;
  std::string fstring;
          
  if(!(iport = static_cast<FieldIPort *>(get_input_port("Field"))))
  {
    error("Could not find 'Field' input port");
    return;
  }

  if (!(iport->get(handle)))
  {   
    warning("No field was found on input port");
    return;
  }

  if (handle.get_rep() == 0)
  {   
    warning("Input field is empty");
    return;
  }

  size_t nprop = handle->nproperties();

  for (size_t p=0;p<nprop;p++)
  {
    if(handle->get_property(handle->get_property_name(p),fhandle))
    {
      if (fhandle.get_rep()) stringlist += "{" + handle->get_property_name(p) + "} ";
    }

    if(handle->get_property(handle->get_property_name(p),fstring))
    {
      stringlist += "{" + handle->get_property_name(p) + "} ";
    }
  }

  guistrings_.set(stringlist);
  get_ctx()->reset();
  
  if (handle->is_property(string1name))
  {
    if (handle->get_property(string1name, fstring))
    {
      fhandle = scinew String(fstring);
    }
    else
    {
      handle->get_property(string1name, fhandle);
    }
    send_output_handle("String1", fhandle);
  }
 
  if (handle->is_property(string2name))
  {
    if (handle->get_property(string2name, fstring))
    {
      fhandle = scinew String(fstring);
    }
    else
    {
      handle->get_property(string2name,fhandle);
    }
    send_output_handle("String2", fhandle);
  }
 
  if (handle->is_property(string3name))
  {
    if (handle->get_property(string3name,fstring))
    {
      fhandle = scinew String(fstring);
    }
    else
    {
      handle->get_property(string3name, fhandle);
    }
    send_output_handle("String3", fhandle);
  }
}


} // end namespace
