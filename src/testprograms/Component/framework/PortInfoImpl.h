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


#ifndef PortInfoImpl_h
#define PortInfoImpl_h

#include <testprograms/Component/framework/cca_sidl.h>

#include <Core/CCA/PIDL/PIDL.h>
#include <Core/CCA/PIDL/URL.h>
#include <Core/CCA/SSIDL/array.h>

#include <Core/Exceptions/InternalError.h>

#include <string>
#include <map>

namespace sci_cca {

using std::string;
using SSIDL::array1;

class PortInfoImpl : public PortInfo {
public:
  PortInfoImpl( const string & name,
		const string & type,
		const array1<string> & properties);
  ~PortInfoImpl();

  virtual  string  getType();
  virtual  string  getName();
  virtual  string  getProperty( const string & name );

private:
  string name_;
  string type_;
  array1<string> properties_;
};

} // namespace sci_cca

#endif 

