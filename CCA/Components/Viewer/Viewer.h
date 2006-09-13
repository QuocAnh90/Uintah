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
 *  Viewer.h
 *
 *  Written by:
 *   Keming Zhang
 *   Department of Computer Science
 *   University of Utah
 *   May 2002
 *
 */

#ifndef CCA_Components_Viewer_h
#define CCA_Components_Viewer_h

#include <Core/CCA/spec/cca_sidl.h>

namespace Viewer {

class Viewer;
class MainWindow;

class ViewPort : public virtual sci::cca::ports::ViewPort {
public:
  virtual ~ViewPort() {}
  virtual int view2dPDE(const SSIDL::array1<double> &nodes,
                        const SSIDL::array1<int> &triangles,
                        const SSIDL::array1<double> &solution);

private:
  MainWindow* mw;
};

// rename to 2D Viewer???

class Viewer : public sci::cca::Component {
public:
  Viewer();
  virtual ~Viewer();
  virtual void setServices(const sci::cca::Services::pointer& svc);

private:
  Viewer(const Viewer&);
  Viewer& operator=(const Viewer&);
  sci::cca::Services::pointer services;
};


}

#endif
