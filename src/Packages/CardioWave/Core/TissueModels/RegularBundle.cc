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

#include <Core/Basis/Constant.h>
#include <Core/Basis/HexTrilinearLgn.h>
#include <Core/Datatypes/StructHexVolMesh.h>
#include <Core/Containers/FData.h>
#include <Core/Datatypes/GenericField.h>

#include <Packages/CardioWave/Core/TissueModels/RegularBundle.h>
#include <vector>

namespace CardioWave {

using namespace SCIRun;

bool TissueModel_RegularBundle::create_mesh(FieldHandle& output)
{
  int numnodes_x = (numelems_x_ics_ + 2*numelems_x_ecs_)*numcellsx_ + 1;
  int numnodes_y = (numelems_y_ics_ + 2*numelems_y_ecs_)*numcellsy_ + 1;
  int numnodes_z = numelems_z_*numcellsz_ + numelems_bath_start_ + numelems_bath_end_ + 1;
  
  double dz = cell_length_/numelems_z_;
  double dx_ics = sqrt(cell_crosssection_)/numelems_x_ics_;
  double dy_ics = sqrt(cell_crosssection_)/numelems_y_ics_;
  double dx_ecs = (sqrt(cell_crosssection_/ics_vol_frac_) - sqrt(cell_crosssection_))/(2*numelems_x_ecs_);
  double dy_ecs = (sqrt(cell_crosssection_/ics_vol_frac_) - sqrt(cell_crosssection_))/(2*numelems_y_ecs_);


  StructHexVolMesh<HexTrilinearLgn<Point> >* omesh = scinew StructHexVolMesh<HexTrilinearLgn<Point> >(numnodes_x,numnodes_y,numnodes_z);
  MeshHandle mesh = dynamic_cast<Mesh *>(omesh);

  if (mesh.get_rep() == 0)
  {
    error("TissueModel_RegularBundle: Could not obtain memory for mesh");
    return (false);
  }

  GenericField<StructHexVolMesh<HexTrilinearLgn<Point> >,ConstantBasis<int>,FData3d<int,StructHexVolMesh<HexTrilinearLgn<Point> > > >* ofield =
     scinew GenericField<StructHexVolMesh<HexTrilinearLgn<Point> >,ConstantBasis<int>,FData3d<int,StructHexVolMesh<HexTrilinearLgn<Point> > > >(omesh);

  output = dynamic_cast<Field*>(ofield);
  if (output.get_rep() == 0)
  {
    error("TissueModel_RegularBundle: Could not obtain memory for field");
    return (false);
  }


  std::vector<double> x(numnodes_x);
  std::vector<int> xe(numnodes_x-1);

  x[0] = 0.0;
  int m = 1; int n = 1;
  for (int p=0; p<numelems_x_ecs_; p++, m++) {x[m] = x[m-1]+dx_ecs; xe[m-1] = 0; }
  for (int r=0; r<(numcellsx_-1); r++, n++)
  {
    for (int p=0; p<numelems_x_ics_;p++, m++) { x[m] = x[m-1]+dx_ics; xe[m-1] = n; }
    for (int p=0; p<2*numelems_x_ecs_;p++, m++) { x[m] = x[m-1]+dx_ecs; xe[m-1] = 0; }
  }
  for (int p=0; p<numelems_x_ics_;p++, m++) { x[m] = x[m-1]+dx_ics; xe[m-1] = n; }
  for (int p=0; p<numelems_x_ecs_; p++, m++) { x[m] = x[m-1]+dx_ecs; xe[m-1] = 0; }

  std::vector<double> y(numnodes_y);
  std::vector<int> ye(numnodes_y-1);

  y[0] = 0.0;
  m = 1; n = 1;
  for (int p=0; p<numelems_y_ecs_; p++, m++) {y[m] = y[m-1]+dy_ecs; ye[m-1] = 0; }
  for (int r=0; r<(numcellsy_-1); r++, n++)
  {
    for (int p=0; p<numelems_y_ics_;p++, m++) { y[m] = y[m-1]+dy_ics; ye[m-1] = n; }
    for (int p=0; p<2*numelems_y_ecs_;p++, m++) { y[m] = y[m-1]+dy_ecs; ye[m-1] = 0; }
  }
  for (int p=0; p<numelems_y_ics_;p++, m++) { y[m] = y[m-1]+dy_ics; ye[m-1] = n; }
  for (int p=0; p<numelems_y_ecs_; p++, m++) { y[m] = y[m-1]+dy_ecs; ye[m-1] = 0; }

    
  for (int i=0; i<numnodes_x; i++)
  {
    for (int j=0; j<numnodes_y; j++)
    {
      for (int k=0; k<numnodes_z; k++)
      {
        omesh->set_point(Point(x[i],y[j],dz*k),LatVolMesh<HexTrilinearLgn<Point> >::Node::index_type(omesh,i,j,k));
      }
    }
  }

  for (int i=0; i<numnodes_x-1; i++)
  {
    for (int j=0; j<numnodes_y-1; j++)
    {
      for (int k=0; k<numnodes_z-1; k++)
      {
        int ze = 0;
        if (k >= numelems_bath_start_ && k < numelems_bath_start_ + (numcellsz_*(numelems_z_))) ze = (k-numelems_bath_start_)/numelems_z_ + 1;
        if (ze && xe[i] && ye[j])
        {
          ofield->set_value(xe[i]+((numcellsx_-1)*ye[j])+(numcellsx_*numcellsy_)*(ze-1),LatVolMesh<HexTrilinearLgn<Point> >::Elem::index_type(omesh,i,j,k));
        }
        else
        {
          ofield->set_value(0,LatVolMesh<HexTrilinearLgn<Point> >::Elem::index_type(omesh,i,j,k));
        }
      }
    }
  }
  
  return (true); 
}


}
