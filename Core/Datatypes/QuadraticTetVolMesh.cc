//  The contents of this file are subject to the University of Utah Public
//  License (the "License"); you may not use this file except in compliance
//  with the License.
//  
//  Software distributed under the License is distributed on an "AS IS"
//  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
//  License for the specific language governing rights and limitations under
//  the License.
//  
//  The Original Source Code is SCIRun, released March 12, 2001.
//  
//  The Original Source Code was developed by the University of Utah.
//  Portions created by UNIVERSITY are Copyright (C) 2001, 1994
//  University of Utah. All Rights Reserved.
//  
//    File   : QuadraticTetVolMesh.cc
//    Author : Martin Cole
//    Date   : Sun Feb 24 14:38:20 2002


#include <Core/Datatypes/QuadraticTetVolMesh.h>
#include <Core/Malloc/Allocator.h>
#include <iostream>
#include <algorithm>

namespace SCIRun {

QuadraticTetVolMesh::QuadraticTetVolMesh() :
  TetVolMesh()
{
}

QuadraticTetVolMesh::QuadraticTetVolMesh(const QuadraticTetVolMesh& copy) :
  TetVolMesh(copy)
{
}

QuadraticTetVolMesh::QuadraticTetVolMesh(const TetVolMesh &tv) :
  TetVolMesh(tv)
{
}

QuadraticTetVolMesh::~QuadraticTetVolMesh()
{
}

void
QuadraticTetVolMesh::begin(Node::iterator &itr) const
{
  itr = 0;
}

void
QuadraticTetVolMesh::end(Node::iterator &itr) const
{
  itr = points_.size() + cells_.size() * 6;
}

void
QuadraticTetVolMesh::size(Node::size_type &s) const
{
  s = points_.size() + cells_.size() * 6;
}

const Point &
QuadraticTetVolMesh::ave_point(const Point &p0, const Point &p1) const
{
  static Point rval;
  rval = Point(p0.x() + p1.x() * 0.5, 
	       p0.y() + p1.y() * 0.5, 
	       p0.z() + p1.z() * 0.5);
  return rval;
}

void
QuadraticTetVolMesh::get_nodes(Node::array_type &array, 
			       Cell::index_type idx) const
{
  TetVolMesh::get_nodes(array, idx);
  const int sz = points_.size();
  const int cindex = idx * 6;
  array.push_back(sz + cindex);
  array.push_back(sz + cindex + 1);
  array.push_back(sz + cindex + 2);
  array.push_back(sz + cindex + 3);
  array.push_back(sz + cindex + 4);
  array.push_back(sz + cindex + 5);
}

void 
QuadraticTetVolMesh::get_point(Point &result, Node::index_type index) const
{ 
  // nodes 1 - 4 are the natural tet extremities.
  // nodes 5 - 10 are center points of edges bewteen the 4 tet nodes
  const int sz = points_.size();
  if (index < sz) {
    TetVolMesh::get_point(result, index);
  } else {
    int cell    = (index - sz) / 6;
    int node_id = (index - sz) % 6;
    Point p0, p1;
    switch (node_id) {
    case 0:
      // maps to node 5 center between nodes 1, 2
      TetVolMesh::get_point(p0, cells_[cell * 4]);
      TetVolMesh::get_point(p1, cells_[cell * 4 + 1]);
      break;
    case 1:
      // maps to node 6 center between nodes 1, 3
      TetVolMesh::get_point(p0, cells_[cell * 4]);
      TetVolMesh::get_point(p1, cells_[cell * 4 + 2]);
      break;
    case 2:
      // maps to node 7 center between nodes 1, 4
      TetVolMesh::get_point(p0, cells_[cell * 4]);
      TetVolMesh::get_point(p1, cells_[cell * 4 + 3]);
      break;
    case 3:
      // maps to node 8 center between nodes 2, 3
      TetVolMesh::get_point(p0, cells_[cell * 4 + 1]);
      TetVolMesh::get_point(p1, cells_[cell * 4 + 2]);
      break;
    case 4:
      // maps to node 9 center between nodes 3, 4
      TetVolMesh::get_point(p0, cells_[cell * 4 + 2]);
      TetVolMesh::get_point(p1, cells_[cell * 4 + 3]);
      break;
    case 5:
      // maps to node 10 center between nodes 2, 2
      TetVolMesh::get_point(p0, cells_[cell * 4 + 1]);
      TetVolMesh::get_point(p1, cells_[cell * 4 + 3]);
      break;
    };
    result = ave_point(p0, p1);
  }
}

void
QuadraticTetVolMesh::get_weights(const Point& p, Node::array_type &l,
		   vector<double> &w)
{
  Cell::index_type idx;
  if (locate(idx, p))
  {
    Node::array_type ra(10);
    get_nodes(ra,idx);
    Point p1, p2, p3, p4;
    TetVolMesh::get_point(p1, ra[0]);
    TetVolMesh::get_point(p2, ra[1]);
    TetVolMesh::get_point(p3, ra[2]);
    TetVolMesh::get_point(p4, ra[3]);
    double x1=p1.x();
    double y1=p1.y();
    double z1=p1.z();
    double x2=p2.x();
    double y2=p2.y();
    double z2=p2.z();
    double x3=p3.x();
    double y3=p3.y();
    double z3=p3.z();
    double x4=p4.x();
    double y4=p4.y();
    double z4=p4.z();

    double a1=+x2*(y3*z4-y4*z3)+x3*(y4*z2-y2*z4)+x4*(y2*z3-y3*z2);
    double a2=-x3*(y4*z1-y1*z4)-x4*(y1*z3-y3*z1)-x1*(y3*z4-y4*z3);
    double a3=+x4*(y1*z2-y2*z1)+x1*(y2*z4-y4*z2)+x2*(y4*z1-y1*z4);
    double a4=-x1*(y2*z3-y3*z2)-x2*(y3*z1-y1*z3)-x3*(y1*z2-y2*z1);
    double iV6=1.0/(a1+a2+a3+a4);
  
    double b2=+(y4*z1-y1*z4)+(y1*z3-y3*z1)+(y3*z4-y4*z3);
    double c2=-(x4*z1-x1*z4)-(x1*z3-x3*z1)-(x3*z4-x4*z3);
    double d2=+(x4*y1-x1*y4)+(x1*y3-x3*y1)+(x3*y4-x4*y3);

    double b3=-(y1*z2-y2*z1)-(y2*z4-y4*z2)-(y4*z1-y1*z4);
    double c3=+(x1*z2-x2*z1)+(x2*z4-x4*z2)+(x4*z1-x1*z4);
    double d3=-(x1*y2-x2*y1)-(x2*y4-x4*y2)-(x4*y1-x1*y4);

    double b4=+(y2*z3-y3*z2)+(y3*z1-y1*z3)+(y1*z2-y2*z1);
    double c4=-(x2*z3-x3*z2)-(x3*z1-x1*z3)-(x1*z2-x2*z1);
    double d4=+(x2*y3-x3*y2)+(x3*y1-x1*y3)+(x1*y2-x2*y1);

    //map real to parent nodes
    double xi = iV6*(a2+b2*p.x()+c2*p.y()+d2*p.z());
    double nu = iV6*(a3+b3*p.x()+c3*p.y()+d3*p.z());
    double gam = iV6*(a4+b4*p.x()+c4*p.y()+d4*p.z());

    l.push_back(ra[0]);
    l.push_back(ra[1]);
    l.push_back(ra[2]);
    l.push_back(ra[3]);
    l.push_back(ra[4]);
    l.push_back(ra[5]);
    l.push_back(ra[6]);
    l.push_back(ra[7]);
    l.push_back(ra[8]);
    l.push_back(ra[9]);

    w.resize(10);
    w.push_back(-1*(1-xi-nu-gam)*(1-2*(1-xi-nu-gam)));
    w.push_back(-1*xi*(1-2*xi));
    w.push_back(-1*nu*(1-2*nu));
    w.push_back(-1*gam*(1-2*gam));
    w.push_back(4*xi*(1-xi-nu-gam));
    w.push_back(4*nu*(1-xi-nu-gam));
    w.push_back(4*gam*(1-xi-nu-gam));
    w.push_back(4*xi*nu);
    w.push_back(4*nu*gam);
    w.push_back(4*xi*gam);
  }
}

  // get grad relative to point p
void
QuadraticTetVolMesh::get_gradient_basis(Cell::index_type ci, const Point& p,
					Vector& g0, Vector& g1, Vector& g2, 
					Vector& g3, Vector& g4, Vector& g5, 
					Vector& g6, Vector& g7, Vector& g8, 
					Vector& g9) const
{
  Vector dxi, dnu, dgam;
  ///////////////////////start convert to xi,nu,gam       
  Point p1, p2, p3, p4;
  TetVolMesh::get_point(p1, cells_[ci * 4]);
  TetVolMesh::get_point(p2, cells_[ci * 4 + 1]);
  TetVolMesh::get_point(p3, cells_[ci * 4 + 2]);
  TetVolMesh::get_point(p4, cells_[ci * 4 + 3]);

  double x1=p1.x();
  double y1=p1.y();
  double z1=p1.z();
  double x2=p2.x();
  double y2=p2.y();
  double z2=p2.z();
  double x3=p3.x();
  double y3=p3.y();
  double z3=p3.z();
  double x4=p4.x();
  double y4=p4.y();
  double z4=p4.z();
  
  double a1=+x2*(y3*z4-y4*z3)+x3*(y4*z2-y2*z4)+x4*(y2*z3-y3*z2);
  double a2=-x3*(y4*z1-y1*z4)-x4*(y1*z3-y3*z1)-x1*(y3*z4-y4*z3);
  double a3=+x4*(y1*z2-y2*z1)+x1*(y2*z4-y4*z2)+x2*(y4*z1-y1*z4);
  double a4=-x1*(y2*z3-y3*z2)-x2*(y3*z1-y1*z3)-x3*(y1*z2-y2*z1);
  double iV6=1.0/(a1+a2+a3+a4);
  
  double b2=+(y4*z1-y1*z4)+(y1*z3-y3*z1)+(y3*z4-y4*z3);
  double c2=-(x4*z1-x1*z4)-(x1*z3-x3*z1)-(x3*z4-x4*z3);
  double d2=+(x4*y1-x1*y4)+(x1*y3-x3*y1)+(x3*y4-x4*y3);
  
  double b3=-(y1*z2-y2*z1)-(y2*z4-y4*z2)-(y4*z1-y1*z4);
  double c3=+(x1*z2-x2*z1)+(x2*z4-x4*z2)+(x4*z1-x1*z4);
  double d3=-(x1*y2-x2*y1)-(x2*y4-x4*y2)-(x4*y1-x1*y4);
  
  double b4=+(y2*z3-y3*z2)+(y3*z1-y1*z3)+(y1*z2-y2*z1);
  double c4=-(x2*z3-x3*z2)-(x3*z1-x1*z3)-(x1*z2-x2*z1);
  double d4=+(x2*y3-x3*y2)+(x3*y1-x1*y3)+(x1*y2-x2*y1);
  
  //map real to parent nodes
  double xi = iV6*(a2+b2*p.x()+c2*p.y()+d2*p.z());
  double nu = iV6*(a3+b3*p.x()+c3*p.y()+d3*p.z());
  double gam = iV6*(a4+b4*p.x()+c4*p.y()+d4*p.z());

  ///////////////////////end convert to xi,nu,gam

  double jac_el = calc_jac_derivs(dxi,dnu,dgam,xi,nu,gam, ci); 
  
  if (jac_el <= 0) cerr << "ERROR: jacobian <= 0\n";
  
  for (int i=0; i< 10; i++) {

    double dphidxi = calc_dphi_dxi(i,xi,nu,gam);
    double dphidnu = calc_dphi_dnu(i,xi,nu,gam);
    double dphidgam = calc_dphi_dgam(i,xi,nu,gam);

    double value1 = (dnu.y()*dgam.z() - dnu.z()*dgam.y())*dphidxi - (dxi.y()*dgam.z() - dxi.z()*dgam.y())*dphidnu + (dxi.y()*dnu.z() - dxi.z()*dnu.y())*dphidgam;

    value1 = (1.0/jac_el)*value1;
   
    double value2 =  -1.0*(dnu.x()*dgam.z() - dnu.z()*dgam.x())*dphidxi + (dxi.x()*dgam.z() - dxi.z()*dgam.x())*dphidnu - (dxi.x()*dnu.z() - dxi.z()*dnu.x())*dphidgam;

    value2 = (1.0/jac_el)*value2;

    double value3 = (dnu.x()*dgam.y() - dnu.y()*dgam.x())*dphidxi - (dxi.x()*dgam.y() - dxi.y()*dgam.x())*dphidnu + (dxi.x()*dnu.y() - dxi.y()*dnu.x())*dphidgam;

    value3 = (1.0/jac_el)*value3;

    switch (i) {
    case 0: g0 = Vector(value1,value2,value3);
      break;
    case 1: g1 = Vector(value1,value2,value3);
      break;
    case 2: g2 = Vector(value1,value2,value3);
      break;
    case 3: g3 = Vector(value1,value2,value3);
      break;
    case 4: g4 = Vector(value1,value2,value3);
      break;
    case 5: g5 = Vector(value1,value2,value3);
      break;
    case 6: g6 = Vector(value1,value2,value3);
      break;
    case 7: g7 = Vector(value1,value2,value3);
      break;
    case 8: g8 = Vector(value1,value2,value3);
      break;
    case 9: g9 = Vector(value1,value2,value3);
      break;
    }
  }
}

  //! gradient for gauss pts 
  double 
    QuadraticTetVolMesh::get_gradient_basis(Cell::index_type ci, int gaussPt, 
					    const Point&, Vector& g0, Vector& g1, 
					    Vector& g2, Vector& g3, Vector& g4, 
					    Vector& g5, Vector& g6, Vector& g7, 
					    Vector& g8, Vector& g9) const
{
  double xi, nu, gam;
  Vector dxi, dnu, dgam;

  switch(gaussPt) {
  case 0: 
    xi = 0.25;
    nu = 0.25;
    gam = 0.25;
    break;
  case 1:
    xi = 0.5;
    nu = 1.0/6.0;
    gam = 1.0/6.0;
    break;
  case 2:
    xi = 1.0/6.0;
    nu = 0.5;
    gam = 1.0/6.0;
    break;
  case 3:
    xi = 1.0/6.0;
    nu = 1.0/6.0;
    gam = 0.5;
    break;
  case 4:
    xi = 1.0/6.0;
    nu = 1.0/6.0;
    gam = 1.0/6.0;
    break;
  default: 
    cerr << "Error in get_gradQuad: Incorrect index for gaussPt. "
	 << "index = " << gaussPt << endl;
  }

  double jac_el = calc_jac_derivs(dxi,dnu,dgam,xi,nu,gam, ci);

  if (jac_el <= 0) cerr << "ERROR: jacobian <= 0\n";

  for (int i=0; i< 10; i++) {
    double dphidxi = calc_dphi_dxi(i,xi,nu,gam);
    double dphidnu = calc_dphi_dnu(i,xi,nu,gam);
    double dphidgam = calc_dphi_dgam(i,xi,nu,gam);

    double value1 = ((dnu.y()*dgam.z() - dnu.z()*dgam.y())*dphidxi - 
		     (dxi.y()*dgam.z() - dxi.z()*dgam.y())*dphidnu + 
		     (dxi.y()*dnu.z() - dxi.z()*dnu.y())*dphidgam);

    value1 = (1.0/jac_el)*value1;
   
    double value2 =  (-1.0*(dnu.x()*dgam.z() - dnu.z()*dgam.x())*dphidxi + 
		      (dxi.x()*dgam.z() - dxi.z()*dgam.x())*dphidnu - 
		      (dxi.x()*dnu.z() - dxi.z()*dnu.x())*dphidgam);

    value2 = (1.0/jac_el)*value2;

    double value3 = ((dnu.x()*dgam.y() - dnu.y()*dgam.x())*dphidxi - 
		     (dxi.x()*dgam.y() - dxi.y()*dgam.x())*dphidnu + 
		     (dxi.x()*dnu.y() - dxi.y()*dnu.x())*dphidgam);

    value3 = (1.0/jac_el)*value3;

    switch (i) {
    case 0: g0 = Vector(value1,value2,value3);
      break;
    case 1: g1 = Vector(value1,value2,value3);
      break;
    case 2: g2 = Vector(value1,value2,value3);
      break;
    case 3: g3 = Vector(value1,value2,value3);
      break;
    case 4: g4 = Vector(value1,value2,value3);
      break;
    case 5: g5 = Vector(value1,value2,value3);
      break;
    case 6: g6 = Vector(value1,value2,value3);
      break;
    case 7: g7 = Vector(value1,value2,value3);
      break;
    case 8: g8 = Vector(value1,value2,value3);
      break;
    case 9: g9 = Vector(value1,value2,value3);
      break;
    }
  }

  return(jac_el);
}

double 
QuadraticTetVolMesh::calc_jac_derivs(Vector &dxi, Vector &dnu, Vector &dgam, 
				     double xi, double nu, double gam, 
				     Cell::index_type ci) const 
{

  double dxdxi = 0.0;
  double dxdnu = 0.0;
  double dxdgam = 0.0;
  double dydxi = 0.0;
  double dydnu = 0.0;
  double dydgam = 0.0;
  double dzdxi = 0.0;
  double dzdnu = 0.0;
  double dzdgam = 0.0;

  Node::array_type node_array;
  get_nodes(node_array, ci);
  ASSERT(node_array.size() == 10);

  for (int k=0; k<10; k++) {
    Point p;
    get_point(p, node_array[k]);
    double deriv = calc_dphi_dxi(k,xi,nu,gam);
    const double pxd = p.x() * deriv;
    const double pyd = p.y() * deriv;
    const double pzd = p.z() * deriv;

    dxdxi += pxd;
    dydxi += pyd;
    dzdxi += pzd;

    deriv = calc_dphi_dnu(k,xi,nu,gam);

    dxdnu += pxd;
    dydnu += pyd;
    dzdnu += pzd;

    deriv = calc_dphi_dgam(k,xi,nu,gam);

    dxdgam += pxd;
    dydgam += pyd;
    dzdgam += pzd;
  }

  dxi = Vector(dxdxi, dydxi, dzdxi);
  dnu = Vector(dxdnu, dydnu, dzdnu);
  dgam = Vector(dxdgam, dydgam, dzdgam);

  return (dxdxi * (dydnu * dzdgam - dydgam * dzdnu) - 
	  dydxi * (dxdnu * dzdgam - dxdgam * dzdnu) + 
	  dzdxi * (dxdnu * dydgam - dxdgam * dydnu));

}

double 
QuadraticTetVolMesh::calc_dphi_dxi(int ptNum, double xi, double nu, 
				   double gam) const
{

  double value;

  switch (ptNum) {
  case 0: value = -3 + 4*xi + 4*nu + 4*gam;
    break;
  case 1: value = -1 + 4*xi;
    break;
  case 2: value = 0;
    break;
  case 3: value = 0;
    break;
  case 4: value = 4 - 4*gam - 8*xi - 4*nu;
    break;
  case 5: value = -4*nu;
    break;
  case 6: value = -4*gam;
    break;
  case 7: value = 4*nu;
    break;
  case 8: value = 0;
    break;
  case 9: value = 4*gam;
    break;
  default: 
    cerr << "Error in calcu_dphi_dxi: Incorrect index for shape function." 
	 << "index = " << ptNum << endl;
  }

  return(value);
}

double 
QuadraticTetVolMesh::calc_dphi_dnu(int ptNum, double xi, double nu, 
				   double gam) const
{

  double value;

  switch (ptNum) {
  case 0: value = -3 + 4*xi + 4*nu + 4*gam;
    break;
  case 1: value = 0;
    break;
  case 2: value = -1 + 4*nu;
    break;
  case 3: value = 0;
    break;
  case 4: value = -4*xi;
    break;
  case 5: value = 4 - 4*gam - 4*xi - 8*nu;
    break;
  case 6: value = -4*gam;
    break;
  case 7: value = 4*xi;
    break;
  case 8: value = 4*gam;
    break;
  case 9: value = 0;
    break;
  default: 
    cerr << "Error in calc_dphi_dnu: Incorrect index for shape function. " 
	 << "index = " << ptNum << endl;
  }
  
  return(value);
}

double 
QuadraticTetVolMesh::calc_dphi_dgam(int ptNum, double xi, double nu, 
				    double gam) const 
{
  double value;

  switch (ptNum) {
  case 0: value = -3 + 4*xi + 4*nu + 4*gam;
    break;
  case 1: value = 0;
    break;
  case 2: value = 0;
    break;
  case 3: value = -1 + 4*gam;
    break;
  case 4: value = -4*xi;
    break;
  case 5: value = -4*nu;
    break;
  case 6: value = 4 - 8*gam - 4*xi - 4*nu;
    break;
  case 7: value = 0;
    break;
  case 8: value = 4*nu;
    break;
  case 9: value = 4*xi;
    break;
  default: 
    cerr << "Error in calc_dphi_dgam: Incorrect index for shape function. "
	 << "index = " << ptNum << endl;
  }

  return(value);
}

} // end namespace SCIRun
