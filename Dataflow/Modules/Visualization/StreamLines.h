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

//    File   : StreamLines.h
//    Author : Michael Callahan
//    Date   : July 2001

#if !defined(StreamLines_h)
#define StreamLines_h

#include <Core/Disclosure/TypeDescription.h>
#include <Core/Disclosure/DynamicLoader.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/CurveField.h>
#include <algorithm>

namespace SCIRun {

class StreamLinesAlgo : public DynamicAlgoBase
{
public:
  virtual FieldHandle execute(MeshHandle seed_mesh_h,
			      VectorFieldInterface *vfi,
			      double tolerance,
			      double stepsize,
			      int maxsteps,
			      int direction,
			      bool remove_colinear_p) = 0;


  //! support the dynamically compiled algorithm concept
  static CompileInfo *get_compile_info(const TypeDescription *smesh,
				       const TypeDescription *sloc);
protected:
  //bool interpolate(VectorFieldInterface *vfi, const Point &p, Vector &v);

  //! This particular implementation uses Runge-Kutta-Fehlberg.
  void FindStreamLineNodes(vector<Point>&, Point, double, double, int, 
  			   VectorFieldInterface *, bool remove_colinear_p);

  //! Compute the inner terms of the RKF formula.
  //bool ComputeRKFTerms(vector<Vector> &, const Point&, double,
  //		       VectorFieldInterface *);

};


template <class SFLD, class SLOC>
class StreamLinesAlgoT : public StreamLinesAlgo
{
public:
  //! virtual interface. 
  virtual FieldHandle execute(MeshHandle seed_mesh_h,
			      VectorFieldInterface *vfi,
			      double tolerance,
			      double stepsize,
			      int maxsteps,
			      int direction,
			      bool remove_colinear_p);
};



template <class SMESH, class SLOC>
FieldHandle
StreamLinesAlgoT<SMESH, SLOC>::execute(MeshHandle seed_mesh_h,
				       VectorFieldInterface *vfi,
				       double tolerance,
				       double stepsize,
				       int maxsteps,
				       int direction,
				       bool rcp)
{
  SMESH *smesh = dynamic_cast<SMESH *>(seed_mesh_h.get_rep());

  const double tolerance2 = tolerance * tolerance;

  CurveMeshHandle cmesh = scinew CurveMesh();
  CurveField<double> *cf = scinew CurveField<double>(cmesh, Field::NODE);

  Point seed;
  Vector test;
  vector<Point> nodes;
  vector<Point>::iterator node_iter;
  CurveMesh::Node::index_type n1, n2;

  // Try to find the streamline for each seed point.
  typename SLOC::iterator seed_iter, seed_iter_end;
  smesh->begin(seed_iter);
  smesh->end(seed_iter_end);
  int cfsize = 0;
  while (seed_iter != seed_iter_end)
  {
    smesh->get_point(seed, *seed_iter);

    // Is the seed point inside the field?
    if (!vfi->interpolate(test, seed))
    {
      //cout << "StreamLines: WARNING: seed point was not inside the field.\n";
      ++seed_iter;
      continue;
    }

    nodes.clear();
    nodes.push_back(seed);

    // Find the negative streamlines.
    if( direction <= 1 ) {
      FindStreamLineNodes(nodes, seed, tolerance2, -stepsize,
			  maxsteps, vfi, rcp);
      std::reverse(nodes.begin(), nodes.end());
    }
    // Append the positive streamlines.
    if( direction >= 1 )
      FindStreamLineNodes(nodes, seed, tolerance2, stepsize,
			  maxsteps, vfi, rcp);

    node_iter = nodes.begin();
    if (node_iter != nodes.end())
    {
      n1 = cmesh->add_node(*node_iter);
      cf->resize_fdata();
      cf->set_value((double)(*seed_iter), n1);
      cfsize++;
      ++node_iter;
      while (node_iter != nodes.end())
      {
	n2 = cmesh->add_node(*node_iter);
	cf->resize_fdata();
	cf->set_value((double)(*seed_iter), n2);
	cfsize++;
	cmesh->add_edge(n1, n2);
	n1 = n2;
	++node_iter;
      }
    }

    ++seed_iter;
  }

  cf->freeze();

  if (cfsize == 0)
  {
    delete cf;
    return 0;
  }
  else
  {
    return cf;
  }
}


} // end namespace SCIRun

#endif // StreamLines_h
