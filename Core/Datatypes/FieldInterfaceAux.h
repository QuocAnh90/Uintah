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

//
//  Written by:
//   Michael Callahan
//   Department of Computer Science
//   University of Utah
//   May 2002
//
//  Copyright (C) 2002 SCI Institute
//
//
//
// This is the templated implementations of the FieldInterface classes.
// It should not need to be widely included.
//


#ifndef Datatypes_FieldInterfaceAux_h
#define Datatypes_FieldInterfaceAux_h

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Tensor.h>
#include <Core/Util/DynamicLoader.h>

namespace SCIRun {

//! Should only be instantiated for fields with scalar data.
template <class F>
class SFInterface : public ScalarFieldInterface {
public:
  SFInterface(const F *fld) :
    fld_(fld)
  {}
  
  virtual bool compute_min_max(double &minout, double &maxout) const;
  virtual bool interpolate(double &result, const Point &p) const;
  virtual bool interpolate_many(vector<double> &results,
				const vector<Point> &points) const;
  virtual double find_closest(double &result, const Point &p) const;
private:

  bool finterpolate(double &result, const Point &p) const;

  const F        *fld_;
};



class ScalarFieldInterfaceMaker : public DynamicAlgoBase
{
public:
  virtual ScalarFieldInterface *make(const Field *field) = 0;
  static CompileInfo *get_compile_info(const TypeDescription *ftd);
};



template <class F>
class SFInterfaceMaker : public ScalarFieldInterfaceMaker
{
public:

  virtual ScalarFieldInterface *make(const Field *field)
  {
    const F *tfield = dynamic_cast<const F *>(field);
    return scinew SFInterface<F>(tfield);
  }
};



template <class F>
bool
SFInterface<F>::finterpolate(double &result, const Point &p) const
{
  typename F::mesh_handle_type mesh = fld_->get_typed_mesh();
  switch(fld_->data_at())
  {
  case F::NODE:
    {
      typename F::mesh_type::Node::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::EDGE:
    {
      typename F::mesh_type::Edge::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::FACE:
    {
      typename F::mesh_type::Face::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::CELL:
    {
      typename F::mesh_type::Cell::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::NONE:
    return false;
  }
  return true;
}


template <class F>
bool
SFInterface<F>::interpolate(double &result, const Point &p) const
{
  return finterpolate(result, p);
}


template <class F>
bool
SFInterface<F>::interpolate_many(vector<double> &results,
				 const vector<Point> &points) const
{
  bool all_interped_p = true;
  results.resize(points.size());
  unsigned int i;
  for (i=0; i < points.size(); i++)
  {
    all_interped_p &=  interpolate(results[i], points[i]);
  }
  return all_interped_p;
}


template <class F>
bool
SFInterface<F>::compute_min_max(double &minout, double &maxout) const
{
  bool result = false;
  minout = 1.0e6;
  maxout = -1.0e6;
  typename F::mesh_handle_type mesh = fld_->get_typed_mesh();
  switch (fld_->data_at())
  {
  case F::NODE:
    {
      typename F::mesh_type::Node::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Node::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	typename F::value_type val;
	if (fld_->value(val, *bi))
	{
	  if (val < minout) minout = val;
	  if (val > maxout) maxout = val;
	  result = true;
	}
	++bi;
      }      
    }
    break;

  case F::EDGE:
    {
      typename F::mesh_type::Edge::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Edge::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	typename F::value_type val;
	if (fld_->value(val, *bi))
	{
	  if (val < minout) minout = val;
	  if (val > maxout) maxout = val;
	  result = true;
	}
	++bi;
      }      
    }
    break;

  case F::FACE:
    {
      typename F::mesh_type::Face::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Face::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	typename F::value_type val;
	if (fld_->value(val, *bi))
	{
	  if (val < minout) minout = val;
	  if (val > maxout) maxout = val;
	  result = true;
	}
	++bi;
      }      
    }
    break;

  case F::CELL:
    {
      typename F::mesh_type::Cell::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Cell::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	typename F::value_type val;
	if (fld_->value(val, *bi))
	{
	  if (val < minout) minout = val;
	  if (val > maxout) maxout = val;
	  result = true;
	}
	++bi;
      }      
    }
    break;
    
  case F::NONE:
    break;
  }
  return result;
}


template <class F>
double
SFInterface<F>::find_closest(double &minout, const Point &p) const
{
  double mindist = 1.0e15;
  typename F::mesh_handle_type mesh = fld_->get_typed_mesh();
  switch (fld_->data_at())
  {
  case F::NODE:
    {
      typename F::mesh_type::Node::index_type index;
      typename F::mesh_type::Node::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Node::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      typename F::value_type val;
      fld_->value(val, index);
      minout = (double)val;
    }
    break;

  case F::EDGE:
    {
      typename F::mesh_type::Edge::index_type index;
      typename F::mesh_type::Edge::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Edge::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      typename F::value_type val;
      fld_->value(val, index);
      minout = (double)val;
    }
    break;

  case F::FACE:
    {
      typename F::mesh_type::Face::index_type index;
      typename F::mesh_type::Face::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Face::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      typename F::value_type val;
      fld_->value(val, index);
      minout = (double)val;
    }
    break;

  case F::CELL:
    {
      typename F::mesh_type::Cell::index_type index;
      typename F::mesh_type::Cell::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Cell::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      typename F::value_type val;
      fld_->value(val, index);
      minout = (double)val;
    }
    break;
    
  case F::NONE:
    break;
  }

  return mindist;
}



//! Should only be instantiated for fields with scalar data.
template <class F>
class VFInterface : public VectorFieldInterface {
public:
  VFInterface(const F *fld) :
    fld_(fld)
  {}
  
  virtual bool compute_min_max(Vector &minout, Vector  &maxout) const;
  virtual bool interpolate(Vector &result, const Point &p) const;
  virtual bool interpolate_many(vector<Vector> &results,
				const vector<Point> &points) const;
  virtual double find_closest(Vector &result, const Point &p) const;

private:
  bool finterpolate(Vector &result, const Point &p) const;

  const F        *fld_;
};



class VectorFieldInterfaceMaker : public DynamicAlgoBase
{
public:
  virtual VectorFieldInterface *make(const Field *field) = 0;
  static CompileInfo *get_compile_info(const TypeDescription *ftd);
};



template <class F>
class VFInterfaceMaker : public VectorFieldInterfaceMaker
{
public:

  virtual VectorFieldInterface *make(const Field *field)
  {
    const F *tfield = dynamic_cast<const F *>(field);
    return scinew VFInterface<F>(tfield);
  }
};



template <class F>
bool
VFInterface<F>::finterpolate(Vector  &result, const Point &p) const
{
  typename F::mesh_handle_type mesh = fld_->get_typed_mesh();
  switch(fld_->data_at())
  {
  case F::NODE:
    {
      typename F::mesh_type::Node::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::EDGE:
    {
      typename F::mesh_type::Edge::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::FACE:
    {
      typename F::mesh_type::Face::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::CELL:
    {
      typename F::mesh_type::Cell::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::NONE:
    return false;
  }
  return true;
}


template <class F>
bool
VFInterface<F>::interpolate(Vector &result, const Point &p) const
{
  return finterpolate(result, p);
}


template <class F>
bool
VFInterface<F>::interpolate_many(vector<Vector> &results,
				 const vector<Point> &points) const
{
  bool all_interped_p = true;
  results.resize(points.size());
  unsigned int i;
  for (i=0; i < points.size(); i++)
  {
    all_interped_p &=  interpolate(results[i], points[i]);
  }
  return all_interped_p;
}


template <class F>
bool
VFInterface<F>::compute_min_max(Vector  &minout, Vector  &maxout) const
{
  static const Vector MaxVector(1.0e6, 1.0e6, 1.0e6);
  static const Vector MinVector(-1.0e6, -1.0e6, -1.0e6);

  bool result = false;
  minout = MaxVector;
  maxout = MinVector;
  typename F::mesh_handle_type mesh = fld_->get_typed_mesh();
  switch (fld_->data_at())
  {
  case F::NODE:
    {
      typename F::mesh_type::Node::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Node::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	typename F::value_type val;
	if (fld_->value(val, *bi))
	{
	  if (val.x() < minout.x()) minout.x(val.x());
	  if (val.y() < minout.y()) minout.y(val.y());
	  if (val.z() < minout.z()) minout.z(val.z());

	  if (val.x() > maxout.x()) maxout.x(val.x());
	  if (val.y() > maxout.y()) maxout.y(val.y());
	  if (val.z() > maxout.z()) maxout.z(val.z());
	  result = true;
	}
	++bi;
      }      
    }
    break;

  case F::EDGE:
    {
      typename F::mesh_type::Edge::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Edge::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	typename F::value_type val;
	if (fld_->value(val, *bi))
	{
	  if (val.x() < minout.x()) minout.x(val.x());
	  if (val.y() < minout.y()) minout.y(val.y());
	  if (val.z() < minout.z()) minout.z(val.z());

	  if (val.x() > maxout.x()) maxout.x(val.x());
	  if (val.y() > maxout.y()) maxout.y(val.y());
	  if (val.z() > maxout.z()) maxout.z(val.z());
	  result = true;
	}
	++bi;
      }      
    }
    break;

  case F::FACE:
    {
      typename F::mesh_type::Face::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Face::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	typename F::value_type val;
	if (fld_->value(val, *bi))
	{
	  if (val.x() < minout.x()) minout.x(val.x());
	  if (val.y() < minout.y()) minout.y(val.y());
	  if (val.z() < minout.z()) minout.z(val.z());

	  if (val.x() > maxout.x()) maxout.x(val.x());
	  if (val.y() > maxout.y()) maxout.y(val.y());
	  if (val.z() > maxout.z()) maxout.z(val.z());
	  result = true;
	}
	++bi;
      }      
    }
    break;

  case F::CELL:
    {
      typename F::mesh_type::Cell::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Cell::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	typename F::value_type val;
	if (fld_->value(val, *bi))
	{
	  if (val.x() < minout.x()) minout.x(val.x());
	  if (val.y() < minout.y()) minout.y(val.y());
	  if (val.z() < minout.z()) minout.z(val.z());

	  if (val.x() > maxout.x()) maxout.x(val.x());
	  if (val.y() > maxout.y()) maxout.y(val.y());
	  if (val.z() > maxout.z()) maxout.z(val.z());
	  result = true;
	}
	++bi;
      }      
    }
    break;
    
  case F::NONE:
    break;
  }
  return result;
}


template <class F>
double
VFInterface<F>::find_closest(Vector &minout, const Point &p) const
{
  double mindist = 1.0e15;
  typename F::mesh_handle_type mesh = fld_->get_typed_mesh();
  switch (fld_->data_at())
  {
  case F::NODE:
    {
      typename F::mesh_type::Node::index_type index;
      typename F::mesh_type::Node::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Node::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      fld_->value(minout, index);
    }
    break;

  case F::EDGE:
    {
      typename F::mesh_type::Edge::index_type index;
      typename F::mesh_type::Edge::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Edge::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      fld_->value(minout, index);
    }
    break;

  case F::FACE:
    {
      typename F::mesh_type::Face::index_type index;
      typename F::mesh_type::Face::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Face::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      fld_->value(minout, index);
    }
    break;

  case F::CELL:
    {
      typename F::mesh_type::Cell::index_type index;
      typename F::mesh_type::Cell::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Cell::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      fld_->value(minout, index);
    }
    break;
    
  case F::NONE:
    break;
  }

  return mindist;
}



//! Should only be instantiated for fields with scalar data.
template <class F>
class TFInterface : public TensorFieldInterface {
public:
  TFInterface(const F *fld) :
    fld_(fld)
  {}

  virtual bool interpolate(Tensor &result, const Point &p) const;
  virtual bool interpolate_many(vector<Tensor> &results,
				const vector<Point> &points) const;
  virtual double find_closest(Tensor &result, const Point &p) const;
  
private:
  bool finterpolate(Tensor &result, const Point &p) const;

  const F        *fld_;
};



class TensorFieldInterfaceMaker : public DynamicAlgoBase
{
public:
  virtual TensorFieldInterface *make(const Field *field) = 0;
  static CompileInfo *get_compile_info(const TypeDescription *ftd);
};



template <class F>
class TFInterfaceMaker : public TensorFieldInterfaceMaker
{
public:

  virtual TensorFieldInterface *make(const Field *field)
  {
    const F *tfield = dynamic_cast<const Field *>(field);
    return scinew TFInterface<F>(tfield);
  }
};



template <class F>
bool
TFInterface<F>::finterpolate(Tensor &result, const Point &p) const
{
  typename F::mesh_handle_type mesh = fld_->get_typed_mesh();
  switch(fld_->data_at())
  {
  case F::NODE:
    {
      typename F::mesh_type::Node::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::EDGE:
    {
      typename F::mesh_type::Edge::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::FACE:
    {
      typename F::mesh_type::Face::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::CELL:
    {
      typename F::mesh_type::Cell::array_type locs;
      vector<double> weights;
      mesh->get_weights(p, locs, weights);

      // weights is empty if point not found.
      if (weights.size() <= 0) return false;

      result = 0;
      for (unsigned int i = 0; i < locs.size(); i++)
      {
	typename F::value_type tmp;
	if (fld_->value(tmp, locs[i]))
	{
	  result += tmp * weights[i];
	}
      }
    }
    break;

  case F::NONE:
    return false;
  }
  return true;
}



template <class F>
bool
TFInterface<F>::interpolate(Tensor &result, const Point &p) const
{
  return finterpolate(result, p);
}


template <class F>
bool
TFInterface<F>::interpolate_many(vector<Tensor> &results,
				 const vector<Point> &points) const
{
  bool all_interped_p = true;
  results.resize(points.size());
  unsigned int i;
  for (i=0; i < points.size(); i++)
  {
    all_interped_p &=  interpolate(results[i], points[i]);
  }
  return all_interped_p;
}


template <class F>
double
TFInterface<F>::find_closest(Tensor &minout, const Point &p) const
{
  double mindist = 1.0e15;
  typename F::mesh_handle_type mesh = fld_->get_typed_mesh();
  switch (fld_->data_at())
  {
  case F::NODE:
    {
      typename F::mesh_type::Node::index_type index;
      typename F::mesh_type::Node::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Node::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      fld_->value(minout, index);
    }
    break;

  case F::EDGE:
    {
      typename F::mesh_type::Edge::index_type index;
      typename F::mesh_type::Edge::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Edge::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      fld_->value(minout, index);
    }
    break;

  case F::FACE:
    {
      typename F::mesh_type::Face::index_type index;
      typename F::mesh_type::Face::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Face::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      fld_->value(minout, index);
    }
    break;

  case F::CELL:
    {
      typename F::mesh_type::Cell::index_type index;
      typename F::mesh_type::Cell::iterator bi; mesh->begin(bi);
      typename F::mesh_type::Cell::iterator ei; mesh->end(ei);
      while (bi != ei)
      {
	Point c;
	mesh->get_center(c, *bi);
	const double dist = (p - c).length2();
	if (dist < mindist)
	{
	  mindist = dist;
	  index = *bi;
	}
	++bi;
      }
      fld_->value(minout, index);
    }
    break;
    
  case F::NONE:
    break;
  }

  return mindist;
}


} // end namespace SCIRun


#endif // Datatypes_FieldInterface_h


