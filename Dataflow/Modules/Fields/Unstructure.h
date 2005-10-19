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


//    File   : Unstructure.h
//    Author : Michael Callahan
//    Date   : September 2001

#if !defined(Unstructure_h)
#define Unstructure_h

#include <Core/Util/TypeDescription.h>
#include <Core/Util/DynamicLoader.h>
#include <Core/Util/ProgressReporter.h>
#include <Core/Basis/HexTrilinearLgn.h>
#include <Core/Basis/QuadBilinearLgn.h>
#include <Core/Basis/CrvLinearLgn.h>
#include <Core/Basis/Constant.h>
#include <Core/Datatypes/LatVolMesh.h>
#include <Core/Datatypes/ImageMesh.h>
#include <Core/Datatypes/ScanlineMesh.h>

namespace SCIRun {

typedef LatVolMesh<HexTrilinearLgn<Point> >               LVMesh;
typedef ImageMesh<QuadBilinearLgn<Point> >                IMesh;
typedef ScanlineMesh<CrvLinearLgn<Point> >                SMesh;

class UnstructureAlgo : public DynamicAlgoBase
{
public:
  virtual FieldHandle execute(ProgressReporter *module, FieldHandle src) = 0;

  //! support the dynamically compiled algorithm concept
  static CompileInfoHandle get_compile_info(const TypeDescription *fsrc,
					    const string &mesh_dst,
					    const string &basis_dst,
					    const string &data_dst);
};


template <class FSRC, class FDST>
class UnstructureAlgoT : public UnstructureAlgo
{
public:
  //! virtual interface. 
  virtual FieldHandle execute(ProgressReporter *module, FieldHandle src);
};


struct SpecialUnstructuredHash
{
  size_t operator()(const LVMesh::Node::index_type &n) const
  { return n.i_ * ((1 << 20) - 1) + n.j_ * ((1 << 10) - 1) + n.k_; }

  size_t operator()(const IMesh::Node::index_type &n) const
  { return n.i_ * ((1 << 10) - 1) + n.j_; }

  size_t operator()(const SMesh::Node::index_type &n) const
  { return n; }

  size_t operator()(const LVMesh::Elem::index_type &n) const
  { return n.i_ * ((1 << 20) - 1) + n.j_ * ((1 << 10) - 1) + n.k_; }

  size_t operator()(const IMesh::Elem::index_type &n) const
  { return n.i_ * ((1 << 10) - 1) + n.j_; }

  size_t operator()(const SMesh::Elem::index_type &n) const
  { return n; }
};

#ifdef HAVE_HASH_MAP
struct SpecialUnstructuredEqual
{
  bool operator()(const LVMesh::Node::index_type &a,
		  const LVMesh::Node::index_type &b) const
  { return a.i_ == b.i_ && a.j_ == b.j_ && a.k_ == b.k_; }

  bool operator()(const IMesh::Node::index_type &a,
		  const IMesh::Node::index_type &b) const
  { return a.i_ == b.i_ && a.j_ == b.j_; }

  bool operator()(const SMesh::Node::index_type &a,
		  const SMesh::Node::index_type &b) const
  { return a == b; }

  bool operator()(const LVMesh::Elem::index_type &a,
		  const LVMesh::Elem::index_type &b) const
  { return a.i_ == b.i_ && a.j_ == b.j_ && a.k_ == b.k_; }

  bool operator()(const IMesh::Elem::index_type &a,
		  const IMesh::Elem::index_type &b) const
  { return a.i_ == b.i_ && a.j_ == b.j_; }

  bool operator()(const SMesh::Elem::index_type &a,
		  const SMesh::Elem::index_type &b) const
  { return a == b; }
};
#else
struct SpecialUnstructuredLess
{
  bool operator()(const LVMesh::Node::index_type &a,
		  const LVMesh::Node::index_type &b) const
  { return a.i_ < b.i_ || a.i_ == b.i_ && ( a.j_ < b.j_ || a.j_ == b.j_ && a.k_ < b.k_); }

  bool operator()(const IMesh::Node::index_type &a,
		  const IMesh::Node::index_type &b) const
  { return a.i_ < b.i_ || a.i_ == b.i_ && a.j_ < b.j_; }

  bool operator()(const SMesh::Node::index_type &a,
		  const SMesh::Node::index_type &b) const
  { return a < b; }

  bool operator()(const LVMesh::Elem::index_type &a,
		  const LVMesh::Elem::index_type &b) const
  { return a.i_ < b.i_ || a.i_ == b.i_ && ( a.j_ < b.j_ || a.j_ == b.j_ && a.k_ < b.k_); }

  bool operator()(const IMesh::Elem::index_type &a,
		  const IMesh::Elem::index_type &b) const
  { return a.i_ < b.i_ || a.i_ == b.i_ && a.j_ < b.j_; }

  bool operator()(const SMesh::Elem::index_type &a,
		  const SMesh::Elem::index_type &b) const
  { return a < b; }
};
#endif


template <class FSRC, class FDST>
FieldHandle
UnstructureAlgoT<FSRC, FDST>::execute(ProgressReporter *module,
				      FieldHandle field_h)
{
  FSRC *ifield = dynamic_cast<FSRC *>(field_h.get_rep());
  typename FSRC::mesh_handle_type mesh = ifield->get_typed_mesh();

  typename FDST::mesh_handle_type outmesh = scinew typename FDST::mesh_type();

  bool pointCloud =
    (outmesh->get_type_description()->get_name() == 
     "PointCloudMesh<ConstantBasis<Point> >");

#ifdef HAVE_HASH_MAP
  typedef hash_map<typename FSRC::mesh_type::Node::index_type,
    typename FDST::mesh_type::Node::index_type,
    SpecialUnstructuredHash, SpecialUnstructuredEqual> node_hash_type;

  typedef hash_map<typename FSRC::mesh_type::Elem::index_type,
    typename FDST::mesh_type::Elem::index_type,
    SpecialUnstructuredHash, SpecialUnstructuredEqual> elem_hash_type;
#else
  typedef map<typename FSRC::mesh_type::Node::index_type,
    typename FDST::mesh_type::Node::index_type,
    SpecialUnstructuredLess> node_hash_type;

  typedef map<typename FSRC::mesh_type::Elem::index_type,
    typename FDST::mesh_type::Elem::index_type,
    SpecialUnstructuredLess> elem_hash_type;
#endif
  
  node_hash_type nodemap;
  elem_hash_type elemmap;
  mesh->synchronize(Mesh::ALL_ELEMENTS_E);

  typename FSRC::mesh_type::Node::iterator bn, en;
  mesh->begin(bn); mesh->end(en);
  while (bn != en)
  {
    ASSERT(nodemap.find(*bn) == nodemap.end());
    Point np;
    mesh->get_center(np, *bn);
    nodemap[*bn] = outmesh->add_point(np);
    ++bn;
  }

  // Point clouds do not have elements so skip.
  if( !pointCloud ) {
    typename FSRC::mesh_type::Elem::iterator bi, ei;
    mesh->begin(bi); mesh->end(ei);
    while (bi != ei) {
      // Add this element to the new mesh.
      typename FSRC::mesh_type::Node::array_type onodes;
      mesh->get_nodes(onodes, *bi);
      typename FDST::mesh_type::Node::array_type nnodes(onodes.size());

      for (unsigned int i=0; i<onodes.size(); i++)
	{
	  ASSERT(nodemap.find(onodes[i]) != nodemap.end());
	  nnodes[i] = nodemap[onodes[i]];
	}

      elemmap[*bi] = outmesh->add_elem(nnodes);
      ++bi;
    }
  }

  // really should copy normals
  outmesh->synchronize(Mesh::NORMALS_E);

  FDST *ofield = scinew FDST(outmesh);

  if (field_h->basis_order() == 1)
  {
    typename node_hash_type::iterator hitr = nodemap.begin();

    while (hitr != nodemap.end())
    {
      typename FSRC::value_type val;
      ifield->value(val, (*hitr).first);
      ofield->set_value(val, (*hitr).second);

      ++hitr;
    }
  }

  // Point clouds do not have elements so skip.
  else if (!pointCloud &&
	   ifield->order_type_description()->get_name() ==
	   ofield->order_type_description()->get_name())
  {
    typename elem_hash_type::iterator hitr = elemmap.begin();
    
    while (hitr != elemmap.end())
    {
      typename FSRC::value_type val;
      ifield->value(val, (*hitr).first);
      ofield->set_value(val, (*hitr).second);

      ++hitr;
    }
  }
  else
  {
    module->warning("Unable to copy data at this field data location.");
  }

  return ofield;
}


} // end namespace SCIRun

#endif // Unstructure_h
