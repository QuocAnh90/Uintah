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


/*
 *  HexToTet.h:  Convert a Hex field into a Tet field using 1-5 split
 *
 *  Written by:
 *   David Weinstein
 *   University of Utah
 *   December 2002
 *
 *  Copyright (C) 1994, 2001 SCI Group
 */

#if !defined(HexToTet_h)
#define HexToTet_h

#include <Core/Util/TypeDescription.h>
#include <Core/Util/DynamicLoader.h>
#include <Core/Basis/NoData.h>
#include <Core/Basis/Constant.h>
#include <Core/Basis/TetLinearLgn.h>
#include <Core/Datatypes/TetVolMesh.h>
#include <Core/Datatypes/GenericField.h>
#include <Core/Util/Assert.h>

namespace SCIRun {

typedef TetVolMesh<TetLinearLgn<Point> > TVMesh;

class HexToTetAlgo : public DynamicAlgoBase
{
public:
  virtual bool execute(ProgressReporter *reporter,
                       FieldHandle src, FieldHandle &dst) = 0;

  //! support the dynamically compiled algorithm concept
  static CompileInfoHandle get_compile_info(const TypeDescription *data_td);
};


template <class FSRC>
class HexToTetAlgoT : public HexToTetAlgo
{
public:
  //! virtual interface. 
  virtual bool execute(ProgressReporter *reporter,
                       FieldHandle src, FieldHandle& dst);
};


template <class FSRC>
bool
HexToTetAlgoT<FSRC>::execute(ProgressReporter *reporter,
                             FieldHandle srcH, FieldHandle& dstH) 
{
  FSRC *hvfield = dynamic_cast<FSRC*>(srcH.get_rep());
  typename FSRC::mesh_type *hvmesh = hvfield->get_typed_mesh().get_rep();
  TVMesh::handle_type tvmesh = scinew TVMesh();

  typename FSRC::mesh_type::Node::size_type hnsize; hvmesh->size(hnsize);
  typename FSRC::mesh_type::Elem::size_type hesize; hvmesh->size(hesize);

  tvmesh->node_reserve((unsigned int)hnsize);

  // Copy points directly, assuming they will have the same order.
  typename FSRC::mesh_type::Node::iterator nbi, nei;
  hvmesh->begin(nbi); hvmesh->end(nei);
  while (nbi != nei)
  {
    Point p;
    hvmesh->get_center(p, *nbi);
    tvmesh->add_point(p);
    ++nbi;
  }

  hvmesh->synchronize(Mesh::FACE_NEIGHBORS_E);

  tvmesh->elem_reserve((unsigned int)hesize * 5);

  vector<typename FSRC::mesh_type::Elem::index_type> elemmap;

  vector<bool> visited(hesize, false);
  vector<bool> nodeisdiagonal(hnsize, false);

  typename FSRC::mesh_type::Elem::iterator bi, ei;
  hvmesh->begin(bi); hvmesh->end(ei);

  const unsigned int surfsize = (unsigned int)pow(hesize, 2.0 / 3.0);
  vector<typename FSRC::mesh_type::Elem::index_type> buffers[2];
  buffers[0].reserve(surfsize);
  buffers[1].reserve(surfsize);
  bool flipflop = true;
  hvmesh->synchronize(Mesh::FACES_E);

  while (bi != ei)
  {
    if (!visited[(unsigned int)*bi])
    {
      buffers[flipflop].clear();
      buffers[flipflop].push_back(*bi);

      typename FSRC::mesh_type::Node::array_type hvnodes;
      hvmesh->get_nodes(hvnodes, *bi);
      nodeisdiagonal[hvnodes[0]] = true;
      nodeisdiagonal[hvnodes[2]] = true;
      nodeisdiagonal[hvnodes[5]] = true;
      nodeisdiagonal[hvnodes[7]] = true;

      while (buffers[flipflop].size() > 0)
      {
	for (unsigned int i = 0; i < buffers[flipflop].size(); i++)
	{
	  if (visited[(unsigned int)buffers[flipflop][i]]) { continue; }
	  visited[(unsigned int)buffers[flipflop][i]] = true;

	  hvmesh->get_nodes(hvnodes, buffers[flipflop][i]);
	  ASSERT(hvnodes.size() == 8);

	  if (nodeisdiagonal[hvnodes[0]] || nodeisdiagonal[hvnodes[2]] ||
              nodeisdiagonal[hvnodes[5]] || nodeisdiagonal[hvnodes[7]])
	  {
            nodeisdiagonal[hvnodes[0]] = true;
            nodeisdiagonal[hvnodes[2]] = true;
            nodeisdiagonal[hvnodes[5]] = true;
            nodeisdiagonal[hvnodes[7]] = true;
	    tvmesh->add_tet((TVMesh::Node::index_type)(hvnodes[0]),
			    (TVMesh::Node::index_type)(hvnodes[1]),
			    (TVMesh::Node::index_type)(hvnodes[2]),
			    (TVMesh::Node::index_type)(hvnodes[5]));

	    tvmesh->add_tet((TVMesh::Node::index_type)(hvnodes[0]),
			    (TVMesh::Node::index_type)(hvnodes[2]),
			    (TVMesh::Node::index_type)(hvnodes[3]),
			    (TVMesh::Node::index_type)(hvnodes[7]));

	    tvmesh->add_tet((TVMesh::Node::index_type)(hvnodes[0]),
			    (TVMesh::Node::index_type)(hvnodes[5]),
			    (TVMesh::Node::index_type)(hvnodes[2]),
			    (TVMesh::Node::index_type)(hvnodes[7]));

	    tvmesh->add_tet((TVMesh::Node::index_type)(hvnodes[0]),
			    (TVMesh::Node::index_type)(hvnodes[5]),
			    (TVMesh::Node::index_type)(hvnodes[7]),
			    (TVMesh::Node::index_type)(hvnodes[4]));

	    tvmesh->add_tet((TVMesh::Node::index_type)(hvnodes[5]),
			    (TVMesh::Node::index_type)(hvnodes[2]),
			    (TVMesh::Node::index_type)(hvnodes[7]),
			    (TVMesh::Node::index_type)(hvnodes[6]));
	  }
	  else
	  {
            nodeisdiagonal[hvnodes[1]] = true;
            nodeisdiagonal[hvnodes[3]] = true;
            nodeisdiagonal[hvnodes[4]] = true;
            nodeisdiagonal[hvnodes[6]] = true;
	    tvmesh->add_tet((TVMesh::Node::index_type)(hvnodes[0]),
			    (TVMesh::Node::index_type)(hvnodes[1]),
			    (TVMesh::Node::index_type)(hvnodes[3]),
			    (TVMesh::Node::index_type)(hvnodes[4]));

	    tvmesh->add_tet((TVMesh::Node::index_type)(hvnodes[1]),
			    (TVMesh::Node::index_type)(hvnodes[2]),
			    (TVMesh::Node::index_type)(hvnodes[3]),
			    (TVMesh::Node::index_type)(hvnodes[6]));

	    tvmesh->add_tet((TVMesh::Node::index_type)(hvnodes[1]),
			    (TVMesh::Node::index_type)(hvnodes[3]),
			    (TVMesh::Node::index_type)(hvnodes[4]),
			    (TVMesh::Node::index_type)(hvnodes[6]));

	    tvmesh->add_tet((TVMesh::Node::index_type)(hvnodes[1]),
			    (TVMesh::Node::index_type)(hvnodes[5]),
			    (TVMesh::Node::index_type)(hvnodes[6]),
			    (TVMesh::Node::index_type)(hvnodes[4]));

	    tvmesh->add_tet((TVMesh::Node::index_type)(hvnodes[3]),
			    (TVMesh::Node::index_type)(hvnodes[4]),
			    (TVMesh::Node::index_type)(hvnodes[6]),
			    (TVMesh::Node::index_type)(hvnodes[7]));
	  }

	  elemmap.push_back(buffers[flipflop][i]);

	  typename FSRC::mesh_type::Cell::array_type neighbors;
	  hvmesh->get_neighbors(neighbors, buffers[flipflop][i]);

	  for (unsigned int i = 0; i < neighbors.size(); i++)
	  {
	    if (!visited[(unsigned int)neighbors[i]])
	    {
	      buffers[!flipflop].push_back(neighbors[i]);
	    }
	  }
	}
	buffers[flipflop].clear();
	flipflop = !flipflop;
      }
    }
    ++bi;
  }
  typedef vector<typename FSRC::value_type>          FDat; 

  if (hvfield->basis_order() == -1) {
    typedef NoDataBasis<typename FSRC::value_type>     DatBasis;
    typedef GenericField<TVMesh, DatBasis, FDat>       TVField;   
    TVField *tvfield =  scinew TVField(tvmesh);
    dstH = tvfield;
  } else if (hvfield->basis_order() == 0) {
    typedef ConstantBasis<typename FSRC::value_type>    DatBasis;
    typedef GenericField<TVMesh, DatBasis, FDat>       TVField;   
    TVField *tvfield =  scinew TVField(tvmesh);
    dstH = tvfield;
    typename FSRC::value_type val;
    for (unsigned int i = 0; i < elemmap.size(); i++)
    {
      hvfield->value(val, elemmap[i]);
      tvfield->set_value(val, (TVMesh::Elem::index_type)(i*5+0));
      tvfield->set_value(val, (TVMesh::Elem::index_type)(i*5+1));
      tvfield->set_value(val, (TVMesh::Elem::index_type)(i*5+2));
      tvfield->set_value(val, (TVMesh::Elem::index_type)(i*5+3));
      tvfield->set_value(val, (TVMesh::Elem::index_type)(i*5+4));
    }

  } else if (hvfield->basis_order() == 1) {
    typedef TetLinearLgn<typename FSRC::value_type>    DatBasis;
    typedef GenericField<TVMesh, DatBasis, FDat>       TVField;   
    TVField *tvfield =  scinew TVField(tvmesh);
    dstH = tvfield;
    typename FSRC::value_type val;
    for (unsigned int i = 0; i < hnsize; i++)
    {
      hvfield->value(val, (typename FSRC::mesh_type::Node::index_type)(i));
      tvfield->set_value(val, (TVMesh::Node::index_type)(i));
    }
  }
  else
  {
    reporter->warning("Could not load data values onto output field.");
    reporter->warning("Use DirectInterp if data values are required.");
  }

  dstH->copy_properties(hvfield);
  return true;
}


class LatToTetAlgo : public DynamicAlgoBase
{
public:
  virtual bool execute(ProgressReporter *reporter,
                       FieldHandle src, FieldHandle &dst) = 0;

  //! support the dynamically compiled algorithm concept
  static CompileInfoHandle get_compile_info(const TypeDescription *data_td);
};


template <class FSRC>
class LatToTetAlgoT : public LatToTetAlgo
{
public:
  //! virtual interface. 
  virtual bool execute(ProgressReporter *reporter,
                       FieldHandle src, FieldHandle& dst);
};


template <class FSRC>
bool
LatToTetAlgoT<FSRC>::execute(ProgressReporter *reporter,
                             FieldHandle srcH, FieldHandle& dstH) 
{
  FSRC *hvfield = dynamic_cast<FSRC*>(srcH.get_rep());

  typename FSRC::mesh_type *hvmesh = hvfield->get_typed_mesh().get_rep();
  TVMesh::handle_type tvmesh = scinew TVMesh();

  typename FSRC::mesh_type::Node::size_type lnsize;
  hvmesh->size(lnsize);
  tvmesh->node_reserve((unsigned int)lnsize);

  // Copy points directly, assuming they will have the same order.
  typename FSRC::mesh_type::Node::iterator nbi, nei;
  hvmesh->begin(nbi); hvmesh->end(nei);
  while (nbi != nei)
  {
    Point p;
    hvmesh->get_center(p, *nbi);
    tvmesh->add_point(p);
    ++nbi;
  }

  typename FSRC::mesh_type::Elem::size_type lesize;
  hvmesh->size(lesize);
  tvmesh->elem_reserve((unsigned int)lesize * 5);

  typename FSRC::mesh_type::Elem::iterator bi, ei;
  hvmesh->begin(bi); hvmesh->end(ei);
  while (bi != ei)
  {
    typename FSRC::mesh_type::Node::array_type hvnodes(8);
    hvmesh->get_nodes(hvnodes, *bi);
    if (!(((*bi).i_ ^ (*bi).j_ ^ (*bi).k_)&1))
    {
      tvmesh->add_tet((TVMesh::Node::index_type)((unsigned int)hvnodes[0]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[1]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[2]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[5]));

      tvmesh->add_tet((TVMesh::Node::index_type)((unsigned int)hvnodes[0]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[2]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[3]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[7]));

      tvmesh->add_tet((TVMesh::Node::index_type)((unsigned int)hvnodes[0]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[5]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[2]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[7]));

      tvmesh->add_tet((TVMesh::Node::index_type)((unsigned int)hvnodes[0]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[5]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[7]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[4]));

      tvmesh->add_tet((TVMesh::Node::index_type)((unsigned int)hvnodes[5]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[2]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[7]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[6]));
    }
    else
    {
      tvmesh->add_tet((TVMesh::Node::index_type)((unsigned int)hvnodes[0]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[1]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[3]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[4]));

      tvmesh->add_tet((TVMesh::Node::index_type)((unsigned int)hvnodes[1]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[2]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[3]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[6]));

      tvmesh->add_tet((TVMesh::Node::index_type)((unsigned int)hvnodes[1]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[3]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[4]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[6]));

      tvmesh->add_tet((TVMesh::Node::index_type)((unsigned int)hvnodes[1]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[5]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[6]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[4]));

      tvmesh->add_tet((TVMesh::Node::index_type)((unsigned int)hvnodes[3]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[4]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[6]),
		      (TVMesh::Node::index_type)((unsigned int)hvnodes[7]));
    }
    ++bi;
  }
  

  typedef vector<typename FSRC::value_type>          FDat; 

  if (hvfield->basis_order() == -1) {
    typedef NoDataBasis<typename FSRC::value_type>     DatBasis;
    typedef GenericField<TVMesh, DatBasis, FDat>       TVField;   
    TVField *tvfield =  scinew TVField(tvmesh);
    dstH = tvfield;
  } else if (hvfield->basis_order() == 0) {
    typedef ConstantBasis<typename FSRC::value_type>    DatBasis;
    typedef GenericField<TVMesh, DatBasis, FDat>       TVField;   
    TVField *tvfield =  scinew TVField(tvmesh);
    dstH = tvfield;
    typename FSRC::value_type val;
    typename FSRC::mesh_type::Cell::iterator cbi, cei;    
    hvmesh->begin(cbi); hvmesh->end(cei);
    while (cbi != cei)
    {
      hvfield->value(val, *cbi);
      unsigned int i = (unsigned int)*cbi;
      tvfield->set_value(val, (TVMesh::Cell::index_type)(i*5+0));
      tvfield->set_value(val, (TVMesh::Cell::index_type)(i*5+1));
      tvfield->set_value(val, (TVMesh::Cell::index_type)(i*5+2));
      tvfield->set_value(val, (TVMesh::Cell::index_type)(i*5+3));
      tvfield->set_value(val, (TVMesh::Cell::index_type)(i*5+4));
      ++cbi;
    }
  } else if (hvfield->basis_order() == 1) {
    typedef TetLinearLgn<typename FSRC::value_type>    DatBasis;
    typedef GenericField<TVMesh, DatBasis, FDat>       TVField;   
    TVField *tvfield =  scinew TVField(tvmesh);
    dstH = tvfield;
    typename FSRC::value_type val;
    hvmesh->begin(nbi); hvmesh->end(nei);
    while (nbi != nei)
    {
      hvfield->value(val, *nbi);
      tvfield->set_value(val, (TVMesh::Node::index_type)(unsigned int)(*nbi));
      ++nbi;
    }
  }
  else
  {
    reporter->warning("Could not load data values onto output field.");
    reporter->warning("Use DirectInterp if data values are required.");
  }
  
  dstH->copy_properties(hvfield);
  return true;
}

}

#endif
