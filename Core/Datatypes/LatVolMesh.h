/*
 *  LatVolMesh.h: Templated Mesh defined on a 3D Regular Grid
 *
 *  Written by:
 *   Michael Callahan
 *   Department of Computer Science
 *   University of Utah
 *   January 2001
 *
 *  Copyright (C) 2001 SCI Group
 *
 */

#ifndef SCI_project_LatVolMesh_h
#define SCI_project_LatVolMesh_h 1

#include <Core/Geometry/Point.h>
#include <Core/Containers/LockingHandle.h>
#include <Core/Datatypes/MeshBase.h>
#include <string>

namespace SCIRun {

using std::string;

class SCICORESHARE LatVolMesh : public MeshBase
{
public:
  struct LatIndex
  {
  public:
    LatIndex() : i_(0), j_(0), k_(0) {}
    LatIndex(unsigned i, unsigned j, unsigned k) : i_(i), j_(j), k_(k) {}
    
    unsigned i_, j_, k_;
  };

  struct CellIndex : public LatIndex
  {
    CellIndex() : LatIndex() {}
    CellIndex(unsigned i, unsigned j, unsigned k) : LatIndex(i,j,k) {}    
  };
  
  struct NodeIndex : public LatIndex
  {
    NodeIndex() : LatIndex() {}
    NodeIndex(unsigned i, unsigned j, unsigned k) : LatIndex(i,j,k) {}    
  };
  
  struct LatIter : public LatIndex
  {
    LatIter(const LatVolMesh *m, unsigned i, unsigned j, unsigned k) 
      : LatIndex(i, j, k), mesh_(m) {}
    
    const LatIndex &operator *() { return *this; }
    
    bool operator ==(const LatIter &a) const
    {
      return i_ == a.i_ && j_ == a.j_ && k_ == a.k_ && mesh_ == a.mesh_;
    }
    
    bool operator !=(const LatIter &a) const
    {
      return !(*this == a);
    }
    
    const LatVolMesh *mesh_;
  };
  
  
  struct NodeIter : public LatIter
  {
    NodeIter(const LatVolMesh *m, unsigned i, unsigned j, unsigned k) 
      : LatIter(m, i, j, k) {}
    
    const NodeIndex &operator *() const { return (const NodeIndex&)(*this); }

    NodeIter &operator++()
    {
      i_++;
      if (i_ >= mesh_->nx_)	{
	i_ = 0;
	j_++;
	if (j_ >= mesh_->ny_) {
	  j_ = 0;
	  k_++;
	}
      }
      return *this;
    }
    
  private:

    NodeIter operator++(int)
    {
      NodeIter result(*this);
      operator++();
      return result;
    }
  };
  
  
  struct CellIter : public LatIter
  {
    CellIter(const LatVolMesh *m, unsigned i, unsigned j, unsigned k) 
      : LatIter(m, i, j, k) {}
    
    const CellIndex &operator *() const { return (const CellIndex&)(*this); }

    CellIter &operator++()
    {
      i_++;
      if (i_ >= mesh_->nx_-1) {
	i_ = 0;
	j_++;
	if (j_ >= mesh_->ny_-1) {
	  j_ = 0;
	  k_++;
	}
      }
      return *this;
    }
    
  private:

    CellIter operator++(int)
    {
      CellIter result(*this);
      operator++();
      return result;
    }
  };
  
  typedef LatIndex index_type;
  
  //! Index and Iterator types required for Mesh Concept.
  typedef NodeIndex       node_index;
  typedef NodeIter        node_iterator;
 
  //typedef EdgeIndex       edge_index;     
  //typedef EdgeIterator    edge_iterator;
  typedef int             edge_index;
  typedef int             edge_iterator;
 
  //typedef FaceIndex       face_index;
  //typedef FaceIterator    face_iterator;
  typedef unsigned        face_index;
  typedef unsigned        face_iterator;
 
  typedef CellIndex       cell_index;
  typedef CellIter        cell_iterator;

  // storage types for get_* functions
  typedef vector<node_index>  node_array;
  typedef vector<edge_index>  edge_array;
  typedef vector<face_index>  face_array;
  typedef vector<cell_index>  cell_array;
  typedef vector<double>      weight_array;

  friend class NodeIter;
  friend class CellIter;

  LatVolMesh()
    : nx_(1),ny_(1),nz_(1),min_(Point(0,0,0)),max_(Point(1,1,1)) {}
  LatVolMesh(unsigned x, unsigned y, unsigned z, Point &min, Point &max) 
    : nx_(x),ny_(y),nz_(z),min_(min),max_(max) {}
  LatVolMesh(const LatVolMesh &copy)
    : nx_(copy.get_nx()),ny_(copy.get_ny()),nz_(copy.get_nz()),
      min_(copy.get_min()),max_(copy.get_max()) {}
  virtual ~LatVolMesh() {}

  node_iterator node_begin() const { return node_iterator(this, 0, 0, 0); }
  node_iterator node_end() const { return node_iterator(this, 0, 0, nz_); }
  index_type    nodes_size() const { return LatIndex(nx_, ny_, nz_); }

  edge_iterator edge_begin() const { return 0; }
  edge_iterator edge_end() const { return 0; }
  int edges_size() const { return nx_*ny_*nz_*3 -
			     (nx_*ny_ + ny_*nz_ + nx_*nz_); }

  cell_iterator cell_begin() const { return cell_iterator(this, 0, 0, 0); }
  cell_iterator cell_end() const { return cell_iterator(this, 0, 0, nz_-1); }
  index_type    cells_size() const { return LatIndex(nx_-1, ny_-1, nz_-1); }

  face_iterator face_begin() const { return 0; }
  face_iterator face_end() const { return 0; }
  int           faces_size() const { return (nx_-1)*(ny_-1)*(nz_-1)*3 + 
				       (nx_-1)*(ny_-1) +
				       (ny_-1)*(nz_-1) + 
				       (nx_-1)*(nz_-1); }


  node_index nodes_size() { return node_index(nx_, ny_, nz_); }
  edge_index edges_size()
  { return nx_*ny_*nz_*3 - nx_*ny_ - ny_*nz_ -  nx_*nz_; } 
  face_index faces_size()
  {
    return (nx_-1)*(ny_-1)*(nz_-1)*3 +
      (nx_-1)*(ny_-1) + (ny_-1)*(nz_-1) + (nx_-1)*(nz_-1);
  }
  cell_index cells_size() { return cell_index(nx_-1, ny_-1, nz_-1); }

  //! get the mesh statistics
  unsigned get_nx() const { return nx_; }
  unsigned get_ny() const { return ny_; }
  unsigned get_nz() const { return nz_; }
  Point get_min() const { return min_; }
  Point get_max() const { return max_; }
  virtual BBox get_bounding_box() const;

  //! set the mesh statistics
  void set_nx(unsigned x) { nx_ = x; }
  void set_ny(unsigned y) { ny_ = y; }
  void set_nz(unsigned z) { nz_ = z; }
  void set_min(Point p) { min_ = p; }
  void set_max(Point p) { max_ = p; }

  //! get the child elements of the given index
  void get_nodes(node_array &, edge_index) const {}
  void get_nodes(node_array &, face_index) const {}
  void get_nodes(node_array &, cell_index) const;
  void get_edges(edge_array &, face_index) const {}
  void get_edges(edge_array &, cell_index) const {}
  void get_faces(face_array &, cell_index) const {}

  //! get the parent element(s) of the given index
  unsigned get_edges(edge_array &, node_index) const { return 0; }
  unsigned get_faces(face_array &, node_index) const { return 0; }
  unsigned get_faces(face_array &, edge_index) const { return 0; }
  unsigned get_cells(cell_array &, node_index) const { return 0; }
  unsigned get_cells(cell_array &, edge_index) const { return 0; }
  unsigned get_cells(cell_array &, face_index) const { return 0; }

  //! similar to get_cells() with face_index argument, but
  //  returns the "other" cell if it exists, not all that exist
  void get_neighbor(cell_index &, face_index) const {}

  //! get the center point (in object space) of an element
  void get_center(Point &, node_index) const;
  void get_center(Point &, edge_index) const {}
  void get_center(Point &, face_index) const {}
  void get_center(Point &, cell_index) const;

  bool locate(node_index &, const Point &) const;
  bool locate(edge_index &, const Point &) const { return false; }
  bool locate(face_index &, const Point &) const { return false; }
  bool locate(cell_index &, const Point &) const;

  void unlocate(Point &result, const Point &p) const { result =  p; };

  void get_point(Point &, node_index) const;

  virtual void io(Piostream&);
  static PersistentTypeID type_id;
  static  const string type_name(int n = -1);
  virtual const string get_type_name(int n = -1) const { return type_name(n); }

private:

  //! the node_index space extents of a LatVolMesh (min=0, max=n-1)
  unsigned nx_, ny_, nz_;

  //! the object space extents of a LatVolMesh
  Point min_, max_;

  // returns a LatVolMesh
  static Persistent *maker() { return new LatVolMesh(); }
};

typedef LockingHandle<LatVolMesh> LatVolMeshHandle;



} // namespace SCIRun

#endif // SCI_project_LatVolMesh_h





