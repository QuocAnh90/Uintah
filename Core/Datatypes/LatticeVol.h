
#ifndef Datatypes_LatticeVol_h
#define Datatypes_LatticeVol_h

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/MeshRG.h>
#include <Core/Containers/LockingHandle.h>
#include <vector>

namespace SCIRun {

template <class Data>
class LatticeVol: public Field 
{
public:
  //! Typedefs to support the Field concept.
  typedef Data            value_type;
  typedef MeshRG          mesh_type;
  typedef vector<Data>    fdata_type;


  LatticeVol();
  LatticeVol(data_location data_at);
  virtual ~LatticeVol();

  //! Required virtual functions from field base.
  virtual MeshBaseHandle get_mesh() { return (MeshBaseHandle)mesh_; } const;

  //! Required interfaces from field base.
  virtual InterpolateToScalar* query_interpolate_to_scalar() const;



  //! Required interface to support Field Concept.
  value_type operator[] (int);
  
  /* this should go in the field base class
  template <class Functor>
  void interpolate(const Point &p, Functor &f);
  */
  
  MeshRGHandle get_rg_mesh() { return mesh_; }

  //! attempt to set data, succeeds if sizes match.
  bool set_fdata(fdata_type *fdata) { fdata_ = fdata; }



  //! Persistent I/O.
  void    io(Piostream &stream);
  static  PersistentTypeID type_id;
  static  const string type_name(int);
  virtual const string get_type_name(int n) const;
private:
  //! A Lattice Mesh.
  MeshRGHandle                mesh_;
  //! Data container.
  fdata_type                  fdata_;
  
};

const double LATTICE_VOL_VERSION = 1.0;

template <class Data>
void LatticeVol<Data>::io(Piostream& stream){

  stream.begin_class(typeName().c_str(), LATTICE_VOL_VERSION);
  Field::io(stream);
  Pio(stream, mesh_.get_rep());
  Pio(stream, fdata_);
  stream.end_class();
}

template <class Data> template <class Functor>
void 
LatticeVol<Data>::interpolate(const Point &p, Functor &f) {

  if (f.wieghts_ == 0)
    f.wieghts_ = new double[8]; // eight nodes for hex
  MeshRG::cell_index ci;
  mesh_->locate_cell(ci, p, f.weights_);

  switch (data_at()) {
  case Field::NODE :
    {
      int i = 0;
      MeshRG::node_array nodes;
      get_nodes_from_cell(nodes, ci);
      MeshRG::node_array::iterator iter = nodes.begin();
      while (iter != nodes.end()) {
	f(*data_, *iter);
	++iter; ++i;
      }
    }
  break;
  case Field::EDGE:
    {
    }
    break;
  case Field::FACE:
    {
    }
    break;
  case Field::CELL:
    {
    }
    break;
  } 
} 

/* this should go into the field base class 

template <class Data>
struct InterpFunctor {
  typedef Data data_type;
  typedef typename Data::value_type value_type;

  InterpFunctor(int num_weights = 0) :
    result_(0),
    weights_(0) 
  {
    if (num_weights > 0) {
      weights_ = new double[num_weights];
    }
  }

  virtual ~InterpFunctor() {
    if (weights_) { delete[] weights; }
  }

  double         *weights_;
  value_type      result_;
};

// sample interp functor.
template <class Data, class Index>
struct LinearInterp : public InterpFunctor<Data> {
  
  LinearInterp(int num_weights) :
    InterpFunctor<Data>(num_weights) {}

  void 
  operator()(const Data &data, Index idx, int widx) {
      result_ += data[idx] * weights_[widx];
      cout << "linear interping :)" << endl;
    }
};

*/  // go into base field class

} // end namespace SCIRun

#endif // Datatypes_TetVol_h
















