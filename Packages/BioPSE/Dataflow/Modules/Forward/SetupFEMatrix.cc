/*
 *  SetupFEMatrix.cc:  Setups the global finite element matrix
 *
 *  Written by:
 *   Ruth Nicholson Klepfer
 *   Department of Bioengineering
 *   University of Utah
 *   October 1994
 *
 *  Modified:
 *   Alexei Samsonov
 *   Department of Computer Science
 *   University of Utah
 *   March 2001    
 *  Copyright (C) 1994, 2001 SCI Group
 */

#include <Dataflow/Network/Module.h>
#include <Dataflow/Ports/MatrixPort.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Dataflow/Ports/FieldPort.h>
#include <Core/Datatypes/TetVol.h>
#include <Packages/BioPSE/Core/Algorithms/NumApproximation/BuildFEMatrix.h>
#include <Core/Geometry/Point.h>
#include <Core/Malloc/Allocator.h>
#include <Core/GuiInterface/GuiVar.h>
#include <Core/Thread/Barrier.h>
#include <Core/Thread/Parallel.h>
#include <Core/Thread/Thread.h>
#include <iostream>

using std::cerr;
using std::endl;

#define PINVAL 0

namespace BioPSE {

using namespace SCIRun;
typedef LockingHandle<TetVol<int> >    CondMeshHandle;
typedef LockingHandle<TetVol<double> > DirBCMeshHandle;
typedef LockingHandle<FieldSet>        FieldSetHanlde;

class SetupFEMatrix : public Module {
  
  //! Private data
  FieldIPort*        iportField_;
  MatrixOPort*       oportRhs_;
  MatrixOPort*       oportMtrx_;
 
  GuiInt             uiUseCond_;
  int                lastUseCond_;
  
  GuiString          uiBCFlag_;
  
  MatrixHandle       hGblMtrx_;
  MatrixHandle       hRhs_;
  int                gen_;
  int                refNode_;

public:
  
  //! Constructor/Destructor
  SetupFEMatrix(const clString& id);
  virtual ~SetupFEMatrix();

  //! Public methods
  virtual void execute();
};

extern "C" Module* make_SetupFEMatrix(const clString& id) {
  return new SetupFEMatrix(id);
}

SetupFEMatrix::SetupFEMatrix(const clString& id): 
  Module("SetupFEMatrix", id, Filter), 
  uiUseCond_("UseCondTCL", id, this),
  lastUseCond_(1),
  uiBCFlag_("BCFlag", id, this),
  refNode_(0)
{
  // Create the input ports
  iportField_ = scinew FieldIPort(this, "Mesh", FieldIPort::Atomic);
  add_iport(iportField_);
  
  // Create the output ports
  oportMtrx_=scinew MatrixOPort(this, "FEM Matrix", MatrixIPort::Atomic);
  add_oport(oportMtrx_);
  
  oportRhs_=scinew MatrixOPort(this, "RHS", MatrixIPort::Atomic);
  add_oport(oportRhs_);
  gen_=-1;
  uiUseCond_.set(1);
}

SetupFEMatrix::~SetupFEMatrix(){
}

void SetupFEMatrix::execute(){
 
  FieldHandle hField;
  if(!iportField_->get(hField)){
    msgStream_ << "Error: Cann't get field" << endl;
    return;
  }
  
  if (hField->generation == gen_ 
      && hGblMtrx_.get_rep() 
      && hRhs_.get_rep()
      && lastUseCond_==uiUseCond_.get()) {
    oportMtrx_->send(hGblMtrx_);
    oportRhs_->send(hRhs_);
    return;
  }
  
  gen_ = hField->generation;
  
  CondMeshHandle hCondMesh;
  if (hField->get_type_name(0)=="TetVol" && hField->get_type_name(1)=="int"){
    
    hCondMesh = dynamic_cast<TetVol<int>* >(hField.get_rep());
    
    if (!hCondMesh.get_rep()){
      msgStream_ << "Error: Unable to cast to TetVol<int>*" << endl;
      return;
    }
  }
  else {
    msgStream_ << "Error: The mesh in conductivity tensor field is not of type TetVol<int>" << endl;
    return;
  }
  
  //-- polling Field for Dirichlet BC
  vector<pair<int, double> > dirBC;
  if (!hField->get("dirichlet", dirBC)){
    msgStream_ << "The Field Set doesn't contain Dirichlet boundary conditions" << endl;
    //! Pinging 0 node to 0
    dirBC.push_back(pair<int, double>(0, 0.0));
  }
  else {
    clString bcFlag = uiBCFlag_.get();
    if (bcFlag=="PinZero"){
      msgStream_ << "Pinging to 0" << endl;
      dirBC.erase(dirBC.begin(), dirBC.end());
      dirBC.push_back(pair<int, double>(0, 0.0));
    }
  }

  //! finding conductivity tensor lookup table
  Array1<Tensor> tens;

  if (uiUseCond_.get()==1 && hCondMesh->get("conductivity_table", tens)){
    msgStream_ << "Using supplied conductivity tensors "  << endl;
  }
  else {
    msgStream_ << "Using identity conductivity tensors" << endl;
    tens.resize(256);
    Array1<double> t(6);
    t[0] = t[3] = t[5] = 1;
    t[1] = t[2] = t[5] = 0;
    tens.initialize(Tensor(t));
  }
  
  lastUseCond_ = uiUseCond_.get();
  
  if(BuildFEMatrix::build_FEMatrix(hCondMesh, dirBC, tens, hGblMtrx_, hRhs_)){
    msgStream_ << "Matrix is ready" << endl;
    msgStream_ << "Size: " << hGblMtrx_->nrows() << "-by-" << hGblMtrx_->ncols() << endl;
  };
  
  //! outputing
  oportMtrx_->send(hGblMtrx_);
  oportRhs_->send(hRhs_);
}

} // End namespace BioPSE
