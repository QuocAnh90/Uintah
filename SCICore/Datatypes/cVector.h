/*
 *  cVector.h : ?
 *
 *  Written by:
 *   Author: ?
 *   Department of Computer Science
 *   University of Utah
 *   Date: ?
 *
 *  Copyright (C) 199? SCI Group
 */

#ifndef cVECTOR_H
#define cVECTOR_H 1

#include <SCICore/share/share.h>

#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <SCICore/Math/Complex.h>

#include <SCICore/CoreDatatypes/Datatype.h>
#include <SCICore/Containers/LockingHandle.h>

namespace SCICore {
namespace CoreDatatypes {

using SCICore::Math::Complex;
using SCICore::Containers::LockingHandle;
using SCICore::PersistentSpace::Piostream;
using SCICore::PersistentSpace::PersistentTypeID;

class cVector;
typedef LockingHandle<cVector> cVectorHandle;

class SCICORESHARE cVector :public Datatype{

  friend class cDMatrix;
  friend class cSMatrix;
  
private:
  Complex *a;
  int Size;
  
public:

// Persistent representation...
    virtual void io(Piostream&);
    static PersistentTypeID type_id;

  
  cVector (int N); //constructor
  cVector (const cVector &B); //copy constructor;
  cVector &operator=(const cVector &B); //assigment operator  
  ~cVector(); //Destructor;
  
  int size() {return Size;};
  double norm();
  Complex &operator()(int i);
  void conj();
  int load(char* filename);
  
  
  void set(const cVector& B);
  void add(const cVector& B);
  void subtr(const cVector& B);
  void mult(const Complex x); 
  
  cVector operator+(const cVector& B) const;
  cVector operator-(const cVector& B) const;
  
  friend SCICORESHARE Complex operator*(cVector& A, cVector& B);
  friend SCICORESHARE cVector  operator*(const cVector& B,Complex x);
  friend SCICORESHARE cVector  operator*(Complex x, const cVector &B);
  friend SCICORESHARE cVector  operator*(const cVector& B,double x);
  friend SCICORESHARE cVector  operator*(double x, const cVector &B);
  
  friend SCICORESHARE ostream &operator<< (ostream &output, cVector &B);
  
};

} // End namespace CoreDatatypes
} // End namespace SCICore

//
// $Log$
// Revision 1.2  1999/08/17 06:39:02  sparker
// Merged in modifications from PSECore to make this the new "blessed"
// version of SCIRun/Uintah.
//
// Revision 1.1  1999/07/27 16:56:35  mcq
// Initial commit
//
// Revision 1.3  1999/05/06 19:56:01  dav
// added back .h files
//
// Revision 1.1  1999/05/05 21:04:53  dav
// added SCICore .h files to /include directories
//
// Revision 1.1  1999/04/25 04:07:25  dav
// Moved files into CoreDatatypes
//
// Revision 1.1.1.1  1999/04/24 23:12:48  dav
// Import sources
//
//

#endif



