
/*
 *  Array3.h: Interface to dynamic 3D array class
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#ifndef SCI_Containers_Array3_h
#define SCI_Containers_Array3_h 1

#include <SCICore/Util/Assert.h>
#include <SCICore/Persistent/Persistent.h>

namespace SCICore {

namespace Tester {
  class RigorousTest;
}

namespace Geometry {
  void Pio();  // This is a dummy declaration to get things to compile.
}

namespace Containers {

using SCICore::PersistentSpace::Piostream;
using SCICore::Tester::RigorousTest;

/**************************************

CLASS
   Array3
   
KEYWORDS
   Array3

DESCRIPTION
    Array3.h: Interface to dynamic 3D array class
  
    Written by:
     Steven G. Parker
     Department of Computer Science
     University of Utah
     March 1994
  
    Copyright (C) 1994 SCI Group
PATTERNS
   
WARNING
  
****************************************/

template<class T>
class Array3 {
    T*** objs;
    int dm1;
    int dm2;
    int dm3;
    void allocate();
public:
    //////////
    //Default Constructor
    Array3();
    
    //////////
    //Copy Constructor
    Array3(const Array3&);

    //////////
    //Constructor
    Array3(int, int, int);
    
    //////////
    //Assignment Operator
    Array3<T>& operator=(const Array3&);
    
    //////////
    //Class Destructor
    ~Array3();
    
    //////////
    //Access the nXnXn element of the array
    inline T& operator()(int d1, int d2, int d3) const
	{
	    ASSERTL3(d1>=0 && d1<dm1);
	    ASSERTL3(d2>=0 && d2<dm2);
	    ASSERTL3(d3>=0 && d3<dm3);
	    return objs[d1][d2][d3];
	}
    
    //////////
    //Returns the number of spaces in dim1	    
    inline int dim1() const {return dm1;}
    //////////
    //Returns the number of spaces in dim2
    inline int dim2() const {return dm2;}
    //////////
    //Returns the number of spaces in dim3
    inline int dim3() const {return dm3;}
    
    //////////
    //Re-size the Array
    void newsize(int, int, int);

    //////////
    //Initialize all elements to T
    void initialize(const T&);

    T* get_onedim();
    void get_onedim_byte( unsigned char *v );

    inline T*** get_dataptr() {return objs;}

    //////////
    //Rigorous Tests
    static void test_rigorous(RigorousTest* __test);
    
    friend void TEMPLATE_TAG Pio TEMPLATE_BOX (Piostream&, Array3<T>&);
    friend void TEMPLATE_TAG Pio TEMPLATE_BOX (Piostream&, Array3<T>*&);

};

} // End namespace Containers
} // End namespace SCICore

////////////////////////////////////////////////////////////
//
// Start of included Array3.cc
//

#include <SCICore/Containers/String.h>
#include <SCICore/Malloc/Allocator.h>

namespace SCICore {
namespace Containers {

template<class T>
Array3<T>::Array3()
{
    objs=0;
}

template<class T>
void Array3<T>::allocate()
{
    objs=scinew T**[dm1];
    T** p=scinew T*[dm1*dm2];
    T* pp=scinew T[dm1*dm2*dm3];
    for(int i=0;i<dm1;i++){
	objs[i]=p;
	p+=dm2;
	for(int j=0;j<dm2;j++){
	    objs[i][j]=pp;
	    pp+=dm3;
	}
    }
}

template<class T>
void Array3<T>::newsize(int d1, int d2, int d3)
{
    if(objs && dm1==d2 && dm2==d2 && dm3==d3)return;
    dm1=d1;
    dm2=d2;
    dm3=d3;
    if(objs){
	delete[] objs[0][0];
	delete[] objs[0];
	delete[] objs;
    }
    allocate();
}

template<class T>
Array3<T>::Array3(const Array3<T>& a)
: dm1(a.dm1), dm2(a.dm2), dm3(a.dm3)
{
    allocate();
}

template<class T>
Array3<T>::Array3(int dm1, int dm2, int dm3)
: dm1(dm1), dm2(dm2),dm3(dm3)
{
    allocate();
}

template<class T>
Array3<T>::~Array3()
{
    if(objs){
	delete[] objs[0][0];
	delete[] objs[0];
	delete[] objs;
    }
}

template<class T>
void Array3<T>::initialize(const T& t)
{
    ASSERT(objs != 0);
    for(int i=0;i<dm1;i++){
	for(int j=0;j<dm2;j++){
	    for(int k=0;k<dm3;k++){
		objs[i][j][k]=t;
	    }
	}
    }
}

template<class T>
T* Array3<T>::get_onedim()
{
  int i,j,k, index;
  T* a = scinew T[dm1*dm2*dm3];
  
  index=0;
  for( i=0; i<dm1; i++)
    for( j=0; j<dm2; j++ )
      for( k=0; k<dm3; k++ )
	a[index++] = objs[i][j][k];
  return a;
}

template<class T>
void
Array3<T>::get_onedim_byte( unsigned char *v )
{
  int i,j,k, index;
  index = 0;
  
  for( k=0; k<dm3; k++ )
    for( j=0; j<dm2; j++ )
      for( i=0; i<dm1; i++)
	v[index++] = objs[i][j][k];
}

#define ARRAY3_VERSION 1

template<class T>
void Pio(Piostream& stream, Containers::Array3<T>& data)
{
    using SCICore::PersistentSpace::Pio;
    using SCICore::Geometry::Pio;

    /*int version=*/stream.begin_class("Array3", ARRAY3_VERSION);
    if(stream.reading()){
	// Allocate the array...
	int d1, d2, d3;
	Pio(stream, d1);
	Pio(stream, d2);
	Pio(stream, d3);
	data.newsize(d1, d2, d3);
    } else {
	Pio(stream, data.dm1);
	Pio(stream, data.dm2);
	Pio(stream, data.dm3);
    }
    for(int i=0;i<data.dm1;i++){
	for(int j=0;j<data.dm2;j++){
	    for(int k=0;k<data.dm3;k++){
		Pio(stream, data.objs[i][j][k]);
	    }
	}
    }
    stream.end_class();
}

template<class T>
void Pio(Piostream& stream, Containers::Array3<T>*& data) {
    if (stream.reading()) {
	data=scinew Array3<T>;
    }
    Containers::Pio(stream, *data);
}

} // End namespace Containers
} // End namespace SCICore

//
// $Log$
// Revision 1.3  1999/08/19 23:18:04  sparker
// Removed a bunch of #include <SCICore/Util/NotFinished.h> statements
// from files that did not need them.
//
// Revision 1.2  1999/08/17 06:38:35  sparker
// Merged in modifications from PSECore to make this the new "blessed"
// version of SCIRun/Uintah.
//
// Revision 1.1  1999/07/27 16:56:11  mcq
// Initial commit
//
// Revision 1.4  1999/07/07 21:10:35  dav
// added beginnings of support for g++ compilation
//
// Revision 1.3  1999/05/06 19:55:42  dav
// added back .h files
//
// Revision 1.1  1999/05/05 21:04:29  dav
// added SCICore .h files to /include directories
//
// Revision 1.1.1.1  1999/04/24 23:12:26  dav
// Import sources
//
//

#endif


