/*
 *  StreamLines.cc:
 *
 *  Written by:
 *   moulding
 *   TODAY'S DATE HERE
 *
 */

#include <Dataflow/Network/Module.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Datatypes/Field.h>
#include <Dataflow/Ports/FieldPort.h>
#include <Core/Containers/Array1.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Datatypes/LatticeVol.h>

#include <iostream>
#include <vector>
#include <math.h>

#include <Packages/Moulding/share/share.h>

namespace Moulding {

using namespace SCIRun;

// LUTs for the RK-fehlberg algorithm 
double a[]   ={16.0/135, 0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55};
double ab[]  ={1.0/360, 0, -128.0/4275, -2197.0/75240, 1.0/50, 2.0/55};
double c[]   ={0, 1.0/4, 3.0/8, 12.0/13, 1.0, 1.0/2}; /* not used */
double d[][5]={{0, 0, 0, 0, 0},
	       {1.0/4, 0, 0, 0, 0},
	       {3.0/32, 9.0/32, 0, 0, 0},
	       {1932.0/2197, -7200.0/2197, 7296.0/2197, 0, 0},
	       {439.0/216, -8.0, 3680.0/513, -845.0/4104, 0},
	       {-8.0/27, 2.0, -3544.0/2565, 1859.0/4104, -11.0/40}};

class MouldingSHARE StreamLines : public Module {
public:
  StreamLines(const clString& id);

  virtual ~StreamLines();

  virtual void execute();

  virtual void tcl_command(TCLArgs&, void*);

private:
  // data members

  FieldIPort*       d_vfport;
  FieldIPort*       d_sfport;

  FieldHandle       d_vfhandle;
  FieldHandle       d_sfhandle;

  // member functions

  // Using Runge-Kutta-Fehlberg, find the points that should
  // be connected together to create a single streamline.
  void FindStreamLineSupport(Array1<Point>&, Point, float, float, int);

  // compute the inner terms of the RKF formula
  int ComputeRKFTerms(Array1<Vector>&,Point&, float);
};

extern "C" MouldingSHARE Module* make_StreamLines(const clString& id) {
  return scinew StreamLines(id);
}

StreamLines::StreamLines(const clString& id)
  : Module("StreamLines", id, Source)
  //, "Visualization", "Moulding")
{
  d_vfport = scinew FieldIPort(this, "Flow field", FieldIPort::Atomic);
  add_iport(d_vfport);

  d_sfport = scinew FieldIPort(this, "Seed field", FieldIPort::Atomic);
  add_iport(d_sfport);
}

StreamLines::~StreamLines()
{
}

int
StreamLines::ComputeRKFTerms(Array1<Vector>& v /* storage for terms */,
			     Point& p          /* previous point */,
			     float s           /* current step size */)
{
  int check = 0;
#if 0
  check |= vlinterpolate->vlinterpolate(p,
					v[0]);
  check |= vlinterpolate->vlinterpolate(p+
					v[0]*d[1][0],
					v[1]);
  check |= vlinterpolate->vlinterpolate(p+
					v[0]*d[2][0]+
					v[1]*d[2][1],
					v[2]);
  check |= vlinterpolate->vlinterpolate(p+
					v[0]*d[3][0]+
					v[1]*d[3][1]+
					v[2]*d[3][2],
					v[3]);
  check |= vlinterpolate->vlinterpolate(p+
					v[0]*d[4][0]+
					v[1]*d[4][1]+
					v[2]*d[4][2]+
					v[3]*d[4][3],
					v[4]);
  check |= vlinterpolate->vlinterpolate(p+
					v[0]*d[5][0]+
					v[1]*d[5][1]+
					v[2]*d[5][2]+
					v[3]*d[5][3]+
					v[4]*d[5][4],
					v[5]);

  for (int loop=0;loop<5;loop++)
    v[loop]*=s;
#endif
  return check;
}
  
void
StreamLines::FindStreamLineSupport(Array1<Point>& v /* storage for points */,
				   Point x          /* initial point */,
				   float t          /* error tolerance */,
				   float s          /* initial step size */,
				   int n            /* max number of steps */)
{
  int loop = 0;
  Array1<Vector> terms;
  Vector error;
  float err;
  int check;

  // initialize the terms
  for (loop=0;loop<6;loop++)
    terms.add(Vector(0,0,0));

  // add the initial point to the list of points found.
  v.add(x);
  std::cerr << "found point = " << x << endl;

  for (loop=0;loop<n;loop++) {

    // compute the next set of terms
    check = ComputeRKFTerms(terms,x,s);

    // compute the approximate local truncation error
    error = terms[0]*ab[0]+terms[1]*ab[1]+terms[2]*ab[2]+
            terms[3]*ab[3]+terms[4]*ab[4]+terms[5]*ab[5];
    err = sqrt(error(0)*error(0)+error(1)*error(1)+error(2)*error(2));

    // is the error tolerable?  Adjust the step size accordingly.
    if (err<t/128.0)
      s*=2;
    else if (err>t) {
      s/=2;
      loop--;         // re-do this step.
      continue;
    }

    // compute and add the point to the list of points found.
    x = x + terms[0]*a[0]+terms[1]*a[1]+terms[2]*a[2]+
            terms[3]*a[3]+terms[4]*a[4]+terms[5]*a[5];
    v.add(x);
    std::cerr << "found point = " << x << endl;
  }
}

void StreamLines::execute()
{
  d_vfport->get(d_vfhandle);
  if(!d_vfhandle.get_rep())
    return;
  
  Array1<Point> support;
  FindStreamLineSupport(support,Point(32,32,32),.01,.05,100);
}

void StreamLines::tcl_command(TCLArgs& args, void* userdata)
{
  Module::tcl_command(args, userdata);
}

} // End namespace Moulding


