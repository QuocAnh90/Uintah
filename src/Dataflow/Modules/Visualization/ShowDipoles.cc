/*
 *  ShowDipoles.cc:  Builds the RHS of the FE matrix for current sources
 *
 *  Written by:
 *   David Weinstein
 *   University of Utah
 *   May 1999
 *
 *  Copyright (C) 1999 SCI Group
 */

#include <Dataflow/Network/Module.h>
#include <Dataflow/Ports/ColumnMatrixPort.h>
#include <Dataflow/Ports/GeometryPort.h>
#include <Dataflow/Ports/MatrixPort.h>
#include <Dataflow/Ports/MeshPort.h>
#include <Dataflow/Ports/SurfacePort.h>
#include <Dataflow/Widgets/ArrowWidget.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Geom/GeomLine.h>
#include <Core/Geom/Switch.h>
#include <Core/Geometry/Point.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Trig.h>
#include <Core/TclInterface/TCLvar.h>
#include <iostream>
using std::cerr;

namespace SCIRun {

class ShowDipoles : public Module {
  MatrixIPort *imat;
  MatrixOPort *omat;
  ColumnMatrixIPort *icol;
  ColumnMatrixOPort *ocol;
  GeometryOPort* ogeom;
  int gen;
  MatrixHandle dipoleMatH;
  TCLstring widgetSizeTCL;
  TCLstring scaleModeTCL;
  TCLint showLastVecTCL;
  TCLint showLinesTCL;
  int which;
  double lastSize;
  clString execMsg;
  Array1<GeomSwitch *> widget_switch;
public:
  ShowDipoles(const clString& id);
  virtual ~ShowDipoles();
  virtual void execute();
  CrowdMonitor widget_lock;
  Array1<int> widget_id;
  Array1<ArrowWidget*> widget;
  int gidx;

  MaterialHandle greenMatl;
  MaterialHandle deflMatl;

  virtual void widget_moved(int last);
  int nDips;
};

extern "C" Module* make_ShowDipoles(const clString& id)
{
  return scinew ShowDipoles(id);
}

ShowDipoles::ShowDipoles(const clString& id) : 
  Module("ShowDipoles", id, Filter), 
  widgetSizeTCL("widgetSizeTCL", id, this),
  widget_lock("ShowDipoles widget lock"),
  scaleModeTCL("scaleModeTCL", id, this),
  showLastVecTCL("showLastVecTCL", id, this),
  showLinesTCL("showLinesTCL", id, this)
{
  // Create the input port
  imat=scinew MatrixIPort(this, "DipoleMatrix", MatrixIPort::Atomic);
  add_iport(imat);

  // Create the output ports
  omat=scinew MatrixOPort(this, "DipoleMatrix", MatrixIPort::Atomic);
  add_oport(omat);
  ogeom=scinew GeometryOPort(this,"Geometry",GeometryIPort::Atomic);
  add_oport(ogeom);
  gen=-1;
  nDips=0;
  lastSize=-1;
  greenMatl=new Material(Color(0.2, 0.8, 0.2));
  gidx=0;
}

ShowDipoles::~ShowDipoles()
{
}

void ShowDipoles::execute()
{
  MatrixHandle mh;
  Matrix* mp;
  if (!imat->get(mh) || !(mp=mh.get_rep())) {
    cerr << "No input in ShowDipoles Matrix port.\n";
    return;
  }
  cerr << "nrows="<<mp->nrows()<<"  ncols="<<mp->ncols()<<"\n";
  if (mp->ncols() != 6) {
    cerr << "Error - dipoles must have six entries.\n";
    return;
  }
  double widgetSize;
  if (!widgetSizeTCL.get().get_double(widgetSize)) {
    widgetSize=1;
    widgetSizeTCL.set("1.0");
  }
     
  if (mh->generation != gen || lastSize != widgetSize) {// load this data in
    if (mp->nrows() != nDips) {
	     
      cerr << "NEW SIZE FOR DIPOLEMATTOGEOM  mp->nrows()="<<mp->nrows()<<" nDips="<<nDips<<"\n";
	     
      // nDips always just says how many switches we have set to true
      // need to fix switch setting first and then do allocations if
      //   necessary
	     
      if (widget_switch.size()) {
	widget[nDips-1]->SetCurrentMode(0);
	widget[nDips-1]->SetMaterial(0, deflMatl);
      }
      if (mp->nrows()<nDips) {
	for (int i=mp->nrows(); i<nDips; i++)
	  widget_switch[i]->set_state(0);
	nDips=mp->nrows();
      } else {
	int i;
	for (i=nDips; i<widget_switch.size(); i++)
	  widget_switch[i]->set_state(1);
	for (; i<mp->nrows(); i++) {
	  widget.add(scinew ArrowWidget(this, &widget_lock, widgetSize));
	  deflMatl=widget[0]->GetMaterial(0);
	  widget_switch.add(widget[i]->GetWidget());
	  widget_switch[i]->set_state(1);
	  widget_id.add(ogeom->addObj(widget_switch[i], clString(clString("Dipole")+to_string(i)), &widget_lock));
	}
	nDips=mp->nrows();
      }
      if (showLastVecTCL.get()) {
	widget[nDips-1]->SetCurrentMode(0);
	widget[nDips-1]->SetMaterial(0, deflMatl);
      } else {
	widget[nDips-1]->SetCurrentMode(2);
	widget[nDips-1]->SetMaterial(0, greenMatl);
      }
    }
    Array1<Point> pts;
    int i;
    clString scaleMode=scaleModeTCL.get();
    double max;
    for (i=0; i<mp->nrows(); i++) {
      double dv=Vector((*mp)[i][3], (*mp)[i][4], (*mp)[i][5]).length();
      if (dv<0.00000001) dv=1;
      if (i==0 || dv<max) max=dv;
    }

    for (i=0; i<mp->nrows(); i++) {
      Point p((*mp)[i][0], (*mp)[i][1], (*mp)[i][2]);
      pts.add(p);
      widget[i]->SetPosition(p);
      Vector v((*mp)[i][3], (*mp)[i][4], (*mp)[i][5]);
      //	     cerr << "widget["<<i<<"] is at position "<<p<<" and dir "<<v<<"\n";
      double str=v.length();
      if (str<0.0000001) v.z(1);
      v.normalize();
      widget[i]->SetDirection(v);
      //	     widget[i]->SetScale(str*widgetSize);
      //	     widget[i]->SetScale(widgetSize);
      double sc=widgetSize;
      if (scaleMode == "normalize") sc*=(str/max);
      else if (scaleMode == "scale") sc*=str;
      widget[i]->SetScale(sc);
      widget[i]->SetLength(2*sc);
    }

    if (gidx) ogeom->delObj(gidx);
    if (showLinesTCL.get()) {
      GeomLines *g=new GeomLines;
      for (i=0; i<pts.size()-2; i++) 
	for (int j=i+1; j<pts.size()-1; j++) 
	  g->add(pts[i], pts[j]);
      GeomMaterial *gm=new GeomMaterial(g, new Material(Color(.8,.8,.2)));
      gidx=ogeom->addObj(gm, clString("Dipole Lines"));
    }

    gen=mh->generation;
    dipoleMatH=mh;
    lastSize=widgetSize;
    ogeom->flushViews();
    omat->send(dipoleMatH);
    //     } else if (execMsg == "widget_moved") {
    //	 cerr << "Can't handle widget_moved callbacks yet...\n";
  } else if (execMsg == "widget_moved") {
    execMsg="";
    Array1<Point> pts;
    int i;
    for (i=0; i<nDips; i++) {
      Point p=widget[i]->GetPosition();
      pts.add(p);
      Vector d=widget[i]->GetDirection();
      double mag=widget[i]->GetScale();
      cerr << "mag="<<mag<<"  widgetSize="<<widgetSize<<"\n";
      d=d*(mag/widgetSize);
      (*mp)[i][0]=p.x();
      (*mp)[i][1]=p.y();
      (*mp)[i][2]=p.z();
      (*mp)[i][3]=d.x();
      (*mp)[i][4]=d.y();
      (*mp)[i][5]=d.z();
    }
    ogeom->delObj(gidx);
    if (showLinesTCL.get()) {
      GeomLines *g=new GeomLines;
      for (i=0; i<pts.size()-2; i++) 
	for (int j=i+1; j<pts.size()-1; j++) 
	  g->add(pts[i], pts[j]);
      GeomMaterial *gm=new GeomMaterial(g, new Material(Color(.8,.8,.2)));
      gidx=ogeom->addObj(gm, clString("Dipole Lines"));
    }
    ogeom->flushViews();
    dipoleMatH=mh;
    omat->send(dipoleMatH);
  } else {
    // just send the same old matrix/vector as last time
    cerr << "sending old stuff!\n";
    omat->send(dipoleMatH);
  }

  //     cerr << "ShowDipoles: Here are the dipoles...\n";
  for (int i=0; i<mp->nrows(); i++) {
    //	 cerr << "   "<<i<<"   ";
    for (int j=0; j<mp->ncols(); j++) {
      //	     cerr << (*mp)[i][j]<<" ";
    }
    //	 cerr << "\n";
  }

}

void ShowDipoles::widget_moved(int last) {
  if(last && !abort_flag) {
    abort_flag=1;
    execMsg="widget_moved";
    want_to_execute();
  }
} 
} // End namespace SCIRun
