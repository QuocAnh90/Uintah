/*
 *  Graph.h:
 *
 *  Written by:
 *   Yarden Livnat
 *   July 20001
 *
 */


#ifndef Graph_h
#define Graph_h

#include <string>
#include <map>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

#include <Core/Thread/ConditionVariable.h>
#include <Core/2d/DrawObj.h>
#include <Core/Datatypes/Datatype.h>
#include <Core/Containers/Array1.h>
#include <Core/Malloc/Allocator.h>
#include <Core/GuiInterface/TCL.h>
#include <Core/GuiInterface/TCLTask.h>
#include <Core/GuiInterface/GuiVar.h>
#include <Core/GuiInterface/TclObj.h>
#include <Core/2d/DrawGui.h>
#include <Core/Thread/Mutex.h>
#include <Core/Geom/Color.h>
#include <Core/2d/OpenGLWindow.h>



namespace SCIRun {

class GraphHelper;

struct ObjInfo {
  string name_;
  DrawGui *obj_;
  bool mapped_;

  ObjInfo( const string &name, DrawGui *d);
  ~ObjInfo() { if (obj_) delete obj_; }

  void draw() { obj_->draw(); }
  void set_window( const string &);
  void set_id( const string & );
};

class Graph : public TclObj, public DrawObj {
  friend GraphHelper;
private:
  ObjInfo *obj_;
  OpenGLWindow *ogl_;
  GuiString gl_window;

  GraphHelper *helper_;
  ConditionVariable has_work_;

public:
  Graph( const string & );
  virtual ~Graph() {}

  void lock() { TCLTask::lock(); }
  void unlock() { TCLTask::unlock(); }

  void add( const string &, DrawGui *);
  virtual void tcl_command(TCLArgs&, void*);
  virtual void set_window( const string &);
  void update();

  virtual void need_redraw();
  virtual void get_bounds( BBox2d &) {}
  virtual void draw() {}
  virtual void io(Piostream& stream);

};


} // End namespace SCIRun

#endif Graph_h


