
/*
 *  ArrowWidget.h
 *
 *  Written by:
 *   James Purciful
 *   Department of Computer Science
 *   University of Utah
 *   Aug. 1994
 *
 *  Copyright (C) 1994 SCI Group
 */


#ifndef SCI_project_Arrow_Widget_h
#define SCI_project_Arrow_Widget_h 1

#include <PSECore/Widgets/BaseWidget.h>

namespace SCICore {
  namespace GeomSpace {
    class GeomPick;
  }
}

namespace PSECore {
namespace Widgets {

  //using SCICore::GeomSpace::GeomPick;

using PSECore::Dataflow::Module;

class ArrowWidget : public BaseWidget {
public:
   ArrowWidget( Module* module, CrowdMonitor* lock, double widget_scale );
   ArrowWidget( const ArrowWidget& );
   virtual ~ArrowWidget();

   virtual void redraw();
   virtual void geom_moved(GeomPick*, int, double, const Vector&, int, const BState&);

   virtual void MoveDelta( const Vector& delta );
   virtual Point ReferencePoint() const;

   void SetPosition( const Point& );
   Point GetPosition() const;
   
   void SetDirection( const Vector& v );
   const Vector& GetDirection() const;

   virtual void widget_tcl( TCLArgs& );

   // Variable indexs
   enum { PointVar };

   // Material indexs
   enum { PointMatl, ShaftMatl, HeadMatl };

protected:
   virtual clString GetMaterialName( const Index mindex ) const;   
   
private:
   Vector direction;
};

} // End namespace Widgets
} // End namespace PSECore

//
// $Log$
// Revision 1.2  1999/08/17 06:38:27  sparker
// Merged in modifications from PSECore to make this the new "blessed"
// version of SCIRun/Uintah.
//
// Revision 1.1  1999/07/27 16:56:05  mcq
// Initial commit
//
// Revision 1.3  1999/05/06 20:17:22  dav
// added back PSECore .h files
//
// Revision 1.1.1.1  1999/04/24 23:12:29  dav
// Import sources
//
//


#endif
