
/*
 *  PythagorasConstraint.h
 *
 *  Written by:
 *   James Purciful
 *   Department of Computer Science
 *   University of Utah
 *   Aug. 1994
 *
 *  Copyright (C) 1994 SCI Group
 */


#ifndef SCI_project_Pythagoras_Constraint_h
#define SCI_project_Pythagoras_Constraint_h 1

#include <SCICore/share/share.h>
#include <PSECore/Constraints/BaseConstraint.h>

namespace PSECore {
namespace Constraints {

class SCICORESHARE PythagorasConstraint : public BaseConstraint {
public:
   PythagorasConstraint( const clString& name,
			 const Index numSchemes,
			 RealVariable* dist1, RealVariable* dist2,
			 RealVariable* hypo );
   virtual ~PythagorasConstraint();
   
protected:
   virtual int Satisfy( const Index index, const Scheme scheme, const Real Epsilon,
			BaseVariable*& var, VarCore& c );
};

} // End namespace Constraints
} // End namespace PSECore

//
// $Log$
// Revision 1.2  1999/08/17 06:38:19  sparker
// Merged in modifications from PSECore to make this the new "blessed"
// version of SCIRun/Uintah.
//
// Revision 1.1  1999/07/27 16:55:56  mcq
// Initial commit
//
// Revision 1.3  1999/05/06 20:17:07  dav
// added back PSECore .h files
//
// Revision 1.1.1.1  1999/04/24 23:12:53  dav
// Import sources
//
//

#endif
