
/*
 *  DistanceConstraint.h
 *
 *  Written by:
 *   James Purciful
 *   Department of Computer Science
 *   University of Utah
 *   Aug. 1994
 *
 *  Copyright (C) 1994 SCI Group
 */


#include <Constraints/DistanceConstraint.h>
#include <Geometry/Vector.h>
#include <Classlib/Debug.h>

static DebugSwitch dc_debug("BaseConstraint", "Distance");

DistanceConstraint::DistanceConstraint( const clString& name,
					const Index numSchemes,
					Variable* p1, Variable* p2,
					Variable* distInX )
:BaseConstraint(name, numSchemes, 3), guess(1, 0, 0)
{
   vars[0] = p1;
   vars[1] = p2;
   vars[2] = distInX;
   whichMethod = 0;

   // Tell the variables about ourself.
   Register();
}

DistanceConstraint::~DistanceConstraint()
{
}


void
DistanceConstraint::Satisfy( const Index index, const Scheme scheme )
{
   Variable& v0 = *vars[0];
   Variable& v1 = *vars[1];
   Variable& v2 = *vars[2];
   Vector v;

   if (dc_debug) {
      ChooseChange(index, scheme);
      printc(cout, scheme);
   }
   
   /* Q <- Sc + Sr * Normalize(P-Sc) */
   switch (ChooseChange(index, scheme)) {
   case 0:
      v = (v0.Get() - v1.Get());
      if (v.length2() == 0.0)
	 v = guess;
      else
	 v.normalize();
      v0.Assign(v1.Get() + (v * v2.Get().x()),
		scheme);
      break;
   case 1:
      v = (v1.Get() - v0.Get());
      if (v.length2() == 0.0)
	 v = guess;
      else
	 v.normalize();
      v1.Assign(v0.Get() + (v * v2.Get().x()),
		scheme);
      break;
   case 2:
      v2.Assign(Point((v1.Get() - v0.Get()).length(), 0, 0),
		scheme);
      break;
   default:
      cerr << "Unknown variable in Distance Constraint!" << endl;
      break;
   }
}

