
/*
 *  Assert.h: Utility for specifying data invariants (Assertions)
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   Feb. 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#ifndef SCI_Classlib_Assert_h
#define SCI_Classlib_Assert_h 1

#ifdef SCI_ASSERTIONS
#ifndef SCI_ASSERTION_LEVEL
#define SCI_ASSERTION_LEVEL 2
#endif

#include "Exceptions.h"

#if SCI_ASSERTION_LEVEL >= 1
#define ASSERTL1(condition) \
	if(!(condition)){ \
		AssertionFailed exc(#condition); \
		EXCEPTION(exc); \
	}
#else
#define ASSERTL1(condition)
#endif

#if SCI_ASSERTION_LEVEL >= 2
#define ASSERTL2(condition) \
	if(!(condition)){ \
		AssertionFailed exc(#condition); \
		EXCEPTION(exc); \
	}
#else
#define ASSERTL2(condition)
#endif

#if SCI_ASSERTION_LEVEL >= 3
#define ASSERTL3(condition) \
	if(!(condition)){ \
		AssertionFailed exc(#condition); \
		EXCEPTION(exc); \
	}
#else
#define ASSERTL3(condition)
#endif

#if SCI_ASSERTION_LEVEL >= 4
#define ASSERTL4(condition) \
	if(!(condition)){ \
		AssertionFailed exc(#condition); \
		EXCEPTION(exc); \
	}
#else
#define ASSERTL4(condition)
#endif

#if SCI_ASSERTION_LEVEL >= 2
#define ASSERTEQ(c1, c2) \
	if(c1 != c2){ \
		AssertionEQFailed exc(#c1, #c2, (int)c1, (int)c2); \
		EXCEPTION(exc); \
	}
#define ASSERTRANGE(c, l, h) \
        if(c < l || c >= h){ \
		AssertionRangeFailed exc(#c, #l, #h, c, l, h); \
		EXCEPTION(exc); \
        }
#else
#define ASSERTEQ(c1, c2)
#define ASSERTRANGE(c, l, h)
#endif

#else

#define ASSERTL1(condition)
#define ASSERTL2(condition)
#define ASSERTL3(condition)
#define ASSERTL4(condition)
#define ASSERTEQ(c1, c2)
#define ASSERTRANGE(c, l, h)

#endif

#define ASSERT(condition) ASSERTL2(condition)

#endif
