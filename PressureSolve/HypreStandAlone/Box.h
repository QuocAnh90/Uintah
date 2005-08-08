#ifndef __BOX_H__
#define __BOX_H__

#include "Macros.h"
#include "Vector.h"
#include "Side.h"

typedef Vector<int> Point;

class Box {
  /*_____________________________________________________________________
    class Box:
    box extents [_lower,_upper].
    _____________________________________________________________________*/

 public:
  /*------------- Box iterators -------------*/

  //  class iterator;
  friend class iterator;

  // Const Box iterator
  class iterator {
  public:
    inline iterator(const Box& box) :
      _box(box)
      {
      }
    virtual ~iterator(void) {}
    void operator ++ (void);
    //    inline void operator ++ (int);
    //    inline void operator -- (void);
    inline Point& operator * (void)
      { return _sub; }
    inline int operator == 
      (const iterator& other) const 
      { return _sub == other._sub; }
    inline int operator !=  
      (const iterator& other) const 
      { return !(*this == other); }

  protected:
    const Box&    _box;     // Reference to attached Box object
    Point   _sub;     // Cell subscript that this iterator represents
  };

  /*============= End class Box::iterator =============*/

  Box(const Counter numDims);
  Box(const Point& lower,
      const Point& upper);
  Box(const Box& b);
  Box& operator = (const Box& b);

  const Counter getNumDims(void) const
    {
      return _lower.getLen();
    }
  const Point&  get(const Side& s) const;
  Point&        get(const Side& s);
  void          set(const Side& s,
                    const Point& value);
  void          set(const Counter d,
                    const Side& s,
                    const int& value);
  iterator      begin(void) const 
    /* Points to lower */
    { 
      iterator iter(*this);
      *iter = _lower;
      return iter;
    }

  iterator      end(void) const
    /* Points to lower - (1,0,...,0) */
    {
      iterator iter(*this);
      *iter = _lower;
      (*iter)[0] -= 1;
      return iter;
    }

  Point         size(void) const;
  Counter       volume(void) const;
  bool          overlaps(const Box&, double epsilon=1.e-6) const;
  bool          contains(const Point& p) const;
  Box           intersect(const Box& b) const;
  bool          degenerate(void) const;
  bool          degenerate(const Counter d) const;
  Box           faceExtents(const Counter d,
                            const Side& s); 
  Box           coarseNbhrExtents(const Vector<Counter>& refRat,
                                  const Counter d,
                                  const Side& s) const;
 protected:
  Point   _lower;   // Lower-left corner of box
  Point   _upper;   // Upper-right corner of box
};

std::ostream& operator << (std::ostream& os, const Box& a);

#endif // __BOX_H__
