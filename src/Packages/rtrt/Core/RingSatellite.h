
#ifndef RINGSATELLITE_H
#define RINGSATELLITE_H 1

#include <Packages/rtrt/Core/Ring.h>
#include <Packages/rtrt/Core/Satellite.h>
#include <stdlib.h>

namespace rtrt {
class RingSatellite;
}

namespace SCIRun {
void Pio(Piostream&, rtrt::RingSatellite*&);
}

namespace rtrt {

class RingSatellite : public Ring
{

 protected:

  Satellite *parent_;

 public:

  RingSatellite(const string &name, Material *mat, const Point &center,
                const Vector &up, double radius, double thickness,
                Satellite *parent=0) 
    : Ring(mat, center, up, radius, thickness), parent_(parent)
  {
    Names::nameObject(name, this);

    if (parent_) 
      cen = parent->get_center();
  }
  virtual ~RingSatellite() {}
  RingSatellite() : Ring() {} // for Pio.

  virtual void uv(UV& uv, const Point&, const HitInfo& hit);

  //! Persistent I/O.
  static  SCIRun::PersistentTypeID type_id;
  virtual void io(SCIRun::Piostream &stream);
  friend void SCIRun::Pio(SCIRun::Piostream&, Ring*&);

  Satellite *get_parent() { return parent_; }
  void set_parent(Satellite *p) { parent_ = p; }

  double get_radius() const { return radius; }
  void set_radius(double r) { radius = r; }

  Point &get_center() { return cen; }
  void set_center(const Point &p) { cen = p; }

  virtual void compute_bounds(BBox& bbox, double offset)
  {
    if (parent_) {
#if 0
      Point center = parent_->get_center();
      bbox.extend(center);
      Point extent = 
        Point(center.x()+parent_->get_orb_radius()+radius+offset,
              center.y()+parent_->get_orb_radius()+radius+offset,
              center.z());
      bbox.extend( extent );
      extent = 
        Point(center.x()-(parent_->get_orb_radius()+radius+offset),
              center.y()-(parent_->get_orb_radius()+radius+offset),
              center.z());
      bbox.extend( extent );

      extent = Point(center.x(),center.y(),
                     center.z()+parent_->get_radius()+offset);
      bbox.extend( extent );

      extent = Point(center.x(),center.y(),
                     center.z()-(parent_->get_radius()+offset));
      bbox.extend( extent );
#else
      bbox.extend(cen,radius+offset);
#endif
    } else {
      bbox.extend(cen, radius+thickness+offset);
    }
  }

  virtual void animate(double t, bool& changed);
};

} // end namespace

#endif
