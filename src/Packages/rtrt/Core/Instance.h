#ifndef INSTANCE_H
#define INSTANCE_H

#include <Packages/rtrt/Core/Object.h>
#include <Packages/rtrt/Core/Ray.h>
#include <Packages/rtrt/Core/Material.h>
#include <Packages/rtrt/Core/PerProcessorContext.h>
#include <Packages/rtrt/Core/HitInfo.h>
#include <Packages/rtrt/Core/InstanceWrapperObject.h>

#include <Core/Geometry/Vector.h>

// Known problems:
// Scales may not work because ray parameter isn't transformed --
// the direction vector scale is adjusted instead.
// Why are all the virtual functions in this file?  We should make it's
//  very own .cc file

namespace rtrt {
  class Instance;
}

namespace SCIRun {
  void Pio(Piostream&, rtrt::Instance*&);
}

namespace rtrt {

class Instance: public Object, public Material
{
public:
  struct InstanceHit {
    Vector normal;
    Object* obj;
  };

  InstanceWrapperObject * o;
  Transform             * currentTransform;
  BBox                    bbox;

  Instance(InstanceWrapperObject* o, Transform* trans) 
    : Object(this), o(o)
  {
    currentTransform = new Transform();
    *currentTransform = *trans;

    if( !currentTransform->inv_valid() ) 
      currentTransform->compute_imat();

    o->compute_bounds(bbox,1E-5);

    bbox.transform_inplace(currentTransform);
  }

  Instance(InstanceWrapperObject* o, Transform* trans, BBox& b) 
    : Object(this), o(o)
  {
    currentTransform = new Transform();
    *currentTransform = *trans;

    if (!currentTransform->inv_valid())
      currentTransform->compute_imat();

    bbox = b.transform(currentTransform);
  }
  Instance() : Object(0), Material() {} // for Pio.

  //! Persistent I/O.
  static  SCIRun::PersistentTypeID type_id;
  virtual void io(SCIRun::Piostream &stream);
  friend void SCIRun::Pio(SCIRun::Piostream&, Instance*&);

  virtual void intersect(Ray& ray, HitInfo& hit, DepthStats* st,
			 PerProcessorContext* ppc)
  {
    double min_t = hit.min_t;
    if (!bbox.intersect(ray, min_t)) return;	  

    Ray tray;

    ray.transform(currentTransform,tray);
    Vector td = tray.direction();
    double scale = td.normalize();
    tray.set_direction(td);

    HitInfo thit;
    if (hit.was_hit) thit.min_t = hit.min_t * scale;

    o->intersect(tray,thit,st,ppc);
	  
    // if the ray hit one of our objects....
    if (thit.was_hit)
      {
	min_t = thit.min_t / scale;
	if(hit.hit(this, min_t)){
	  InstanceHit* i = (InstanceHit*)(hit.scratchpad);
	  Point p = tray.origin() + thit.min_t*tray.direction();	  
	  i->normal = thit.hit_obj->normal(p,thit);
	  i->obj = thit.hit_obj;
	}
      }	      
  }
    
  virtual Vector normal(const Point&, const HitInfo& hit)
  {
    InstanceHit* i = (InstanceHit*)(hit.scratchpad);
    Vector n;
    currentTransform->project_normal(i->normal, n);
    n.normalize();
    return n;
  }

  virtual void compute_bounds(BBox& b, double /*offset*/)
  {
    b.extend(bbox);
  }

  virtual void preprocess(double maxradius, int& pp_offset, int& scratchsize)
  {
    o->preprocess(maxradius,pp_offset,scratchsize);
  }

  virtual void shade(Color& result, const Ray& ray,
		     const HitInfo& hit, int depth,
		     double atten, const Color& accumcolor,
		     Context* cx) {
    InstanceHit* i = (InstanceHit*)(hit.scratchpad);
    Material *mat = i->obj->get_matl();
    mat->shade(result, ray, hit, depth, atten, accumcolor, cx);
  }

  virtual void animate(double time, bool& changed) {
    o->animate(time, changed);
  }

  virtual bool interior_value(double& ret_val, const Ray &ref, const double _t)
  {

    Ray tray;

    ref.transform(currentTransform,tray);
    Vector td = tray.direction();
    double scale = td.normalize();
    tray.set_direction(td);
    double _t2 = _t * scale;

    return o->interior_value(ret_val, tray, _t2);
  }

}; // end class Instance

} // end namespace rtrt

#endif
