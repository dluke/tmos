
#ifndef __SHAPE_HPP__
#define __SHAPE_HPP__


#include <iostream>
#include <tuple>
#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "frame.hpp"
#include "lvec.hpp"
#include "chain.hpp"

#include "boost/optional/optional.hpp"

// abstract base class for shapes
// i.e. sphere, plane, cylinder, 

typedef std::tuple<double, double, Vector3d> overlap;
// distance along axis, m
// 
// Contact point on the capsule surface to calculate torque

std::ostream& operator<<(std::ostream& os, const overlap& over);

// must be abstract base class
class Shape
{
  protected: Shape() : frame(Frame()) {}

  // making frame private and writing getters and setters for various
  // values is technically the better solution but involves writing more code
  public:
  Frame frame;

  virtual void translate(Vector3d r) {
    this->frame.translate(r);
  }

  virtual void rotate(Matrix3d M) {
    this->frame.rotate(M);
  }

  // because all shapes atleast have an reference pt called origin
  virtual Vector3d get_origin() {
    return frame.origin;
  }
  virtual Vector3d set_origin(Vector3d orv) {
    this->frame.origin = orv;
  }
  virtual void set_rotation(Vector3d p) {
    this->frame.set_rotation(p);
  }

  virtual Frame get_frame() { return frame; }

  //abstract methods
  virtual bool contains(Vector3d pt) = 0;

  // required for Straight Pili attachment
  virtual boost::optional<Vector3d> intersects(Lvec& sh) = 0;
  // required for WLC pili attachment
  virtual boost::optional<Vector3d> intersects(Chain& ch);
  
  //return the normal vector at a point on the surface
  virtual Vector3d normal(Vector3d pt) = 0;

  virtual std::string __str__();
  
};

#endif

