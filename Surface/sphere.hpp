
#ifndef __SPHERE_HPP__
#define __SPHERE_HPP__

#include <cmath>
#include <memory>
#include <vector>
#include <utility>
#include <tuple>
#include <string>
#include <iostream>
#include <boost/optional.hpp>
#include <exception>

#include "frame.hpp"
#include "shape.hpp"
#include "vector3d.hpp"
#include "lvec.hpp"
#include "axisangle.hpp"

#include "plane.hpp"

using std::abs;
using std::cos;
using std::acos;
using std::sin;
using std::sqrt;
using std::vector;


// This is an 'oriented' sphere because shape has a frame object which is oriented
// AnyPlane has become the baseclass for any surface which the cell sits on so we have to subclass from it
class Sphere: public AnyPlane
{
public:
  // inherits a frame from Shape
  // sphere radius R
  double R;

  using AnyPlane::AnyPlane;

  //
  Sphere() : AnyPlane() { R = 0; }
  // sphere at origin with radius R
  Sphere(double R) : AnyPlane(), R(R) {} 

  // 
  Sphere(Vector3d origin, double R) : AnyPlane(origin), R(R) {}
  

  // pt is in sphere
  virtual bool contains(Vector3d pt);

  // pt is on sphere
  bool on_surface(Vector3d pt);
  
  // radial vector from origin to pt
  Vector3d radv(Vector3d pt);

  // normal vector at pt
  virtual Vector3d normal(Vector3d pt);

  // use spherical polar coordinates to get a point 
  Vector3d spherical(double theta, double phi);

  using AnyPlane::intersects;
  virtual boost::optional<Vector3d> intersects(Lvec& lv);

  // capsule intersection
  virtual vector<overlap> overlap_vector(Capsule& body);

  double circ() { return 2 * M_PI * R; }
  double area() { return 4 * M_PI * R*R; }
  
  virtual std::string __str__();
};

// a sphere but the surface normal is inverted
// In other words a spherical cavity
// can have two surface contacts
class Shell: public Sphere
{
  public:
  using Sphere::Sphere;
  
  Vector3d normal(Vector3d pt);    

  virtual Vector3d distance_vector(Vector3d pt) override;

  virtual vector<overlap> overlap_vector(Capsule& body) override;

};



Vector3d mirroru(Vector3d v, const Vector3d u);

// Capsule or Spherocylinder is part of a class of Sphere-Swept-Volumes (SSV)
// for which luckly collision detection is not too complicated
// (unlike with cylinder where collisions are complected because of hard edges)
// Like all subclasses of Shape this is oriented so it has a head and tail sphere
class Capsule: public Sphere
{
  public:
  // length of the swept line. Not total length from tip to tip which is length + 2R
  double length;
  // will be default initialised! must only call get_Rk after set_Rk!
  Rk bodyRk; 

  Capsule() : Capsule(Vector3d(), e_z, 0.5, 2.) {}

  Capsule(Vector3d origin, Vector3d orient, double R, double length)
    : Sphere(R), length(length)
  {
    frame = Frame(origin, orient);
  }
  // default R, length and axis rotation. useful for debugging geometry
  Capsule(Vector3d origin, Vector3d orient) : Capsule(origin, orient, 0.5, 2.) {}

  // specify the frame exactly // axis is e3
  Capsule(Vector3d origin, Vector3d axis, Vector3d e1, double R, double length)
    : Sphere(R), length(length)
  {
    Vector3d e2 = axis.cross(e1);
    frame = Frame(origin, e1, e2);
  }

  Capsule(Frame _f, double R, double length) : Sphere(R), length(length) {
    frame = _f;
  }


  ///////////////////////////////////////////
  
  void set_Rk(Rk p) {
    this->bodyRk = p;
  }
  Rk get_Rk(void) {
    return this->bodyRk;
  }
  void init_Rk(void) {
    this->set_Rk(Rk(get_axisangle()));
  }

  ///////////////////////////////////////////

  Vector3d get_axis() {
    return frame.e3;  // Check this remains a unit vector, despite numerical drift.
  }

  Vector3d get_vector() {
    return -length * get_axis();
  }

  Vector3d get_caxis() {
    return 0.5 * length * get_axis();
  }

  //
  Vector3d get_endpt() {
    return frame.origin - get_caxis();
  }
  Vector3d get_headpt() {
    return frame.origin + get_caxis();
  }

  Vector3d get_centerpt() {
    return frame.origin;
  }


  // head -> tail
  Lvec get_lvec();

  // construct an axis angle for the current axis using e_z as reference
  Vector3d get_axisangle();
  // call once before minimisation
  
  void set_rotation(Vector3d p);

  bool contains(Vector3d pt)
  { throw std::runtime_error("Not Implemented"); } 

  // 
  Vector3d normal(Vector3d pt);
  // assume pt is on capsule head
  Vector3d head_normal(Vector3d pt);

  virtual boost::optional<Vector3d> intersects(Lvec& lv) 
  { throw std::runtime_error("Not Implemented"); } 
  virtual boost::optional<Vector3d> intersects(Chain& ch)
  { throw std::runtime_error("Not Implemented"); } 

  //bool intersects(std::shared_ptr<Plane> plane);

  // distance vector is similar to Lvec.distance vector however
  // here we make the assumption that the plane does not intersect with 
  // the core line segment of the capsule
  
  // later in the MD step we will assume the timestep is sufficienly small
  // that distance vector always points to one end of the Capsule core.
  // We should be able to ignore the critical case of capsule exactly horizontal

  // In general we want to use the overlap distance to compute a potential
  // so return it along with a unit vector indicating how to reverse
  // the overlap direction
  // If we want use polymorphism like this, the overlap methods need to be 
  // part of the surface classes, so move this to plane.hpp
  //vector<overlap> overlap_vector(std::shared_ptr<Plane> plane);

  boost::optional<overlap> intersects(Sphere& sph);


  // overlap vector starting be considering only the Head and tail sphere
  //vector<overlap> overlap_vector(Step& step);

  // point on the head of the capsule with 
  //
  // 0 < theta < \pi/2
  // 0 < \phi < 2\pi
  // used? dep?
  Vector3d alt_spherical(double theta, double phi);

  double long_length();
  double true_length();
  Vector3d mirror(Vector3d v);
  // return a pilus anchor point just inside the surface
  std::pair<Vector3d, Vector3d> surface_pt(double d, double phi, double inlength);
  Vector3d spt(double d, double phi);
  int get_segment(double d);
  double halfal() { return length/2. + this->circ()/4.; }
  // for t = 0, m = length/2., for t = 1 ,m = -length/2.
  double get_m(double t) { return this->length * (0.5-t); }
  double get_t(double m) { return 0.5 - m/this->length; }
  Vector3d get_pt(double m) { return get_origin() + m * get_axis(); }

  // sphere area
  double spharea() { return 4 * M_PI * R*R; }

  double sgn_distance_to_surface(Vector3d pt);

  std::string __str__();

  static int HEAD;
  static int BODY;
  static int TAIL;
};

std::ostream& operator<<(std::ostream& os, Capsule& cap);

namespace Caps {
  Capsule headat(Vector3d, Vector3d);
  Capsule tailat(Vector3d, Vector3d);
  Capsule body_from_lv(Lvec);
}

#endif
