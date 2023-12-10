
#ifndef __PLANE_HPP__
#define __PLANE_HPP__

#include <string>
#include <iostream>

#include <boost/optional.hpp>
#include <boost/function.hpp>

#include <Eigen/Dense>

#include "frame.hpp"
#include "shape.hpp"
#include "matrix3d.hpp"
#include "lvec.hpp"
#include "chain.hpp"

using std::cout;
using std::endl;

typedef boost::function<double(double)> contactf;

class Capsule; // forward definition
class ACell;

// our place object is actually a halfspace
//
//class Plane // tmp : decoupled Plane from Shape
// geometry methods here are implemented taking advantage of plane with normal e_z, origin 0
// Need a new object to override with generic implementations for a plane with arbitrary normal
class Plane: public Shape
{
  public:
  Plane() { this->frame = Frame(); }

  Plane(const Plane& pl) { frame = pl.frame; }

  Vector3d get_origin() { return frame.origin; } 


  // for generic planes embedded in 3d space
  Plane(Vector3d origin, Vector3d normal) {
    this->frame = Frame(origin, normal);
  }

  Plane(Vector3d origin) : Plane(origin, e_z) {}


  // use N() for constant normal. i.e. for Plane
  virtual Vector3d N() { return frame.e3; }
  // use normal() if the normal is position dependent
  virtual Vector3d normal(Vector3d pt) { return this->N(); }

  virtual double distance(Vector3d pt) { return (pt - get_origin()).dot(N()); }

  // project pt onto plane along the shortest path (along normal axis)
  virtual Vector3d project(Vector3d pt) { return pt - distance(pt) * N(); }
  virtual Vector3d distance_vector(Vector3d pt) { return distance(pt) * N(); }

  virtual bool contains(Vector3d pt) { return this->distance(pt) < 0.; }
  virtual bool containsany(vector<Vector3d> pts);

  virtual boost::optional<Vector3d> intersects(Lvec& lv);
  virtual boost::optional<Vector3d> intersects(Chain& ch);

  // bacterium body -> surface collision
  virtual vector<overlap> overlap_vector(Capsule& body);
  virtual int get_num_contacts(Capsule& body);
  virtual double surface_energy(ACell& cell, contactf energyf, bool cost_anchor_intersection=0);
  virtual Eigen::RowVectorXf surface_grad(Capsule& body, contactf contact);
  virtual Eigen::RowVectorXf surface_grad_part(Capsule& body, contactf contact, overlap over);

  // Allow subclasses to store temporary information between
  // calls to surface_energy and surface_grad, etc.
  virtual void reset() {/* do nothing */ }

  virtual std::string report_touching() { return "NOTIMPLEMENTED"; }
  virtual std::string __str__();
  virtual std::string _class();

  protected:
  virtual double _surface_energy(ACell& cell, contactf energyf);

};

std::ostream& operator<<(std::ostream& os, Plane pl);


// AnyPlane should be the baseclass and Plane should derive, I don't have time to refactor.
class AnyPlane: public Plane 
{
  public:
  using Plane::Plane;

  // override methods which are not entirely general
  virtual vector<overlap> overlap_vector(Capsule& body);
  int get_num_contacts(Capsule& body);
  virtual Eigen::RowVectorXf surface_grad(Capsule& body, contactf contact);

  //
  virtual std::string _class();

  protected:
  virtual double _surface_energy(ACell& cell, contactf energyf);
};

// These partial plane classes have horribly duplicated code because of limitations of Vector3d
// and efficiency concerns


class PartialPlane: public AnyPlane
{
  public:
/* generic partial plane that is infinite in one direction and has limits in another */
  double lmin, lmax; // limits in the limited direction

  PartialPlane(Vector3d origin, Vector3d normal, Vector3d e1, double lmin, double lmax)
    : lmin(lmin), lmax(lmax)
  {
    // let e2 be the infinite direction, e1 is the finite direction
    Vector3d e2 = normal.cross(e1);
    this->frame = Frame(origin, e1, e2, normal);
  }

  // for a point in the plane, project onto limited direction and get coordinate x
  // such that lmin < x < lmax means the point is contained
  double proj_limited(Vector3d pt) { return (pt - get_origin()).dot(frame.e1); }
  bool on_plane(Vector3d pt);

  boost::optional<Vector3d> intersects(Lvec& lv) override;
  virtual vector<overlap> overlap_vector(Capsule& body) override;

  // get a point in the limited direction
  Vector3d get_lmpt(void) { return get_origin() + lmax*frame.e1; } 

};

class PartialPlaneX: public AnyPlane
{
  public:
  double xmin, xmax;
  bool is_upper;
  using AnyPlane::AnyPlane;

  PartialPlaneX(Vector3d origin, Vector3d normal, double xmin, double xmax, bool is_upper=false)
    : AnyPlane(origin, normal), xmin(xmin), xmax(xmax), is_upper(is_upper)
  {}

  boost::optional<Vector3d> intersects(Lvec& lv) override;
  virtual vector<overlap> overlap_vector(Capsule& body) override;
  
  std::string _class();
  std::string __str__() override;

};

class PartialPlaneZ: public AnyPlane
{
  public:
  double zmin, zmax;
  using AnyPlane::AnyPlane;
  PartialPlaneZ(Vector3d origin, Vector3d normal, double zmin, double zmax)
    : AnyPlane(origin, normal), zmin(zmin), zmax(zmax)
  {}

  boost::optional<Vector3d> intersects(Lvec& lv) override;
  virtual vector<overlap> overlap_vector(Capsule& body) override;
  // special
  vector<overlap> end_overlap(Capsule& body);

  std::string _class();
  std::string __str__() override;
};



// for debugging create a surface object which does not contribute energy or gradient
// and cannot be attached to
class NullPlane: public Plane
{
  using Plane::Plane;

  boost::optional<Vector3d> intersects(Lvec& lv) override { return boost::none; }
  boost::optional<Vector3d> intersects(Chain& ch) override { return boost::none; }

  virtual double surface_energy(ACell& cell, contactf energyf,
      bool cost_anchor_intersection = 0) { return 0.; }
  virtual Eigen::RowVectorXf surface_grad(Capsule& body, contactf contact) 
  { return Eigen::RowVectorXf::Zero(6); }

  std::string __str__() override;
};


#endif
