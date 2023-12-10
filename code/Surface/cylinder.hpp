
#ifndef __CYLINDER_HPP__
#define __CYLINDER_HPP__

#include <string>
#include <map>

#include "frame.hpp"
#include "matrix3d.hpp"
#include "shape.hpp"
#include "plane.hpp"
#include "lvec.hpp"
#include "chain.hpp"

// watch out for circular references
#include "sphere.hpp"
#include "cell.hpp"

#include <boost/optional.hpp>

using std::cout;
using std::endl;

using std::sqrt;
using std::floor;
using std::ceil;
using std::atan2;

const double float_tol = 0.00000001;

// an infinite cylinder with infinite direction frame.e2
// Should have a small radius and by considered as a 90 degree arc
// for the purpose of pili attachment

// All calculations done in xz plane with y=0 
// Inherit intersects(Chain&) method from Shape + a Frame we won't use
class Circle: public Shape 
{
  public:
  double r;
  Vector3d origin;

  Circle(Vector3d origin, double r) 
    : origin(origin.xzproject()), r(r)
  {}
  
  virtual double distance(Vector3d pt);
  //virtual Vector3d distance_vector(Lvec sh);
  virtual boost::optional<Vector3d> intersects(Lvec& sh);

  //-- Useful methods for working with arcs
  // Use atan2 to get the angle corresponding to a surface point 
  virtual double calc_theta(Vector3d pt);
  // An extra check
  virtual double on_surface(Vector3d pt);

  // Boring virtual methods need implementations
  virtual Vector3d normal(Vector3d surfacept) { 
    return surfacept - origin;
  }
  virtual bool contains(Vector3d pt) { return distance(pt) < r; }

  // get/set
  double get_r() { return r; }
  Vector3d get_origin() { return origin; }
  
};

// all calculations for this object can be done in xz projection
class XZarc: public Circle
{
  protected:
  double thmin, thmax; // theta min and max in radians for an arc

  public:

  XZarc(Vector3d origin, double r,double thmin, double thmax)
    : Circle(origin, r), thmin(thmin), thmax(thmax)
  {}

  virtual boost::optional<Vector3d> intersects(Lvec& shxyz);
  // keep distance_vector(Lvec) and distance(Vector3d) from superclass

  double get_thmin() { return thmin; }
  double get_thmax() { return thmax; }

};

// An actual surface object since it derives from Plane
// again we fix the infinite direction in e_y

// use this class to make sure vectors are projected to xz plane
class InfCylinder: public AnyPlane
{
  public:
  XZarc arc;
  bool solid; // if solid is False, return no overlaps or body intersections
  // we only construct 90 degree arcs and the whole arc should be in upper or lower plane
  bool is_lower_hemisphere; 
  bool is_upper_hemisphere; 

  InfCylinder(Vector3d origin, double r, double thmin, double thmax, bool solid)
    : arc(XZarc(origin.xzproject(), r, thmin, thmax)), solid(solid)
  {
    this->frame = Frame(origin);
    if (thmin <= 0 && thmax <= 0) {
      is_lower_hemisphere = true;
    }
    else {
      is_lower_hemisphere = false;
    }
    // 
    if (thmin >= 0 && thmax >= 0) {
      is_upper_hemisphere = true;
    }
    else {
      is_upper_hemisphere = false;
    }
  }

  Vector3d get_longaxis() { return frame.e2; }
  double get_r() { return arc.get_r(); }
  double get_thmin() { return arc.get_thmin(); }
  double get_thmax() { return arc.get_thmax(); }

  double distance(Vector3d pt) override 
  { return arc.distance(pt.xzproject()); }

  virtual boost::optional<Vector3d> intersects(Lvec& lv) override;

  virtual vector<overlap> overlap_vector(Capsule& body) override;
  virtual Eigen::RowVectorXf surface_grad(Capsule& body, contactf contact) override;
  Eigen::RowVectorXf surface_grad_part(Capsule& body, contactf contact, overlap over) override;

  double operator()(double x); 
  
  virtual std::string __str__();

  protected:
  bool _range_not_satisfied(double theta);
};


std::ostream& operator<<(std::ostream& os, InfCylinder cyl);

/////////////////////////////////////////////////////////////////////////////////

// steps of height h, separation sep
class InfSteps: public AnyPlane
{
  public:
  // parameters
  double _height, _sep;
  double _smallr; // smallr should be less than bacteria radius = 0.5
  double _close_n; 
  double h1, h2;

  // STL containers for holding shared pointers to surface objects
  std::map<int,std::shared_ptr<PartialPlaneX> > _xplanes;
  std::map<int,std::shared_ptr<PartialPlaneZ> > _zplanes;
  std::map<int,std::shared_ptr<InfCylinder> > _pcorners; // prev corners
  std::map<int,std::shared_ptr<InfCylinder> > _ncorners; // next corners
  boost::optional<std::vector<std::shared_ptr<Plane>>> _now_touching; // geometry objects we are touching after call to overlap_vector

  // Map is an appropriate data type to store objects that we can index by [...,-2,-1,0,1,2,..]

  // constructors
  InfSteps(Vector3d origin, Vector3d normal, double height, double sep, double smallr)
    : AnyPlane(origin, normal), _height(height), _sep(sep), _smallr(smallr)
  {
    _close_n = _smallr/_sep;
    // limits of vertical planes
    h1 = _smallr;
    h2 = _height - _smallr; 
  }
    
  InfSteps(double height, double sep, double smallr) 
    : InfSteps(Vector3d(), e_z, height, sep, smallr) 
  {}

  // work with a parameter n for which floor(n) is the integer representing this surface segment
  // _| |_| |_| ...
  // 0 1 2 3 4  ...
  float _get_n(double X); 
  double _get_X(float n); 
  double _get_X(int n); 
  bool is_lower(int n);
  bool is_close_prev(float n);
  bool is_close_next(float n);

  // This is composite surface
  // We want methods to construct all the nearby Planes and cylinders

  
  // geometry methods
  virtual boost::optional<Vector3d> intersects(Lvec& lv) override;

  virtual vector<overlap> overlap_vector(Capsule& body) override;
  virtual Eigen::RowVectorXf surface_grad(Capsule& body, contactf contact) override;
  int get_num_contacts(Capsule& body) override;

  // function which takes x coordinate and return step profile
  double operator()(double x); 


  // simple getters
  double get_sep() { return _sep; }
  double get_height() { return _height; }
  double get_smallr() { return _smallr; }

  //static const std::string __class__;
  std::string report_touching();
  std::string __str__() override;

  protected:
  //virtual double _surface_energy(ACell& cell, contactf energyf);

  // For these methods, if the desired object doesn't exist, create it
  std::shared_ptr<PartialPlaneX> _get_partial_plane(int n);
  std::shared_ptr<PartialPlaneZ> _get_partial_vertical(int n); 
  std::shared_ptr<InfCylinder> _get_prev_corner_at(int n);
  std::shared_ptr<InfCylinder> _get_next_corner_at(int n);
  vector<std::shared_ptr<Plane>> _get_forwards_geometry(int n);
  vector<std::shared_ptr<Plane>> _get_forwards_planes(int n);
  vector<std::shared_ptr<Plane>> _get_backwards_geometry(int n);


};

#endif
