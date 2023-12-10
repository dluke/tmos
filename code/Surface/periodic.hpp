
#ifndef __PERIODIC_HPP__
#define __PERIODIC_HPP__

#include <cmath>
#include <vector>
#include <utility>
#include <tuple>
#include <string>
#include <exception>
#include <boost/optional.hpp>


// NLopt library
#include <nlopt.hpp>

#include "vector3d.hpp"
#include "lvec.hpp"
#include "plane.hpp"
#include "sphere.hpp"


using namespace Base;

//using namespace nlopt;
//using namespace Shape;

using std::abs;
using std::sin;
using std::cos;
using std::sqrt;
using std::vector;

const double pi = M_PI;
const double pih = M_PI/2.;

 
class SinePlane;

struct clcl_data{
  double c1, c2, l1, l2;
  SinePlane* sp;
  int itercount;
};

struct f_capleq {
  clcl_data data;
  int itercount;

  double operator()(const std::vector<double> &mx, std::vector<double> &grad);

  static double wrap(
      const std::vector<double> &x, std::vector<double> &grad, void *data) 
  { return (*reinterpret_cast<f_capleq*>(data))(x, grad); }

};


class SinePlane: public AnyPlane
{

  public:
  // the function for the surface is Asin(x/B)
  double A,B,invB;
  double xperiod;
  int iperiod, iconst, inorm;
  // let frame->e1 be the periodic dimension
  // let frame->e2 be the constant direction
  // let frame->e3 be the normal direction
  // and let these e1,e2,e3 be cardinal directions (simplifing projections)

  int nsign; // the sign of the surface normal

  f_capleq capleq;
 
  // constructors

  SinePlane(Frame frame, double A, double B, int nsign);
  SinePlane(Frame frame, double A, double B) : SinePlane(frame, A, B, 1) {};
  SinePlane(double A, double B) : SinePlane(Frame(), A, B) {}
  SinePlane() : SinePlane(1., 1.) {}

  ///////////////////////////////////////////////////////////////////////////////

  // mathematical form for the curve
  double form(double t) { return A * sin(invB*t); }
  double dform(double t) { return A * invB * cos(invB*t); }
  
  // the form but in 3d
  Vector3d sp_form(double t);
  
  virtual bool contains(Vector3d pt);

  // normal vector at pt
  Vector3d normal_form(double t);
  Vector3d normal(Vector3d pt);
  
  // geometric functions
  using AnyPlane::intersects;
  double closest_t(Vector3d pt);
  Vector3d distance_vector(Vector3d pt);
  boost::optional<Vector3d> intersects(Lvec& lv);
  vector<overlap> overlap_vector(Capsule& body);
  
  // debugging 
  std::string report_touching();

  private:
  double _closest_t(double, double);
  std::vector<Vector3d> _now_contact;
  std::vector<double> _contact_m;
  
};


class PartialSinePlane: public SinePlane
{
};

#endif
                 
