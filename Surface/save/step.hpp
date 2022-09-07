
#ifndef __STEP_HPP__
#define __STEP_HPP__

#include <memory>
#include <utility>
#include <tuple>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include "vector3d.hpp"
#include "plane.hpp"

#include <boost/math/tools/minima.hpp>

using std::exp;
using std::vector;

// for returning value and derivative
typedef std::pair<double, double> Fresult;

// define the intersection function for line and step
Fresult step_intersection_generator(double h, double B, 
    double rx, double ry, double lx, double ly, double t);

// generator for distance-to-pt function
//
double step_distance_generator(double h, double B,
    double rx, double ry, double x);

// typedef for pt distance partial function
typedef boost::function<double(double, double, double)> D2function;
// typedef for line intersection partial function
typedef boost::function<Fresult(double, double, double, double, double)> bflinter;


class Step: public Plane
{
  public:
  // Surface with Step, position identifiers
  //
  static const int STEP_LOWER;
  static const int ONSTEP;
  static const int STEP_UPPER;

  bflinter Finter;
  D2function D2f_rx_rz;
  Plane lower, upper;
  // tmp minimiser setup
  int bits;
  int maxit;


  Step(double h, double B) : h(h), B(B) 
  {
    // define the Frame
    //this->frame = Frame(Vector3d(), e_x, e_y, e_z);
    this->frame = Frame();
    this->Finter = boost::bind(step_intersection_generator, h, B,
        _1, _2, _3, _4, _5);
    this->D2f_rx_rz = boost::bind(step_distance_generator, h, B, _1, _2, _3);
    // construct planes on both sides
    lower = Plane(Vector3d(), e_z);
    upper = Plane(Vector3d(0,0,h), e_z);
    setup_brents();
  }

  void setup_brents() {
    //https://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/roots/brent_minima.html
    // setup belongs to the algorithm but we can dump it here for now
    this->bits = std::numeric_limits<double>::digits;
    this->bits /= 2;
    this->maxit = 50;
  }

  std::pair<double, double> get_smoothbounds()
  {
    // tmp
    return std::make_pair(-3, 3);
  }

  std::pair<double, double> get_Dbounds(double r) 
  {
    std::pair<double, double> sbounds = get_smoothbounds();
    return std::make_pair(sbounds.first - r, sbounds.second + r);
  }

  // define wether to use plane or step for distance calculation 
  // for distance to a pt
  int segment_for_pt(double x);

  // define the step function
  double sigmoid(double x) const
  {
    return h/(1 + exp(-B*x));
  }
  double operator()(double x) const { this->sigmoid(x); }

  // compute the distance to a point is good enough for circle curve intersection
  //double distance_ptx(Vector3d pt);
  Vector3d distance_vector(Vector3d pt) override;


  boost::optional<Vector3d> intersects(Lvec lv) override;

  //bool contains(Vector3d pt);
  
  Vector3d get_step_direction() { return frame.e1; }
  Vector3d get_const_direction() { return frame.e2; }

  std::string __str__() override;


  private:
  double h, B; // h height, B steepness
};



#endif
