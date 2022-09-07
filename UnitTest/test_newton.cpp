
#include <iostream>
//#include "vector3d.hpp"
      
#include <utility>
#include <functional>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <functional>
#include <cmath>
   
using std::cout;
using std::endl;
using std::exp;


// doesn't work
// h, B constants, 
//std::function<double(double)> make_sigmoid(double h, double B)
//{
  //return [](double x) { return h/(1 + exp(-B*x)); };
//}

double h = 1;
double B = 4;
//double rx = -0.1; double ry = 0.5;
double rx = -3.4; double ry = 0.1;


// define function of one variable
struct D2f
{
  double operator()(double const& x)
  {
    double hexp = h/(1 + exp(-B*x));
    return (x-rx)*(x-rx) + (hexp-ry)*(hexp-ry);
  }
};


double prx = -0.1;
double pry = 0.5;
//double lx = 5; // hit
double lx = 0.;  // miss
double ly = 0.2;

typedef std::pair<double, double> Fresult;

Fresult step_intersection(double t)
{
  double expexp = exp(-B*(prx + t*lx));
  double F = (pry + t*ly) * expexp + t*ly - h + pry;
  double dFdt = ly * expexp - (pry + t*ly)*B*lx*expexp + ly;
  return Fresult(F, dFdt);
}

double bstepinter(double t)
{
  double expexp = exp(-B*(prx + t*lx));
  return (pry + t*ly) * expexp + t*ly - h + pry;
}


Fresult getpt(double prx, double pry, double lx, double ly, double t) { 
  return Fresult(prx + t * lx, pry + t * ly);
}

double sigmoid(double h, double B, double x)
{
  return h/(1 + exp(-B*x));
}


const boost::uintmax_t maxit = 50;

int main() {
  // construct sigmoid
  //
  // try brents algorithm
  int bits = std::numeric_limits<double>::digits;
  // half precision
  bits /= 2;
  boost::uintmax_t it = maxit;

  double limnx = -3.5; double limmx = 3.5;
  std::pair<double, double> r = boost::math::tools::brent_find_minima(D2f(), limnx, limmx, bits, it);

  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout << "x at minimum = " << r.first << ", f(" << r.first << ") = " << std::sqrt(r.second) << std::endl;
  std::cout << " after " << it << " iterations. " << std::endl;
  std::cout << std::endl;

  // ok now test line vector intersection
  double tinter;
  //std::pair<double, double> tinter;
  double delta = 0.0;
  double tlim1 = 0-delta; double tlim2 = 1+delta;
  Fresult pt;

  it = 50;
  tinter = boost::math::tools::newton_raphson_iterate(
      step_intersection, 0., tlim1, tlim2, bits, it);
  std::cout << "NR method t = " << tinter << " for guess t0 = 0" << endl;
  std::cout << " after " << it << " iterations. " << std::endl;
  pt = getpt(prx, pry, lx, ly, tinter);
  std::cout << "pt " << pt.first << " " << pt.second << std::endl;
  std::cout << "sigmoid " << pt.first << " " << sigmoid(h, B, pt.first) << std::endl;
  std::cout << std::endl;

  it = 50;
  tinter = boost::math::tools::newton_raphson_iterate(
      step_intersection, 1., tlim1, tlim2, bits, it);
  std::cout << "NR method t = " << tinter << " for guess t0 = 1" << endl;
  std::cout << " after " << it << " iterations. " << std::endl;
  pt = getpt(prx, pry, lx, ly, tinter);
  std::cout << "pt " << pt.first << " " << pt.second << std::endl;
  std::cout << "sigmoid " << pt.first << " " << sigmoid(h, B, pt.first) << std::endl;
  std::cout << std::endl;

  // Tryout TOMS algorithm
  //https://www.boost.org/doc/libs/1_60_0/libs/math/doc/html/math_toolkit/roots/roots_noderiv/TOMS748.html


}
