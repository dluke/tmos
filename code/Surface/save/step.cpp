
#include "step.hpp"
// newton method
#include <boost/format.hpp>
#include <boost/math/tools/roots.hpp>
#include <exception>
#include <iostream>

using std::cout;
using std::endl;

const int Step::STEP_LOWER = -1;
const int Step::ONSTEP = 0;
const int Step::STEP_UPPER = 1;


// intersection and distance function generators
Fresult step_intersection_generator(double h, double B, 
    double rx, double ry, double lx, double ly, double t) {
  double expexp = exp(-B*(rx + t*lx));
  double F = (ry + t*ly) * expexp + t*ly - h + ry;
  double dFdt = ly * expexp - (ry + t*ly)*B*lx*expexp + ly;
  return Fresult(F, dFdt);
}

double step_distance_generator(double h, double B,
    double rx, double ry, double x) {
    double hexp = h/(1 + exp(-B*x));
    return (x-rx)*(x-rx) + (hexp-ry)*(hexp-ry);
}


// define distance to point function

//vector<double> initial_guess(double rx, double ry) {
  //return 0.;
//}

int Step::segment_for_pt(double x) {
  double xl, xh;
  std::tie(xl, xh) = get_smoothbounds();
  // shift the partitioning by the assumed cell radius
  if (x < xl)
    return STEP_LOWER;
  else if (x > xh)
    return STEP_UPPER;
  else
    return ONSTEP;
}

// return the x coordinate along the step for the shortest distance pt
// assume that pt.z > Step(pt.x), and that segment_for_pt function
// is sufficient in it's task
Vector3d Step::distance_vector(Vector3d pt) {
  // tmp check called
  std::cout << "Calling Step::distance_vector" << std::endl;
  // tmp use Brent's algorithm
  // unless we are on the upper or lower planes
  double dminx;
  int part = segment_for_pt(pt.x);
  if (part == STEP_LOWER) {
    std::cout << "on lower step" << std::endl;
    return Vector3d(0, 0, pt.z - this->lower.get_origin().z);
  }
  else if (part == STEP_UPPER) {
    std::cout << "on upper step" << std::endl;
    return Vector3d(0, 0, pt.z - this->upper.get_origin().z);
  }
  else
  {
    boost::uintmax_t it;
    it = this->maxit;

    std::cout << "on step" << std::endl;
    // construct the function for this point
    boost::function<double(double)> D2f = boost::bind(D2f_rx_rz, pt.x, pt.z, _1);
    // use brents algorithm to get x 
    std::pair<double,double> bounds = get_Dbounds(0.5); // tmp
    std::pair<double, double> sx = 
      boost::math::tools::brent_find_minima(D2f, bounds.first, bounds.second, bits, it);
    cout << "brent method bounds are " << bounds.first << " " << bounds.second << endl;
    cout << "x pt on step is " << sx.first << endl;
    cout << "number of iterations " << it << endl;
    // sx.second = the actual distance which is useful to return for efficiency
    // rough check for convergence here
    if (it >= 50) { 
      throw std::runtime_error("Brent's method failed to converge after 50 iterations");
    }
    // Use x to get the vector pt - (x, Step(x))
    double sz = this->sigmoid(sx.first);
    return Vector3d(pt.x-sx.first, 0, pt.z-sz);
  }
}


boost::optional<Vector3d> Step::intersects(Lvec lv) 
{
  // check this method is called
  std::cout << "calling Lvec->Step intersection method" << endl;
  //
  // bind the intersection function for this line segment
  boost::function<std::pair<double, double> (double) >  Fwgrad = boost::bind(
      this->Finter, lv.source.x, lv.source.z, lv.line.x, lv.line.z, _1);

  // where to define max iterations and precision?
  boost::uintmax_t it;
  int bits = std::numeric_limits<double>::digits;
  // half precision
  bits /= 2;
  // on sigmoid check precision
  double float_tol = 0.0000001;

  double startt = 0.;
  double endt = 1.;
  double tinter1, tinter2;
  it = 50;
  tinter1 = boost::math::tools::newton_raphson_iterate(
      Fwgrad, startt, 0., 1., bits, it);

  Vector3d pt1 = lv.get_pt(tinter1);
  // need to check if we really hit the surface and where
  double sz1 = this->sigmoid(pt1.x);
  if ( abs( pt1.z  - sz1 ) < float_tol ) {
    // for now lets accep this one
    return pt1;
  }
  
  // otherwise we can try the other end
  it = 50;
  tinter2 = boost::math::tools::newton_raphson_iterate(
      Fwgrad, startt, 0., 1., bits, it);
  Vector3d pt2 = lv.get_pt(tinter2);
  // need to check if we really hit the surface and where
  double sz2 = this->sigmoid(pt2.x);
  if ( abs( pt2.z  - sz2 ) < float_tol ) {
    // for now lets accept this one aswell if we find it
    return pt2;
  }

  // otherwise no intersection
  // current max steps used for no intersection : 2 * 50
  return boost::none;


}


std::string Step::__str__() {
  std::string form = "Step: h(%f) B(%f)\norigin%s\nnormal%s\nLower %s\nUpper %s";
  return ( boost::format(form) 
      % this->h
      % this->B
      % frame.origin.__str__() 
      % frame.e3.__str__() 
      % this->lower.__str__()
      % this->upper.__str__()
      ).str();
}

