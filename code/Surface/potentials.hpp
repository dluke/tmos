

#ifndef __POTENTIALS_HPP__
#define __POTENTIALS_HPP__ 

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/math/tools/polynomial.hpp>

#include <iostream>
#include <cmath>
#include <vector>

using std::cout;
using std::endl;
using std::pow;
using std::abs;
using std::vector;


using namespace boost::math::tools; // for polynomial


// potentials could be part of a surface class which wraps Shape classes but 
// for the moment just write them here
//

//http://www.sklogwiki.org/SklogWiki/index.php/Pseudo_hard_sphere_potential
//
//Pseudo hard sphere WCA potential
//


double wca_force_abs(double eps, double rmin, double r);
double wca_energy_abs(double eps, double rmin, double r);

////////////////////////////////////////////////////////////////////////////////

polynomial<double> construct_veplj();
polynomial<double> construct_veps_();

polynomial<double> deriv(polynomial<double> p);

// class can handle constructing surface potential functions from polynomials
class Spotential
{
  public:
  double eps; 
  double re; // the position of minimum of the potential well should be ~ 0.490
  double R; // Potential goes to zero at R
  double a, b; // convenience 
  bool repulsive;

  Spotential(double eps, double re, double R)
    : eps(eps), re(re), R(R) 
  {
    a = re;
    b = R;
  }

  // repulsive only functions
  // no dependence on re or eps_attract
  double wca_energy(double r)
  {
    return wca_energy_abs(this->eps, this->R, r);
  }

  double wca_force(double r)
  {
    return wca_force_abs(this->eps, this->R, r);
  }

  // repulsive/attractive functions
  double sforce(double r)
  {
    if (r < a) {
      return -lj_grad(r);
    }
    else if (r > b) {
      return 0.;
    }
    else {
      // return force
      return -1 * ( lj_grad(r) * smoothing_function(r) + lj_energy(r) * smoothing_deriv(r) );
    }
  }

  double senergy(double r)
  {
    return lj_energy(r) * smoothing_function(r);
  }

  double smoothing_function(double r);
  double smoothing_deriv(double r); 
  double lj_energy(double r);
  double lj_grad(double r);
};

double smooth(double r, double a, double b);

#endif
