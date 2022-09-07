#include "potentials.hpp"


double smooth(double r, double a, double b) {
  if (b <= a) {
    return 0.;
  }
  if (r < a) {
    return 1.;
  }
  else if (r > b) {
    return 0.;
  }
  else {
    double t = (r - a)/(b - a);
    return 1 + t*t * (2*t - 3);
  }
} 

double smooth_deriv(double r, double a, double b) {
  if (r < a) {
    return 0;
  }
  else if (r > b) {
    return 0;
  }
  else {
    double t = (r - a)/(b - a);
    return 6*t*(t-1);
  }
}
 
//////////////////////////////////////////////////////////////////////////////////
//

double wca_force_abs(double eps, double rmin, double r) {
  //std::cout << "get wca force with eps = " << eps << std::endl;
  if (r >= rmin) {
    return 0;
  }
  else {
    double rm6 = pow(rmin/r, 6);
    return eps * (12/r) * (rm6*rm6 - rm6);
  }
}

double wca_energy_abs(double eps, double rmin, double r)
{
  if (r >= rmin) {
    return 0;
  }
  else {
    double rm6 = pow(rmin/r, 6);
    return eps * (rm6*rm6 - 2*rm6) + eps;
  }
}

double Spotential::lj_energy(double r) {
  double rm6 = pow(re/r, 6);
  return eps * (rm6*rm6 - 2*rm6);
}

double Spotential::lj_grad(double r) {
  double rm6 = pow(re/r, 6);
  return -eps * (12/r) * (rm6*rm6 - rm6);
}



// TODO repair smoothing function based on /home/dan/usb_twitching/pili/src/test_python/test_potentials.py
double Spotential::smoothing_function(double r) 
{
  return smooth(r, a, b);
}

double Spotential::smoothing_deriv(double r)
{
  return smooth_deriv(r, a, b);
}

