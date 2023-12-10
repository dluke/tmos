

#include "sample.hpp"

using std::acos;
using std::sinh;
using std::log;
using std::exp;

const double pi = M_PI;

int init_generator(long pyseed) {
  mtgen.seed(pyseed);
  return 0;
}


double r2()
{
  static std::uniform_real_distribution<double> dis(0.0, 1.0);
  return dis(mtgen);
}

//  Sample Von-Mises Fisher Distribution

Vector3d uniform_v() 
{
  double theta = 2 * M_PI * r2();
  return Vector3d(cos(theta), sin(theta), 0);
}

double rW(double kappa) 
{
  //# inverse cumulative distrubution function 
  double c = 2*sinh(kappa)/kappa;
  return 1/kappa * log(exp(-kappa) + kappa * c * r2());
}

Vector3d sample_vmf(double kappa) 
{
  Vector3d v = uniform_v();
  double W = rW(kappa);
  double wf = sqrt(1-W*W);
  return Vector3d(wf*v.x, wf*v.y, W);
}

// use L use that R = 1.
double vmf_stretch_function(double L, double x)
{
  return L/M_PI * x + (M_PI-L)/2.;
      //return L/np.pi  * x + (np.pi-L)/2.
}

// introduce the cell radius R at this point so we don't need a complicated shrinking operation
// L is the body length, i.e. 2.
Vector3d modified_vmf(double kappa, double L, double R)
{
  Vector3d unitX = sample_vmf(kappa);
  Vector3d X = R * unitX;
  double theta = acos(unitX.dot(e_z));
  if (theta > M_PI/2.) {
    double md = vmf_stretch_function(2*L, theta) - M_PI/2.;
    X.z = 0;
    X = R * X.unit(); // X now on cap and perpendicular to axis
    X.z = -R * md; // R appears again to modify distance md 
  }
  X.z += L/2.; // shift so body center is at 0
  return X; 
}

// for chain generation

Eigen::RowVectorXd sign_row(int size) 
{
  Eigen::RowVectorXd sgn(size);
  for (int i = 0; i < size; i++)
  {
    if (r2() > 0.5)
      sgn[i] = -1;
    else
      sgn[i] = 1;
  }
  return sgn;
}

Eigen::RowVectorXd phi_row(int size) 
{
  static std::uniform_real_distribution<double> rphi(-pi, pi);
  Eigen::RowVectorXd phi(size);
  for (int i = 0; i < size; i++)
  {
    phi[i] = rphi(mtgen);
  }
  return phi;
}


//::uniform_int_distribution<int> init_uniform_rint(int min, int max) 

std::uniform_int_distribution<int> init_uniform_rint(int min, int max) 
{
  // .?
  std::uniform_int_distribution<int> uni;
    *new std::uniform_int_distribution<int>(min,max);
  return uni;
}

