


#ifndef __WLC_HPP__
#define __WLC_HPP__

#include <memory>
#include <vector>
#include <Eigen/Dense>
#include "sample.hpp"
#include "vector3d.hpp"
#include "frame.hpp"
#include "chain.hpp"

// todo - remove this dependency on pybind11
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

//typedef Eigen::Ref<MatrixType, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> EigenDRef;
//typedef py::EigenDRef EigenDRef;

using std::vector;


// if we use this as the type for PiliWLC3d constructor,  does it need to be constructible?
class WLCgeneratorBase 
{
  public:
  WLCgeneratorBase(double ka, double a, int nmax) : 
    ka(ka), a(a), nmax(nmax) {} 

  // virtual std::unique_ptr<Chain> sample(int n, Vector3d axis) = 0;

  // sample azimuthal angles
  virtual std::unique_ptr<Eigen::VectorXd> sample_phi(void);
  // sample polar angles
  virtual std::unique_ptr<Eigen::VectorXd> sample_theta(void) {};

  
  double get_ka() { return ka; }
  double get_a() { return a; } 
  
  protected:
  int nmax;
  double ka, a;
};

class KPgenerator: public WLCgeneratorBase
{
  public:
  
  KPgenerator() : WLCgeneratorBase(0.0,1.0,10), T(303) {}
  KPgenerator(double ka, double a, int nmax, double T) : WLCgeneratorBase(ka,a,nmax), T(T) {
    beta = 1./(kb*T);
    Z = (exp(beta*ka/a) - 1)/beta;
    Lp = beta*ka;
  }

  double one_theta(void);
  std::unique_ptr<Eigen::VectorXd> sample_theta(void) override;

  protected:
  const double kb = 1.38064852 * (10e-23 * 10e18);
  // temperature
  double T;
  double beta;
  // partition functiion, not including a constant factor of 2*\pi*a^3/k_a
  double Z; 
  // persistance length
  double Lp; 
};

// only constructed by python 
class WLCgenerator: public WLCgeneratorBase
{
  public:

  vector<py::EigenDRef<Eigen::MatrixXd> > vcdata;
  //using WLCgeneratorBase::WLCgeneratorBase;
  WLCgenerator(double ka, double a, int nmax) : 
    WLCgeneratorBase(ka, a, nmax) 
  {
    // nmax = nsegments - 1
    vcdata.reserve(nmax);
    //vcdata.resize(nmax);
  }

  void append(py::EigenDRef<Eigen::MatrixXd> cdata);
  void setnth(int n, py::EigenDRef<Eigen::MatrixXd> cdata);
  py::EigenDRef<Eigen::MatrixXd> getnth(int n);

  // sample polar angles
  std::unique_ptr<Eigen::VectorXd> sample_theta(void);

  // build up a chain along the specified axis
  std::unique_ptr<Chain> sample(int n, Vector3d axis);

};

#endif
