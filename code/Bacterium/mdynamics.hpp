
#ifndef __MDYNAMCIS_HPP__
#define __MDYNAMCIS_HPP__

#include <array>
#include <utility>
#include <boost/optional.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <Eigen/Dense>

// use Eigen for minimisation functions
//#include <Eigen/Dense>

#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "sphere.hpp"
#include "pili.hpp"
#include "cell.hpp"
#include "cell3d.hpp"

using std::pair;
using std::sqrt;
using std::make_pair;

struct Summary {
  int nsteps;
  float rms;
  Capsule prior_state;
  Capsule accepted_state;
  std::shared_ptr<Pili> prior_pilus;
  std::shared_ptr<Pili> accepted_pilus;
};


// Is Energy Minimisation not Molecular Dynamics
class MDintegrate
{
  protected: MDintegrate(int maxsteps, double target) 
    : maxsteps(maxsteps), target(target), rms(99999999) {}
  

  public:
  Eigen::RowVectorXf x_, g_; // position and gradient vectors

  // step until maxsteps is reached or we hit target force
  virtual pair<int, double> step(Cell3d& cell);

  // setup minimiser for a run
  virtual double reset(Cell3d& cell) = 0;

  // return the magnitude of the force
  virtual double one_step(Cell3d& cell) = 0;

  Eigen::RowVectorXf get_x() { return x_; }
  double get_rms() { return rms; }
  Summary get_summary() { return summary; }
  
  // currently not used
  bool check_target() { this->rms = g_.norm()/sqrt(6); return (rms <= target); }

  int maxsteps;
  double target, rms;
  bool isprepped;
  Summary summary;
};

// Fire method

class Equilibrate3d: public MDintegrate
{
  public:
  Equilibrate3d(double dt_start, int maxsteps, double target, int Nmin, 
      double f_inc, double f_dec, double alpha_start, double f_alpha, double dt_max, double maxstep)
    : MDintegrate(maxsteps, target), dt_start(dt_start), Nmin(Nmin),
      f_inc(f_inc), f_dec(f_dec), alpha_start(alpha_start), f_alpha(f_alpha), dt_max(dt_max), _maxstep(maxstep)
  {
    _reset_params();
    // must call reset with a cell object before step()
  }

  public:
  void _reset_params();
    
  // call reset to initialise the optimser before a run
  double reset(Cell3d& cell) override;
  double one_step(Cell3d& cell) override;
  // helper method
  void _ForwardEuler_integration(Cell3d& cell);


  Eigen::RowVectorXf _dx, _v; // delta position and velocity 

  double dt, dt_start, ifnorm, vnorm;
  int Nmin, st; // st is temporary step counter 
  double f_inc, f_dec, alpha_start, f_alpha, dt_max, alpha;
  double _maxstep; // the maximum step size // should be << 0.005
  

  // for debugging 
  Eigen::RowVectorXf get_dx(void) { return _dx; }

};



#endif
