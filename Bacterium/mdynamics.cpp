
#include "mdynamics.hpp"

#include <exception>

#include <iostream>
using std::cout;
using std::endl;


pair<int, double> MDintegrate::step(Cell3d& cell)
{
  // given a cell setup the integrator for this run
  double rms;
  rms = this->reset(cell);
  int i;
  for (i = 0; i < maxsteps; i++) // don't need var i
  {
    if (rms <= target) { 
      isprepped = 0; // ensure call to reset() for next time
      break;
    }
    rms = this->one_step(cell); // Equilibrate3d object should call on overridden method
  }
  // currently must enforce force < target condition in python loop
  isprepped = 0;
  summary.nsteps = i;
  summary.rms = rms;
  summary.accepted_state = Capsule(cell.get_body());
  summary.accepted_pilus = cell.last_pilus_ref->clone();
  return make_pair(i, rms);
}


/////////////////////////////////////////////////////////////////////////////////

double Equilibrate3d::reset(Cell3d& cell) 
{
  // reset method
  x_ = cell.get_state(); // get_state is only called once to initialise minimiser
  // clear summary object
  summary = Summary();
  summary.prior_state = Capsule(cell.get_body());
  summary.prior_pilus = cell.last_pilus_ref->clone();

  cell.get_body().set_Rk(Rk(Vector3d(x_[3], x_[4], x_[5])));
  g_ = cell.grad(); // initialise function gradient
  rms = g_.norm()/sqrt(6);
  if (rms <= target) {
    return rms; // exit early
  }
  _dx = Eigen::RowVectorXf::Zero(6);
  _v = Eigen::RowVectorXf::Zero(6);
  ifnorm = 1./g_.norm();
  vnorm = _v.norm();
  _reset_params();
  return rms;
}

void Equilibrate3d::_reset_params() {
  dt = dt_start;
  st = 0; // temporary step counter
  alpha = alpha_start;
  isprepped = 1;
}

//////////////////////////////////////////////////////////////////////////
// 3d implementation of FIRE
// copy from pele but adapt to my cell3d object

void Equilibrate3d::_ForwardEuler_integration(Cell3d& cell) {

  _v -= dt * g_;     //update velocity, assumes all masses are 1
  _dx = dt * _v;     //build displacement vector
  
  // if set max displacement do it here
  // --
  double normdx = _dx.norm();

  if (normdx > _maxstep) {
      _dx *= (_maxstep / normdx); //resize displacement vector if greater than _maxstep
  }

  x_ += _dx;

  cell.set_state(x_); 
  // Must set pstate in addition to set_state! 
  cell.get_body().set_Rk(Rk(Vector3d(x_[3], x_[4], x_[5])));
  // currently relies on body.get_state()
  g_ = cell.grad(); 
  rms = g_.norm()/sqrt(6); // updated rms
}


double Equilibrate3d::one_step(Cell3d& cell) { 
  if (!isprepped) {
    throw std::runtime_error("Must initialise integrator with a cell object by calling reset(cell)");
  }

  st += 1;

  /*equation written in this conditional statement _v = (1- _a)*_v + _a * funit * vnorm*/

  _v = (1. - alpha) * _v - alpha * g_ * ifnorm * vnorm;

  /*run MD*/
  this->_ForwardEuler_integration(cell);

  double P = -1 * _v.dot(g_);

  if (P >= 0) {
    if (st > Nmin) {
      dt = std::min(dt * f_inc, dt_max);
      alpha *= f_alpha;
    }
    ifnorm = 1./g_.norm();
    vnorm = _v.norm();
  }
  else {
    dt *= f_dec;
    alpha = alpha_start;
    st = 0;
    _v = Eigen::RowVectorXf::Zero(6);
    // pele implements stepback here
  }

  return rms;
}
