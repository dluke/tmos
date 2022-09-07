
#ifndef __PILI_HPP__
#define __PILI_HPP__

#include <memory>
#include <cassert>
#include <iostream>
#include <numeric>
#include <unordered_set>
#include <iterator>

#include <exception>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>

#include <boost/optional.hpp>
#include "sample.hpp"
#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "potentials.hpp"
#include "lvec.hpp"
#include "chain.hpp"
#include "shape.hpp"
#include "sphere.hpp"
#include "wlc.hpp"


using std::abs;
using std::pow;
using std::vector;
using std::get;
using std::tuple;
using std::string;
using std::cout;
using std::endl;
using std::begin;
using std::end;
using std::unordered_set;
using std::asin;
using std::acos;


enum action {Extension, Retraction, Ret_on, Ret_off, 
  Ext_on, Ext_off, 
  Resample, Detach, Dissolve, Spawn};

struct fltq {
  Vector3d force;
  double torque;
};

//force, target
struct fpt {
  Vector3d force;
  Vector3d torque;
};

//torque in e1 e2 e3 axes
//dep
struct torque3d {
  Vector3d te1;
  Vector3d te2;
  Vector3d te3;
};



class Pili 
{
public:
  
  int idx;
  // leq is the intrinsic length that changes by increments
  // for free pili use leq, for bound pili pl is the anchor -> attachment length
  double pl, leq;
  double min_length;
  // Defined in cell body frame
  Vector3d axisEq, anchor;
  Vector3d lab_axisEq, lab_anchor;
  // attachement in lab frame
  // attachment -- point on surface with pili attached
  Vector3d attachment; 

  int cycles; // count the number of full extension -> retraction cycles of this machine
  bool new_cycle;
  int isbound, lastistaut; 
  double bound_time, free_time, lifetime, detach_time;
  bool is_fully_retracted;
  bool dissolve_on_retraction = false; 
  // update this parameter after extension,resampling,retraction and after MD step
  bool intersects_surface = false;
  // if detach can fail then track the last intersection length
  double last_inter_len;
  // if detach is forced to succeed then track the shortening length
  double last_shortening = 0;

  int ret_motor;
  int ext_motor;
  
  // break rules of polymorphism by defining a type string. This is only used for vtk output.
  virtual string get_typestr() { return "2d"; }

  /* 3d constructor */ 
  Pili(int idx, double pl, double leq, Vector3d anchor, Vector3d axisEq,
      int ret_motor, int ext_motor) 
    : idx(idx), pl(pl), leq(leq), anchor(anchor), axisEq(axisEq), 
    ret_motor(ret_motor), ext_motor(ext_motor)
  { 
    common_init();
    // assert(min_length > d_bound);
    assert(!(ret_motor && ext_motor));
  }
  
  void common_init();

  // convenience methods for teta (axis) and tetaEq (axisEq)
  double tetaEq() { return theta(axisEq); }

  // updates to pilus are independent of the surface and head position
  virtual void update(Vector3d anchor);

  bool istaut();

  double sym_force_l();
  double force_l();
  
  Vector3d get_attachment();

  double tau_eovers(double fmag);
  double get_rate_s(double fmag);
  double calc_k_dt(double fmag);
  double attach_rate();
  double get_anchor_angle();
  bool check_rcondition();
  double get_bound_retract_fraction();
  double get_bound_retract_rate();
  vector<double> get_rates();

  virtual double elongate();
  virtual double shrink();
  virtual double _shrink();
  virtual double shrink_by(double);
  virtual void switch_se();
  virtual bool has_motor();
  virtual void ret_motor_on();
  virtual void ret_motor_off();
  virtual void ext_motor_on();
  virtual void ext_motor_off();
    
  virtual void _attach_upkeep(double tme);
  virtual void detach(double tme);
  virtual void sample(void) {};
  virtual double get_k_resample(void) { return 0.; }

  // polymorphic copy 
  virtual std::shared_ptr<Pili> clone();

  virtual std::string __str__();

  // 
  //tuple<double, double> energies();
  double energy();
  
  // 3D specific methods
  // in principle a version of the 2d code could have a surface
  // otherwise this is dummy method to support passing arbitrary surfaces
  // to arbitrary subclasses of pili
  virtual boost::optional<Vector3d> check_attach(Vector3d lab_anchor, Capsule body, Shape* surface)
  { throw std::runtime_error("not implemented"); }
  virtual double shorten_to_detach(Capsule body, Shape* surface)
  { throw std::runtime_error("not implemented"); }
  virtual int attach(double tme, Capsule body, Shape* surface)
  { throw std::runtime_error("not implemented"); }
  
  // generic paramters describing a pilus
  //double ks;
  static double ks;
  static double force_min_length;
  static double max_length;
  static double inside_length;
  static double d_free;
  static double d_bound;
  static double k_adh;
  static double adh_accept;
  static double f_stall;
  static double kf_sh;
  static double kf_ex;
  static double kb_sh;
  static double kb_ex;
  static double anchor_angle_threshold;
  static double anchor_angle_smoothing_fraction;
  static double f_release;
  static double dwell_time;
  static double tau_e;
  static double free_tau_e;
  static double force_threshold;

  static double k_ext_off;
  static double k_ret_off;
  static double k_ext_on;
  static double k_ret_on;
  static double eff_k_ret_on;
  static double k_dissolve;
  static double k_resample;
  static double detach_grace_length; 
  static double detach_grace_time; 

  static double sh_to_elong_rate; 
  static double free_tau_eovers;
  static double low_force_proportion;
  static double high_force_proportion;

  static bool allow_bound_extension;
  static bool allow_bound_ext_motor;
  static bool force_bound_retraction;
  static bool force_detachment;
  static bool enforce_stalling;
  // these two are dep
  static bool simplify_free_pili_etor;
  static int rcondition;


  static bool force_max_length;

};

// base class for 3d Pili. Implements an attach method which always returns no attachment
// useful for testing pili dynamics in the absense of a surface
class Pili3d: public Pili
{
  public:
  using Pili::Pili; // inherit a set of constructors

  double get_last_inter_len() { return last_inter_len; }

  void update(Vector3d lab_anchor) override;

};

class PiliRod3d: public Pili3d
{
  public:

  string get_typestr() override { return "rod"; }

  PiliRod3d(int idx, double pl, double leq, Vector3d anchor, Vector3d axisEq, 
      int ret_motor, int ext_motor) 
    : Pili3d(idx, pl, leq, anchor, axisEq, ret_motor, ext_motor) {
  }

  boost::optional<Vector3d> check_attach(Vector3d lab_anchor, Capsule body, Shape* surface);
  int attach(double tme, Capsule body, Shape* surface);

  std::shared_ptr<Pili> clone() override;

  std::string __str__() override;

};

class PiliWLC3d: public Pili3d
{
  public:
  string get_typestr() override { return "wlc"; }

  double ka, a, lasta;
  int n;

  std::shared_ptr<WLCgeneratorBase> wlc;

  // Stores reference to a row and not a copy (how to check?)
  std::unique_ptr<Eigen::VectorXd> chain_theta;
  std::unique_ptr<Eigen::VectorXd> chain_phi;
  // todo use shared ptr
  std::unique_ptr<Chain> chain_instance;

  // Manage a copy of the last chain instance
  std::shared_ptr<Chain> new_chain_instance(void);
  std::shared_ptr<Chain> transformed_chain_instance(Capsule body);
  Chain& get_chain_instance() { return *chain_instance; }
  void set_chain_instance(std::unique_ptr<Chain> ch) { chain_instance.reset(ch.release()); }


  PiliWLC3d(int idx, double pl, double leq, Vector3d anchor, Vector3d axisEq, int ret_motor, int ext_motor,
      std::shared_ptr<WLCgeneratorBase> wlc)
    : Pili3d(idx, pl, leq, anchor, axisEq, ret_motor, ext_motor), wlc(wlc)
  {
    this->ka = wlc->get_ka();
    this->a = wlc->get_a();
    this->n = (int)(leq/a) + 1;
    // track the length of the last element
    lasta = leq - (n-1)*a;
    // initialise chain
    this->sample();
  }

  PiliWLC3d(int idx, double pl, double leq, Vector3d anchor, Vector3d axisEq, int ret_motor, int ext_motor) 
    : PiliWLC3d(idx, pl, leq, anchor, axisEq, ret_motor, ext_motor, std::make_shared<KPgenerator>())
  {}


  boost::optional<Vector3d> check_attach(Vector3d lab_anchor, Capsule body, Shape* surface);
  double shorten_to_detach(Capsule body, Shape* surface);
  int attach(double tme, Capsule body, Shape* surface);
  int _attach(double tme, Vector3d lab_anchor, Vector3d target);

  double get_lasta() { return lasta; }

  // new elongate and shrink methods to update n segments and last segment length
  double elongate() override;
  double shrink() override;
  double shrink_by(double) override;

  void sample(void) override;
  // double get_k_resample(void) override { return (isbound) ? 0. : k_resample; }
  
  // resample even bound pili
  double get_k_resample(void) override { return k_resample; }
  // now during detachment we need to update n, lasta
  void detach(double tme);

  std::string __str__() override;

  friend void swap(PiliWLC3d& lhs, PiliWLC3d& rhs);
  PiliWLC3d(const PiliWLC3d& other);
  
  PiliWLC3d& operator=( PiliWLC3d tmp);
  std::shared_ptr<Pili> clone() override;

};

#endif
