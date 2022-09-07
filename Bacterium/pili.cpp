
#include "pili.hpp"

#include <boost/format.hpp>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Core>
// #include <Eigen/src/core/ArithmeticSequence>

using std::swap;


// The values of these static variables must be initialised or the python modules
// fail to import. However python is responsible for setting all parameters.

//defaults, but in reality defaults are defined in parameters.xml (check?)
double Pili::ks = 0.;
double Pili::force_min_length = 0.; 
double Pili::force_threshold= 0.; 
//
double Pili::max_length = 0.; 
double Pili::inside_length = 0.; 
double Pili::d_free  = 0.01;
double Pili::d_bound = 0.005;
double Pili::k_adh = 5.;
double Pili::adh_accept = 1.;
double Pili::f_stall = 100.;
double Pili::f_release= 95.;
double Pili::kf_sh = 0.8/Pili::d_free;
double Pili::kf_ex = 0.8/Pili::d_free;
double Pili::kb_sh = 0.7/Pili::d_bound;
double Pili::kb_ex = 0.7/Pili::d_bound;
double Pili::anchor_angle_threshold = M_PI/2;
double Pili::anchor_angle_smoothing_fraction = 0.0;

// Koch motor binding/unbinding rates in units of s^-1
double Pili::k_ext_off = 1.0/1.6;
double Pili::k_ret_off = 1.0/9.1;
double Pili::k_ext_on = 1.0/2.4;
double Pili::k_ret_on = 1.0/0.4;
double Pili::eff_k_ret_on = -1.0;
double Pili::k_dissolve = 0.;
double Pili::k_resample = 1.;
double Pili::detach_grace_length = 0;
double Pili::detach_grace_time = 0;

//
// pretty much arbitrary attachment timescale..
double Pili::dwell_time = 0.1;

////////////////////////////////////////////////////////////////////////////////
// bound extension timescale
double Pili::tau_e = 0.05;  // set by Clausen
// free extension timescale
double Pili::free_tau_e = 2.;  // free parameter
// free extension timescale / free retraction timescale
double Pili::free_tau_eovers = 1.;
// Clausen the proportions for elongation at low and high force.
double Pili::low_force_proportion = 3./100;
double Pili::high_force_proportion = 25./100;
//
double Pili::sh_to_elong_rate = 0.;


// 
bool Pili::allow_bound_extension = true;
bool Pili::allow_bound_ext_motor = true;
bool Pili::force_bound_retraction = false;
bool Pili::force_detachment = true;
bool Pili::simplify_free_pili_etor = false;
int Pili::rcondition = 1; 
bool Pili::enforce_stalling = false; 


// Important for WLC pili in the current implementation
bool Pili::force_max_length = false;

/////////////////////////////////////////////////////////////////////////////////

double pi = M_PI;

// initialisation
void Pili::common_init() {
  is_fully_retracted = (leq < Pili::inside_length) ? true : false;
  detach_time = -detach_grace_length;
  cycles = 1;
  new_cycle = false;
  free_time = 0.;
  bound_time = -1.;
  lifetime = 0.;
  isbound = 0;
  lastistaut = 0;
}

bool Pili::istaut() { return isbound && (pl > leq); }

// signed force
double Pili::sym_force_l() {
  double fmag = ks * (pl -leq)/(leq + Pili::force_min_length);
  return fmag;
}

// tension only force
double Pili::force_l() {
  assert (leq > 0.);
  // Set force to zero unless the pilus is taut.
  // In other words we don't consider the bending force here.

  if (!istaut())
  //if (false)
    return 0.;
  else
    return sym_force_l();
}

double Pili::elongate() {
  if (!this->isbound)
  {
    pl += d_free; leq += d_free;
  }
  else
  {
    leq += d_bound;
  }
  if (leq > Pili::inside_length + Pili::d_free)
  {
    is_fully_retracted = false;
    if (new_cycle) {
      cycles += 1;
      new_cycle = false;
    }
  }
  return leq;
}

double Pili::shrink_by(double delta) {
  if (!this->isbound)
  {
    pl -= delta; leq -= delta;
  }
  else
  {
    leq -= delta; 
  }
  this->_shrink();
  return leq;
}

double Pili::shrink() {

  if (!this->isbound)
  {
    pl -= d_free; leq -= d_free;
  }
  else
  {
    leq -= d_bound; 
  }
  this->_shrink();
  return leq;
}

double Pili::_shrink() {
  // prevent lengths less than inside_length - d_free
  if (leq < min_length) 
  { 
    // warning, pili may actually be bound to the surface close under the cell
    // in this case setting pl = leq seems resonable since the pilus must be tiny
    leq = min_length;
    pl = leq;
  }
  // count the extension/retraction cycles
  if (leq < Pili::inside_length && !is_fully_retracted) 
  {
    is_fully_retracted = true;
    new_cycle = true;
  }
}

void Pili::switch_se()
{
  if (ret_motor) {
    ret_motor = 0;
    ext_motor = 1;
  }
  else if (ext_motor) {
    ext_motor = 0;
    ret_motor = 1;
  }
}

bool Pili::has_motor() {
  return (ext_motor || ret_motor);
}

void Pili::ret_motor_on() {
  assert(!has_motor());
  ret_motor = 1;
}

void Pili::ret_motor_off() {
  assert(has_motor());
  ret_motor = 0;
}

void Pili::ext_motor_on() {
  assert(!has_motor());
  ext_motor = 1;
}
void Pili::ext_motor_off() {
  assert(has_motor());
  ext_motor = 0;
}


// sh_to_elong_rate is fixed by rtau and this other rate is variable
// i.e. probaility of switching to elongation is constant but
// probability of switching to retraction is force dependent
const double f1 = 8;
double Pili::tau_eovers(double fmag)
{
  if (fmag < f_stall) {
    double fexp = f1/Pili::f_stall;
    double Abase = Pili::low_force_proportion/pow(Pili::high_force_proportion, fexp);
    double A = pow(Abase, 1./(1-fexp));
    double B = log(Pili::high_force_proportion/A);
    return A * exp(B*fmag/f_stall);
  }
  // force maximum condition on switching rates
  else {
    // Set constant value for force > f_stall (we don't have data)
    return Pili::high_force_proportion;
  }
}

// add python test
double Pili::get_rate_s(double fmag)
{
  return tau_eovers(fmag)/tau_e;
}

// no longer force dependent
double Pili::calc_k_dt(double fmag) 
{
  return 1.0/dwell_time;
}
  
double Pili::energy() {
  double ens;
  if (!istaut())
    ens = 0.;
  else
    ens = ks/2. * (pl-leq)*(pl-leq)/(leq + Pili::min_length);
  return ens;
}

Vector3d Pili::get_attachment() {
  return this->attachment;
}

// recalculate pl 
void Pili::update(Vector3d anchor) {
  if (this->isbound) {
    Vector3d pili = this->attachment - anchor;
    this->pl = pili.len();
  }
  assert(false);
}

// Should be trivial since adh_accept = 1.
// no longer used in 3d
double Pili::attach_rate() {
  return k_adh * adh_accept;
}

double clamp(double n, double lower, double upper) {
  return n <= lower ? lower : n >= upper ? upper : n;
}

double Pili::get_anchor_angle() {
  if (!isbound) {
    return 0;
  }
  Vector3d pil = (attachment - lab_anchor).unit();
  double t = clamp(lab_axisEq.dot(pil), -1, 1);
  return std::acos(t);
}

bool Pili::check_rcondition() {
  if (rcondition == 1)
    return isbound;
  else if (rcondition == 2)
    return isbound && istaut();
  else
    throw std::runtime_error("Bad rcondition");
}

double Pili::get_bound_retract_fraction() {
  if (!isbound) {
    return 0.;
  }
  double theta = get_anchor_angle();
  if (theta <= anchor_angle_threshold) {
    return 1;
  }
  else {
    if (anchor_angle_smoothing_fraction == 0) {
      return 0;
    }
    else {
      double smx = M_PI/2 - anchor_angle_threshold;
      double a = anchor_angle_threshold;
      double b = anchor_angle_threshold + anchor_angle_smoothing_fraction * smx;
      return smooth(theta, a, b);
    }
  }
}

double Pili::get_bound_retract_rate() {
  return get_bound_retract_fraction() * kb_sh;
}

// break up this method for simplicity?
vector<double> Pili::get_rates()
{
  vector<double> prates {0.,0.,0.,0.,0.,0.,0.,0.,0.};

  // force magnitude for rates calculation
  //double fmag = force_l(); 

  double current_k_dt = 1/dwell_time;

  // elong, shrink, attach, release
  if (this->isbound == 1)
  {
    // check if detach is possible
    prates[action::Detach] = current_k_dt;
    if (ret_motor) {
      double _rate = get_bound_retract_rate();
      assert(!isnan(_rate));
      prates[action::Retraction] = _rate;
    }
    else if (ext_motor) {
      // extension is allowed if allow_bound_extension is true and the pilus is taut
      // if the pilus is not taut then we can't extend further (except by bending which is not yet modelled)
      prates[action::Extension] = (allow_bound_extension && istaut()) ? kb_ex : 0.; 
    }
  }
  else
  {
    // unbound
    prates[action::Resample] = this->get_k_resample();
    if (ret_motor) {
      prates[action::Retraction] = kf_sh;
    }
    else if (ext_motor) {
      prates[action::Extension] = kf_ex;
    }
  }
  // prevent shrinking if taut or at min length
  bool stall_condition = (Pili::enforce_stalling && abs(this->force_l()) >= Pili::f_stall);
  bool length_condition = leq == min_length;
  if (stall_condition || length_condition) {
    prates[action::Retraction] = 0.;
  } 

  if (ext_motor) {
    prates[action::Ext_off] = k_ext_off;
  }
  else if (ret_motor) {
    if (this->isbound) {
      prates[action::Ret_off] = force_bound_retraction ? 0. : k_ret_off;
    }
    else {
      prates[action::Ret_off] = k_ret_off;
    }
  }
  else if (!ret_motor && !ext_motor) {
    if (this->isbound) {
      prates[action::Ext_on] = allow_bound_ext_motor ? k_ext_on : 0.;
    }
    else {
      prates[action::Ext_on] = k_ext_on;
    }
    prates[action::Ret_on] = k_ret_on;
  }
  else {
    assert(false); // impossible
  }

  prates[action::Dissolve] = this->dissolve_on_retraction ? 0. : k_dissolve;

  return prates;
}

void Pili::_attach_upkeep(double tme)
{
  isbound = 1;
  bound_time = tme;

  if (!Pili::allow_bound_ext_motor) {
    ext_motor = 0;
  }
  if (Pili::force_bound_retraction){
    ext_motor = 0;
    ret_motor = 1;

  }
}

void Pili::detach(double tme)
{
  isbound = 0;
  pl = leq;
  detach_time = tme;
  lastistaut = 0;
  //
  free_time = tme;
}

std::shared_ptr<Pili> Pili::clone() {
  return std::make_shared<Pili>(*this);
}


void swap(Pili& lhs, Pili& rhs) {
  swap(lhs.idx,rhs.idx);
  swap(lhs.pl,rhs.pl);
  swap(lhs.leq,rhs.leq);
  swap(lhs.axisEq,rhs.axisEq);
  swap(lhs.anchor,rhs.anchor);
  swap(lhs.attachment,rhs.attachment);
  swap(lhs.cycles,rhs.cycles);
  swap(lhs.isbound,rhs.isbound);
  swap(lhs.lastistaut,rhs.lastistaut);
  swap(lhs.is_fully_retracted,rhs.is_fully_retracted);
  swap(lhs.bound_time,rhs.bound_time);
  swap(lhs.free_time,rhs.free_time);
  swap(lhs.lifetime,rhs.lifetime);
  swap(lhs.dissolve_on_retraction,rhs.dissolve_on_retraction);
  swap(lhs.ret_motor,rhs.ret_motor);
  swap(lhs.ext_motor,rhs.ext_motor);
}

std::string Pili::__str__() {
  return "Object Pili -- __str__ not implemented";
}

////////////////////////////////////
//piliRod3d

void Pili3d::update(Vector3d lab_anchor) {
  Vector3d pili = this->attachment - lab_anchor;
  this->pl = pili.len();
}


void swap(Pili3d& lhs, Pili3d& rhs) {
  swap(static_cast<Pili&>(lhs), static_cast<Pili&>(rhs));
}

////////////////////////////////////
//piliRod3d
//


// 3d attachment
// when kmcstep calls this we just check for an attachment. if we find one we 
// update this->pl and call _attach_upkeep()

boost::optional<Vector3d> PiliRod3d::check_attach(Vector3d lab_anchor, Capsule body, Shape* surface)
{
  // Get lab position and direction of the pili to check for collision with the lab surface
  Vector3d lab_axisEq = body.frame.to_lab_rt(this->axisEq); 
  Lvec lv(lab_anchor, leq*lab_axisEq);
  boost::optional<Vector3d> target = surface->intersects(lv);
  if (target) {
    this->last_inter_len = lv.len() - lv.len_to_inter(); 
  }
  this->intersects_surface = bool(target);
  return target;
}

// take the capsule so we can obtain the real position of the pili on the cell head
// why does this take Shape* as argument but cell3d holds Plane* object
int PiliRod3d::attach(double tme, Capsule body, Shape* surface) {
  Vector3d lab_anchor = body.frame.to_lab(this->anchor);
  boost::optional<Vector3d> target = this->check_attach(lab_anchor, body, surface);

  if (!target) {
    return 0; // attachment failed
  }
  else {
    // set this in pili.update?
    this->attachment = *target;
    this->pl = (this->attachment - lab_anchor).len();
    _attach_upkeep(tme);
    return 1;
  }
}

std::string pilstring(std::string classname, Pili& p) {
  std::string form = 
    "%s: idx(%d)\n->pl(%12.6f)\n->leq(%12.6f)\n->anchor%s\n->axisEq%s\n->ret_motor(%d)\n->ext_motor(%d)\n->isbound(%d)";
  return ( boost::format(form) 
      % classname
      % p.idx
      % p.pl
      % p.leq
      % p.anchor.__str__()
      % p.axisEq.__str__()
      % p.ret_motor
      % p.ext_motor
      % p.isbound
      ).str();
}

std::string PiliRod3d::__str__() {
  return pilstring("PiliRod3d", *this);
}

std::shared_ptr<Pili> PiliRod3d::clone() {
  // copy wlc shared pointer because this object is constructed by python
  return std::shared_ptr<Pili>(new PiliRod3d(*this));
}

/////////////////////////////////////////////////////////////////////////////////


void PiliWLC3d::sample() {
  // build a new chain in body frame
  this->n = (int)(leq/a) + 1;
  lasta = leq - (n-1)*a;
  this->chain_theta = wlc->sample_theta();
  this->chain_phi = wlc->sample_phi();
  // std::unique_ptr<Chain> sample_ch = this->wlc->sample(this->n, this->axisEq);
  // sample_ch->set_lasta(get_lasta());
  // this->set_chain_instance(std::move(sample_ch));
}

// dep
// merge this method with the next one
std::shared_ptr<Chain> PiliWLC3d::new_chain_instance(void)
{
  // construct chain, this is not a new chain
  Chain* ch;
  int nbend = this->n - 1;
  if (nbend == 0) {
    ch = new Chain(this->axisEq, this->a);
  }
  else if (nbend > 0) {
    ch = new Chain(this->axisEq, this->a, (*chain_theta)(Eigen::seq(0,nbend-1)), 
      (*chain_phi)(Eigen::seq(0,nbend-1)) );
  }
  return std::shared_ptr<Chain>(ch);
}

std::shared_ptr<Chain> PiliWLC3d::transformed_chain_instance(Capsule body)
{
  Vector3d lab_anchor = body.frame.to_lab(this->anchor);
  //  this not a new chain, just a new instance
  std::shared_ptr<Chain> ch = new_chain_instance();
  ch->set_lasta(get_lasta());
  ch->rotate(body.frame.get_rmatrix());
  ch->set_source(lab_anchor);
  ch->compute_targets();
  return ch;
}

double PiliWLC3d::shorten_to_detach(Capsule body, Shape* surface) {
  std::shared_ptr<Chain> ch = transformed_chain_instance(body);
  boost::optional<Vector3d> target = surface->intersects(*ch);
  if (target) {
    double full_length = ch->get_length();
    double new_length = ch->shorten(ch->inter_n, ch->inter_t);
    double shortening = full_length - new_length;
    return shortening;
  }
  return 0;

}

boost::optional<Vector3d> PiliWLC3d::check_attach(Vector3d lab_anchor, Capsule body, Shape* surface)
{
  //  this not a new chain, just a new instance
  std::shared_ptr<Chain> new_ch = transformed_chain_instance(body);
  Chain& ch = *new_ch;

  // intersection point is computed
  boost::optional<Vector3d> target = surface->intersects(ch);
  if (target) {
    // if no intersection then last_inter_len has no meaning
    this->last_inter_len = ch.get_length() - ch.len_to_inter();
  }
  this->intersects_surface = bool(target);
  return target;
}

int PiliWLC3d::attach(double tme, Capsule body, Shape* surface) {
  Vector3d lab_anchor = body.frame.to_lab(this->anchor);
  boost::optional<Vector3d> target = this->check_attach(lab_anchor, body, surface);
  if (!target) {
    return 0; // attachment failed
  }
  else {
    this->_attach(tme, lab_anchor, *target);
    return 1;
  }
}

// to emulate attachment just set 
// this->attachment, this->pl and call _attach_upkeep(time)
// also need to update the bound and free lists in cell3d object
int PiliWLC3d::_attach(double tme, Vector3d lab_anchor, Vector3d target) {
  this->attachment = target;
  this->pl = (this->attachment - lab_anchor).len();
  _attach_upkeep(tme);
  return 1;
}


double PiliWLC3d::elongate() {
  Pili3d::elongate();
  if (!this->isbound)
  {
    n = int(leq/a) + 1;
    lasta = leq - (n-1)*a;
  }
  // forget about segments for bound pili
  return leq;
}

double PiliWLC3d::shrink_by(double delta) {
  Pili3d::shrink_by(delta);
  n = int(leq/a) + 1;
  lasta = leq - (n-1)*a;

  return leq;
}

double PiliWLC3d::shrink() {
  Pili3d::shrink();
  // additional bookkeeping
  n = int(leq/a) + 1;
  lasta = leq - (n-1)*a;
  // this check should be unnecessary (?)
  if (lasta == 0.) {
    // length is exact multiple of a.
    n--;
    lasta = a;
  }

  // TODO check n is safely satisfies n > 0
  // if (n < 1) { n = 1; }
  
  return leq;
}


void PiliWLC3d::detach(double tme)
{
  Pili3d::detach(tme);
  // should move shrink() call to pili3d
  if (detach_grace_length > 0) {
    this->shrink_by(Pili::detach_grace_length);
  }
  else {
    this->shrink();
  }
  this->n = int(leq/a) + 1;
  this->lasta = leq - (n-1)*a;
}


std::string PiliWLC3d::__str__() {
  return pilstring("PiliWLC3d", *this);
}

std::shared_ptr<Pili> PiliWLC3d::clone() {
  // copy wlc shared pointer because this object is constructed by python
  return std::shared_ptr<Pili>(new PiliWLC3d(*this));
}

// copy constructor
PiliWLC3d::PiliWLC3d(const PiliWLC3d& other) 
  : Pili3d(other), ka(other.ka), a(other.a), n(other.n), lasta(other.lasta) {
  // copy wlc shared pointer because this object is constructed by python
  wlc = other.wlc;
  // how to copy a std::unique_ptr<Eigen::VectorXd> chain_theta;
  chain_theta = std::unique_ptr<Eigen::VectorXd>( new Eigen::VectorXd(*other.chain_theta) );
  chain_phi = std::unique_ptr<Eigen::VectorXd>( new Eigen::VectorXd(*other.chain_phi) );

  
  // use Chain copy constuctor
  // check copy is safe TODO
  // chain_instance = std::make_unique<Chain>(*other.chain_instance);
}

//http://www.cplusplus.com/articles/y8hv0pDG/p
PiliWLC3d& PiliWLC3d::operator=( PiliWLC3d tmp ) {
  swap(*this, tmp);
  return *this;
} 

//friend void swap(PiliWLC3d& lhs, PiliWLC3d& rhs);
void swap(PiliWLC3d& lhs, PiliWLC3d& rhs) {
  // check this todo
  swap(static_cast<Pili3d&>(lhs), static_cast<Pili3d&>(rhs));
  swap(lhs.ka, rhs.ka);
  swap(lhs.a, rhs.a);
  swap(lhs.n, rhs.n);
  swap(lhs.lasta, rhs.lasta);
  swap(lhs.wlc, rhs.wlc);
  // swap(lhs.chain_instance, rhs.chain_instance);
}



//////////////////////////////////////////////////////////////////////////////////////


