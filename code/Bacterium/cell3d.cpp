

#include "cell3d.hpp"

#include <algorithm>
#include <boost/format.hpp>

/////////////////////////////////////////////////////////////////////////
//3d cell subclass

// surface interaction constant
double Cell3d::eps = 1;  
bool Cell3d::repulsive_only = false; 
// size of the attractive well, distance in micrometres
double Cell3d::attract_width = 0.01; 
// check the anchors never intersect the surface
bool Cell3d::cost_anchor_intersection = false;
std::string Cell3d::distrib_type = "gaussian";

//////////////////////////////////////////////////////////////////////////////////

vector<Vector3d> Cell3d::get_bound_anchors() {
  vector<Vector3d> bpts;
  for (auto pilus : this->pili)
  {
    if (pilus->isbound)
      bpts.push_back(this->get_lab_anchor(*pilus)); 
  }
  return bpts;
}

Cell3d& Cell3d::operator=( Cell3d tmp) {
  using std::swap;
  swap(*this, tmp);
  return *this;
}

void swap(Cell3d& lhs, Cell3d& rhs) {
  using std::swap;
  swap(static_cast<ACell&>(lhs), static_cast<ACell&>(rhs));
  swap(lhs.body, rhs.body);
  swap(lhs.surface, rhs.surface);
}


// needs carefully testing
// replace with sampling of Kent distribution
double _sample_spherical(Capsule body, double pvar)
{
    double umax = 1. - pvar/(body.circ()/4.);
    std::uniform_real_distribution<double> nd(umax, 1.);
    double u = nd(mtgen);
    double theta = std::acos(u);
    return theta;
}


//http://mathworld.wolfram.com/SpherePointPicking.html
anchorbits _uniform_pili_body_distrib(Capsule& body, double pilivar, double polarisation)
{

  const double _tol = 0.0001;
  if (ACell::pilivar > body.halfal() + _tol)
  { throw std::runtime_error("ACell::pilivar should be less than half the axis length"); }

  double phi = r2() * 2 * M_PI;
  double theta;
  double d;

  int seg = body.get_segment(ACell::pilivar);

  if (seg == body.HEAD) {
    // spherical part only
    theta = _sample_spherical(body, ACell::pilivar);
    d = theta * body.R;
  }
  else if ( seg == body.BODY )
  {
    double blength = ACell::pilivar - body.circ()/4.;
    double spharea = body.spharea()/2.;
    double cyarea = blength * 2 * M_PI;
    if ( r2() < spharea/(cyarea + spharea) ) {
      // sphere part
      theta = _sample_spherical(body, body.circ()/4.);
      d = theta * body.R;
    }
    else {
      // cylinder part
      d = r2() * blength + body.circ()/4.;
    }
  }
  else { assert(false); }
  // m * body.axis is the component of anchor along the body axis
  double m = (d > body.circ()/4.) ? body.length/2. : d - body.circ()/4.;

  // polarisation can be in range [0., 1.]
  if ( (r2()*2 - 1) > polarisation) {
    d = body.long_length() - d;
    m = -m;
    theta = -(M_PI - theta);
    // phi no change after mirror
  }
  std::pair<Vector3d, Vector3d> ret = body.surface_pt(d, phi, Pili::inside_length);

  // tmp
  if ( !(abs( abs(body.sgn_distance_to_surface(ret.first)) - Pili::inside_length) < _tol) )
   throw std::runtime_error(
       (boost::format("Failed to generate a pili anchor close to the body surface. distance %f"
                      "\nAnchor point %s")
        % abs(body.sgn_distance_to_surface(ret.first))
        % ret.first.__str__()
       ).str()
       );

  return anchorbits{ret.first, ret.second};
}

anchorbits _gaussian_pili_body_distrib(Capsule& body, double pilivar, double polarisation) 
{
  Vector3d X = modified_vmf(pilivar, body.length, body.R-Pili::inside_length);
  if ( (r2()*2 - 1) > polarisation) {
    // flip in z axis
    X.z = -X.z;
  }
  Capsule refbody{Vector3d(), e_z, body.R, body.length};
  Vector3d cpt = refbody.get_lvec().closest_pt(X);
  Vector3d axis = (X - cpt).unit();
  return anchorbits{X, axis};
}


anchorbits Cell3d::pili_body_distrib()
{
  if (Cell3d::distrib_type == "gaussian") {
    //cout << "using gaussian distribution" << endl;
    return _gaussian_pili_body_distrib(get_body(), pilivar, polarisation);
  }
  else if (Cell3d::distrib_type == "uniform") {
    //cout << "using uniform distribution" << endl;
    return _uniform_pili_body_distrib(get_body(), pilivar, polarisation);
  }
  else {
    assert(false);
  }
}


Vector3d Cell3d::pili_body_position()
{
  return pili_body_distrib().anchor;
}

std::shared_ptr<Pili> Cell3d::spawn_pilus()
{
  anchorbits abits = pili_body_distrib();
  double pl = Cell3d::pili_min_length;
  std::shared_ptr<PiliRod3d> pilus (new PiliRod3d(
        pili_counter, pl, pl, abits.anchor, abits.axis, 0, spawn_extension_state));
  pilus->min_length = pl;
  return pilus;
}

Vector3d Cell3d::get_headpt()  {
  return this->body.get_headpt();
}
Vector3d Cell3d::real_pos() 
{
  return this->body.get_headpt();
}
Vector3d Cell3d::get_trail() {
  return this->body.get_endpt();
}
Vector3d Cell3d::real_centre() {
  return this->body.get_centerpt();
}
Vector3d Cell3d::get_vector() {
  return this->body.get_vector();
}
Vector3d Cell3d::get_axis()
{
  return this->body.get_axis();
}

CellEvent* Cell3d::attach(Pili& pilus, double tme) {
  if (tme - pilus.detach_time < Pili::detach_grace_time) {
    // We only just detached the surface so prevent any attachment 
    // return no_attach cell event
    return nullptr;
  }
  int success = pilus.attach(tme, this->body, get_surface());
  if (success)
  {
    CellEvent* e = new CellEvent(tme, pilus.idx, string(""), 
        std::move(this->clone()), string("attach"));
    return e;
  } 
  return nullptr;
}

CellEvent* Cell3d::detach(Pili& pilus, double tme) 
{
  // check if the pilus still intersects the surface
  boost::optional<Vector3d> target;

  if (Pili::force_detachment) {
    // failed detachment not permitted
    // hack the pilus down to size if necessary
    double shortening = pilus.shorten_to_detach(body, surface);
    pilus.shrink_by(shortening + pilus.d_bound);
    pilus.last_shortening = shortening;
    target = boost::none;
  }
  else {
    // shrink the pilus and attempt to detach it that way, detachment allowed to fail
    Vector3d lab_anchor = this->body.frame.to_lab(pilus.anchor);
    double shrink_length = (Pili::detach_grace_length > 0.) ? Pili::detach_grace_length : Pili::d_bound;
    // create temporary pilus 
    auto ptemp = pilus.clone();
    // shorten the length before checking attachment
    ptemp->shrink_by(shrink_length);
    target = ptemp->check_attach(lab_anchor, body, surface);
    pilus.last_inter_len = ptemp->last_inter_len;
  }

  pilus.intersects_surface = bool(target);

  CellEvent* e;
  if (target)
  {
    // no detach occurs because the pilus is still intersecting the surface 
    // we make no change to the system except update pilus.last_inter_len
    // It's important to log these failed detachments because we want to minimise them
    // so that the \tau_dwell parameter and the real \tau_dwell converge as closely as possible
    e = new CellEvent(tme, pilus.idx, string("release"), this->clone(), string("no_release"));
  }
  else {
    // actually perform detachment
    pilus.detach(tme);
    // record detach event
    if (pilus.is_fully_retracted) {
      // dissolve before attachstep
      e = new CellEvent(tme, pilus.idx, string("release"), this->clone(), string("dissolve"));
    }
    else {
      // usual case, pilus exists in detached state
      e = new CellEvent(tme, pilus.idx, string("release"), this->clone(), string("release"));
    }
    
  }
  return e;
}

vector<fpt> Cell3d::pili_force()
{
  Vector3d cpt = this->body.get_centerpt();
  // vector container has move semantics
  vector<fpt> vfpt;

  for (auto p : pili)
  {
    Pili& pilus = *p;
    if (!pilus.isbound) {
      continue;
    }
    // compute axis by obtaining anchor point in lab frame
    Vector3d lab_anchor = get_lab_anchor(pilus);
    Vector3d axis = (pilus.attachment - lab_anchor).unit();
    // compute torque
    Vector3d force = pilus.force_l() * axis;
    Vector3d torque = (lab_anchor - cpt).cross(force);
    fpt ft = { force, torque };
    vfpt.push_back(ft);
  }
  return vfpt;
}

Eigen::RowVectorXf Cell3d::pili_grad() {
  // setup output vectors
  Vector3d fgrad;
  vector<double> pk = {0,0,0};
  double phi, thetaconst;

  // pili indepenent calculation
  Capsule body = this->get_body();
  Rk bodyRk = body.get_Rk();

  for (auto p : pili)
  {
    Pili& pilus = *p;
    // setup
    // get pili anchor in lab frame, \theta' and \phi
    if (!pilus.isbound) {
      continue;
    }
    double fg = pilus.force_l();
    if (fg == 0.) 
      continue; // gradient is 0, skip extra computation
    Vector3d lab_anchor = get_lab_anchor(pilus);
    // axis = r_a - r_b
    Vector3d axis = (lab_anchor - pilus.attachment).unit();
    // gradient of position
    fgrad += fg * axis; // (-1) in force and axis cancel

    // gradient of body axis
    for (int k=0; k<3; k++) {
      pk[k] += fg * dot(bodyRk.Rkv[k]*pilus.anchor, axis);
    }
  }
  Eigen::RowVectorXf grad(6);
  grad << fgrad.x, fgrad.y, fgrad.z, pk[0], pk[1], pk[2];
  return grad;
}

Eigen::RowVectorXf Cell3d::surface_grad() {
  return this->get_surface()->surface_grad(this->get_body(), this->sforce);
}


// use Eigen to avoid this clunky code for adding std::arrays
Eigen::RowVectorXf Cell3d::grad() {
  Eigen::RowVectorXf pgrad = this->pili_grad();
  Eigen::RowVectorXf sgrad = this->surface_grad();
  return pgrad + sgrad;
}

void Cell3d::set_state(Eigen::RowVectorXf state) {
  Vector3d cpt = Vector3d(state[0], state[1], state[2]);
  Vector3d p = Vector3d(state[3], state[4], state[5]);
  this->get_body().set_origin(cpt);
  this->set_rotation(p);
  // update pili lengths with the new position
  this->update_pili();
}

Eigen::RowVectorXf Cell3d::get_state(void) {
  Vector3d rc = real_centre();
  Vector3d ax = get_body().get_axisangle();
  Eigen::RowVectorXf state(6);
  state << rc.x,rc.y,rc.z,ax.x,ax.y,ax.z;
  return state;
}

double Cell3d::state_energy(Eigen::RowVectorXf state) {
  this->set_state(state);
  return this->energy();
}
Eigen::RowVectorXf Cell3d::state_gradient(Eigen::RowVectorXf state) {
  this->set_state(state);
  return this->grad();
}

// may dep
vector<fpt> Cell3d::surface_force()
{
  //cout << "surface forces" << endl;
  // surface interaction force is determined by overlaps 
  vector<fpt> sfpt;
  vector<overlap> inter = this->get_surface()->overlap_vector(this->get_body());
  Vector3d force;
  for (overlap over : inter) {
    // changed meaning of overlap object
    // check this calculation
    double m = get<0>(over);
    double r = get<1>(over);
    Vector3d rv = r * get<2>(over);
    Vector3d r_a = m * get_body().get_axis() - rv;

    if (repulsive_only) {
      force = abs(wca_force_abs(this->eps, body.R, r)) * rv.unit();
    }
    else {
      force = this->sforce(r) * rv.unit();
    }
    Vector3d torque = (r_a - this->body.get_centerpt()).cross(force);

    fpt ft = {force, torque};
    sfpt.push_back(ft);
  }
  return sfpt;
}

fpt Cell3d::total_ft()
{
  vector<fpt> headft = this->pili_force();

  //cout << "get surface forces" << endl;
  vector<fpt> surfaceft = this->surface_force();

  // add up the forces, they will displace the capsule head
  Vector3d ftotal = Vector3d();
  Vector3d ttotal = Vector3d();
  for (fpt f : headft) {
    ftotal += f.force;
    ttotal += f.torque;
  }
  //cout << "total pili force " << ftotal << endl;
  //cout << "total pili torque" << ttotal << endl;
  for (fpt f : surfaceft) {
    ftotal += f.force;
    ttotal += f.torque;
  }
  //cout << "total pili+surface force " << ftotal << endl;
  //cout << "total pili+surface torque" << ttotal << endl;
  fpt ft = {ftotal, ttotal};
  return ft;
}


 
Capsule& Cell3d::get_body() { return this->body; }

void Cell3d::update_pos(Vector3d dxy) 
{
  this->body.translate(dxy);
}

// fix efficiency (?)
// todo cut or fix this method
void Cell3d::update_axis(Matrix3d deltaR)
{
  // get rmatrix from the body
  Matrix3d R = this->body.frame.get_rmatrix();
  // This R matrix could become nearly singular ("euler angle problem")

  // use R matrix to update axis component of body frame
  Matrix3d Rdash = R + deltaR;
  
  // check othonormality
  Rdash.svdnormalise();

  //this->set_rotation(Rdash);
  // additionally enforce orthonormality
  this->body.frame.cross_norm();
}

void Cell3d::set_rotation(Vector3d p)
{
  this->get_body().set_rotation(p);
}

void Cell3d::update_pili() {
  // just iterate through bound pili list
  for (auto pilus : pili)
  {
    if (pilus->isbound)
      // pili3d types which subclass pili should have an update method like this
      pilus->update(get_lab_anchor(*pilus));
  }
}


void Cell3d::update(Vector3d dxy, Matrix3d deltaR) {
  // order of pos and angle updates shouldn't matter
  this->update_pos(dxy);
  this->update_axis(deltaR);
  this->update_pili();
}


//
std::string Cell3d::__str__() {
  
  // introduce a separating line between multiline child objects
  std::string form = "Cell3d: idx(%d)\n->%s\n--\n->%s";
  std::string cell_string = ( boost::format(form) 
      % idx
      % body.__str__()
      % get_surface()->__str__()
      ).str();

  for (auto& p : this->pili) {
    cell_string += "\n--\n";
    cell_string += p->__str__();
  }
  // bonus, print the anchor orientations in the lab frame
  cell_string += "\n---anchor orientations---\n";
  for (auto& p : this->pili) {
    cell_string += std::to_string(p->idx) + " " + get_lab_axisEq(*p).__str__() + "\n";
  }

  return cell_string;
}


//

double Cell3d::energy_surface()
{
  return get_surface()->surface_energy(*this, this->senergy);
}

double Cell3d::energy()
{
  return energy_pili() + energy_surface();
}

//////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Pili> CellWLC3d::spawn_pilus() {
  anchorbits abits = pili_body_distrib();
  double pl = this->pili_min_length;
  std::shared_ptr<PiliWLC3d> pilus (new PiliWLC3d(pili_counter,
        pl, pl, abits.anchor, abits.axis, 0, spawn_extension_state, wlc));
  pilus->min_length = pl;
  return pilus;
}

void swap(CellWLC3d& lhs, CellWLC3d& rhs) {
  using std::swap;
  swap(static_cast<Cell3d&>(lhs), static_cast<Cell3d&>(rhs));
  swap(lhs.wlc, rhs.wlc);
}

CellWLC3d& CellWLC3d::operator=( CellWLC3d tmp) {
  using std::swap;
  swap(*this, tmp);
  return *this;
}

