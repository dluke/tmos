
#include "plane.hpp"
#include <boost/format.hpp>
#include "sphere.hpp"

// watch out for circular references
#include "cell.hpp"

using std::abs;

double INF = 1e10;

std::string Plane::__str__() {
  std::string form = "Plane: origin%s normal%s";
  return ( boost::format(form) % frame.origin.__str__() % frame.e3.__str__() ).str();
}


std::ostream& operator<<(std::ostream& os, Plane& pl)
{
  os << pl.__str__();
  return os;
}

// generic
boost::optional<Vector3d> Plane::intersects(Lvec& lv) {
  double t = -(N().dot(lv.source - get_origin())) / (N().dot(lv.line));
  if (t <= 0 || t > 1) { 
    return boost::none;
  }
  else {
    lv.inter_t = t;
    Vector3d ipt = lv.get_pt(t);
    return ipt;
  }
}

boost::optional<Vector3d> Plane::intersects(Chain& ch) {
  vector<Lvec> lc = ch.line_components();
  for (int i=0; i < lc.size(); i++ )
  {
    Lvec lv = lc[i];
    boost::optional<Vector3d> inter = this->intersects(lv);
    if (inter) {
      ch.inter_n = i;
      ch.inter_t = lv.inter_t;
      return inter;
    }
  }
  return boost::none;
}


int Plane::get_num_contacts(Capsule& body)
{
  int ncontacts = 0;
  if (abs(body.get_headpt().z) < body.R) {
    ncontacts += 1;
  }
  if (abs(body.get_endpt().z) < body.R) {
    ncontacts += 1;
  }
  return ncontacts;
}


// this method is avoided by direct implementation of distance calculation
vector<overlap> Plane::overlap_vector(Capsule& body)
{
  vector<overlap> overs;
  double svd = body.get_headpt().z;
  double endvd = body.get_endpt().z;
  double m = body.length/2.;
  if (svd < body.R) {
    overlap over = std::make_tuple(m, svd, N());
    overs.push_back(over);
  }
  if (endvd < body.R) { // throw the critical case in here
    overlap over = std::make_tuple(-m, endvd, N());
    overs.push_back(over);
  }
  return overs;
}

double Plane::surface_energy(ACell& cell, contactf energyf, bool cost_anchor_intersection) 
{
  if (cost_anchor_intersection && this->containsany(cell.get_bound_anchors())) {
    return INF;
  }
  return _surface_energy(cell, energyf);
}


double Plane::_surface_energy(ACell& cell, contactf energyf) {
  Capsule& body = cell.get_body();
  double energy = 0;
  vector<Vector3d> bodyends{body.get_headpt(), body.get_endpt()};
  for (Vector3d vpt : bodyends) {
    energy += energyf(abs(vpt.z));
  }
  return energy;
}


Eigen::RowVectorXf Plane::surface_grad(Capsule& body, contactf contact) {
  Vector3d fgrad;
  double fg; 
  vector<double> pk = {0,0,0};

  Rk bodyRk = body.get_Rk();

  double ljg;
  // head contact
  ljg = -contact(abs(body.get_headpt().z));
  fgrad += ljg * N();
  for (int k=0; k<3; k++) {
    pk[k] += dot(ljg * body.length/2. * e_z, bodyRk.Rkvez[k]);
  }
  // tail contact
  ljg = -contact(abs(body.get_endpt().z));
  fgrad += ljg * N();
  for (int k=0; k<3; k++) {
    pk[k] += dot(ljg * -body.length/2. * e_z, bodyRk.Rkvez[k]);
  }

  Eigen::RowVectorXf grad(6);
  grad << fgrad.x, fgrad.y, fgrad.z, pk[0], pk[1], pk[2];
  return grad;
}

// generic
Eigen::RowVectorXf Plane::surface_grad_part(Capsule& body, contactf contact, overlap over) 
{
  vector<double> pk = {0,0,0};
  Rk bodyRk = body.get_Rk();
  double m = get<0>(over);
  double r = get<1>(over);
  Vector3d rhat = get<2>(over);

  double fg = -contact(r);
  Vector3d fgrad = fg*rhat;
  for (int k=0; k<3; k++) {
      pk[k] = fg * m * dot(rhat, bodyRk.Rkvez[k]);
  }
  Eigen::RowVectorXf grad(6);
  grad << fgrad.x, fgrad.y, fgrad.z, pk[0], pk[1], pk[2];
  return grad;
}




// move to base class
bool Plane::containsany(vector<Vector3d> pts) {
  for (Vector3d pt : pts) {
    if (this->contains(pt))
      return true;
  }
  return false;
}


//////////////////////////////////////////////////////////////////////////////////
// General Plane

vector<overlap> AnyPlane::overlap_vector(Capsule& body)
{
  // assert that body is outside surface

  Vector3d sv =  this->distance_vector(body.get_headpt());
  Vector3d endv = this->distance_vector(body.get_endpt());
  // distance from overlap point/s to origin
  double svd = sv.len();
  double endvd = endv.len();
  // check lengths are not zero
  vector<overlap> overs;
  double R = body.R;
  double m = body.length/2.;

  if (svd < R) {
    overlap over = std::make_tuple(m, svd, sv.unit());
    overs.push_back(over);
  }
  if (endvd < R) { // throw the critical case in here
    overlap over = std::make_tuple(-m, endvd, endv.unit());
    overs.push_back(over);
  }
  return overs;
}

// generic method
int AnyPlane::get_num_contacts(Capsule& body)
{
  //for (auto over : overlap_vector(body)) {
    //cout << over << endl;
  //}
  return overlap_vector(body).size();
}

// generic method
double AnyPlane::_surface_energy(ACell& cell, contactf energyf) {
  Capsule& body = cell.get_body();
  double energy = 0;
  vector<overlap> inter = this->overlap_vector(body);
  for (overlap over : inter) {
    double r = get<1>(over);
    energy += energyf(r);
  }
  return energy;
}



Eigen::RowVectorXf AnyPlane::surface_grad(Capsule& body, contactf contact) 
{
  vector<overlap> inter = this->overlap_vector(body); 
  Eigen::RowVectorXf grad = Eigen::RowVectorXf::Zero(6);

  for (overlap over : inter) {
    grad += this->surface_grad_part(body, contact, over);

  }

  return grad;
}

// Generic Partial Plane

bool PartialPlane::on_plane(Vector3d pt) 
{ 
  double x = proj_limited(pt); 
  return (x > lmin && x < lmax);
}

boost::optional<Vector3d> PartialPlane::intersects(Lvec& lv) 
{
  boost::optional<Vector3d> inter = AnyPlane::intersects(lv);
  if (inter) {
    double x1 = this->proj_limited(*inter);
    if (x1 >= lmin && x1 <= lmax) {
      return inter;
    }
  }
  return boost::none;
}


vector<overlap> PartialPlane::overlap_vector(Capsule& body)
{
  vector<overlap> overs;
  Vector3d hp = body.get_headpt();
  Vector3d tp = body.get_endpt();
  // project onto plane then project onto limited directino
  double hx1 = this->proj_limited(this->project(hp));
  double tx1 = this->proj_limited(this->project(tp));
  bool range_condition_h = hx1 > lmin && hx1 < lmax;
  bool range_condition_t = tx1 > lmin && tx1 < lmax;
  if (!range_condition_h && !range_condition_t) {
    return overs; // exit early
  }
  // the distance between head pt and plane
  double svd;
  svd = this->distance(hp);
  double m = body.length/2.;
  if (range_condition_h && svd < body.R)
  {
    overlap over = std::make_tuple(m, svd, N());
    overs.push_back(over);
  }

  svd = this->distance(tp);
  if (range_condition_t && svd < body.R) {
    overlap over = std::make_tuple(-m, svd, N());
    overs.push_back(over);
  }
  return overs;
}


// Specific Partial Planes


boost::optional<Vector3d> PartialPlaneX::intersects(Lvec& lv) {
  // first compute contact then check range_conditions
  boost::optional<Vector3d> inter = AnyPlane::intersects(lv);
  if (inter) {
    // range condition
    if ( (*inter).x >= xmin && (*inter).x <= xmax ) {
      return inter;
    }
  }
  return boost::none;
}

boost::optional<Vector3d> PartialPlaneZ::intersects(Lvec& lv) {
  // first compute contact then check range_conditions
  boost::optional<Vector3d> inter = AnyPlane::intersects(lv);
  if (inter) {
    // range condition
    if ( (*inter).z >= zmin && (*inter).z <= zmax ) {
      return inter;
    }
  }
  return boost::none;
}


vector<overlap> PartialPlaneX::overlap_vector(Capsule& body)
{
  vector<overlap> overs;
  Vector3d hp = body.get_headpt();
  Vector3d tp = body.get_endpt();
  bool range_condition_h = hp.x > xmin && hp.x < xmax;
  bool range_condition_t = tp.x > xmin && tp.x < xmax;
  if (!range_condition_h && !range_condition_t) {
    return overs; // exit early
  }
  double svd = hp.z - get_origin().z;
  double m = body.length/2.;
  if (range_condition_h && svd < body.R)
  {
    overlap over = std::make_tuple(m, svd, N());
    overs.push_back(over);
  }

  svd = tp.z - get_origin().z;
  if (range_condition_t && svd < body.R) {
    overlap over = std::make_tuple(-m, svd, N());
    overs.push_back(over);
  }
  return overs;
}

vector<overlap> PartialPlaneZ::overlap_vector(Capsule& body)
{
  vector<overlap> overs;
  Vector3d hp = body.get_headpt();
  Vector3d tp = body.get_endpt();
  bool range_condition_h = hp.z > zmin && hp.z < zmax;
  bool range_condition_t = tp.z > zmin && tp.z < zmax;


  if (!range_condition_h && !range_condition_t) {
    return overs; // exit early
  }
  // our vertical planes are going to be both +x facing  and -x facing alternately
  double xnormal = N().x;
  assert(xnormal == 1. || xnormal == -1.);
  double svd = xnormal * (hp.x - get_origin().x);
  double m = body.length/2.;
  if (range_condition_h && svd < body.R)
  {
    overlap over = std::make_tuple(m, svd, N());
    overs.push_back(over);
  }

  svd = xnormal * (tp.x - get_origin().x);
  if (range_condition_t && svd < body.R) {
    overlap over = std::make_tuple(-m, svd, N());
    overs.push_back(over);
  }
  return overs;
}

// FOR SPECIAL CASE
vector<overlap> PartialPlaneZ::end_overlap(Capsule& body)
{
  Vector3d toppt = get_origin();
  toppt.z = zmax;

  // should be a distance in N() direction ...
  Lvec bodylv = body.get_lvec();
  double t = bodylv.xzproject().closest_t(toppt);

  if (t == 0. || t == 1.) {
    return {}; // should have found endpoint contacts already
  }

  double m = (0.5-t) * body.length;
  Vector3d lpt = bodylv.get_pt(t);
  double svd = N().x * (lpt.x - toppt.x); // body point - surface point
  assert(svd > 0.);
  if (svd < body.R) {
    return {std::make_tuple(m, svd, N())};
  }
  return {};
}


std::string PartialPlaneX::__str__() {
  std::string form = "PartialPlaneX: origin%s normal%s xmin(%f) xmax(%f)";
  return ( boost::format(form) % get_origin().__str__() % N().__str__() % xmin % xmax).str();
}


std::string PartialPlaneZ::__str__() {
  std::string form = "PartialPlaneZ: origin%s normal%s zmin(%f) zmax(%f)";
  return ( boost::format(form) % get_origin().__str__() % N().__str__() % zmin % zmax).str();
}

std::string NullPlane::__str__() {
  std::string form = "NullPlane";
  return form;
}

std::string Plane::_class() {return "Plane";}
std::string AnyPlane::_class() {return "AnyPlane";}
std::string PartialPlaneX::_class() {return "PartialPlaneX";}
std::string PartialPlaneZ::_class() {return "PartialPlaneZ";}
