

#include "cylinder.hpp"

#include <iostream>
#include <exception>
#include <cassert> // needed?
// for overlaps
#include <utility>
#include <tuple>

#include <boost/format.hpp>

double Circle::distance(Vector3d pt) {
  double xd = pt.x - origin.x;
  double zd = pt.z - origin.z;
  return sqrt(xd*xd + zd*zd);
}

//Vector3d Circle::distance_vector(Lvec sh) {
  //// return the distance vector from the circle surface (sh should be projected to xz)
  //Vector3d dv = -sh.distance_vector(origin);
  //double dvl = dv.len();
  //dv.scale( (dvl-r)/dvl );
  //return dv;
//}

double Circle::on_surface(Vector3d pt) 
{
  return (distance(pt) - r) < float_tol; 
}

double Circle::calc_theta(Vector3d pt) {
  assert(this->on_surface(pt)); // paranoid check
  Vector3d topt = pt - origin;
  return atan2(topt.z, topt.x);
}

//https://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm

boost::optional<Vector3d> Circle::intersects(Lvec& shxyz)
{
  // to xyplane
  Lvec sh = shxyz.xzproject();
  // but use shxyz.get_pt() for return value

  Vector3d d = sh.line;
  Vector3d f = sh.source- this->origin;
  double a = d.dot(d);
  double b = 2*f.dot(d);
  double c = f.dot(f) - r*r;
  double discr = b*b - 4*a*c;
  if (discr < 0)
  {
    // Miss
    return boost::none;
  }
  else
  {
    double sqdet = sqrt(discr);
    double t1 = (-b - sqdet)/(2*a);
    double t2 = (-b + sqdet)/(2*a);

    if ( t1 >= 0 && t1 <= 1 )
    {
      // Impale // Poke
      return shxyz.get_pt(t1);
    }

    if ( t2 >= 0 && t2 <= 1 )
    {
      // ExitWound. 
      return shxyz.get_pt(t2);
    }
  }
  // FallShort // Past // CompletelyInside
  return boost::none;
}

//////////////////////////////////////////////////////////////////////////////////
// XZarc

boost::optional<Vector3d> XZarc::intersects(Lvec& sh)
{
  boost::optional<Vector3d> inter = Circle::intersects(sh);
  if (!inter) { // No intersection, just pass it on.
    return inter;   
  }
  else { // Check intersection is on arc
    double theta = this->calc_theta(*inter);
    if (theta >= thmin && theta <= thmax) {
      return inter;
    }
    else {
      // Found and intersection with the circle but not the arc
      return boost::none;
    }
  }
  // never gets here
}

///////////////////////////////////////////////////////////////////
// InfCylinder


boost::optional<Vector3d> InfCylinder::intersects(Lvec& lv) {
  return arc.intersects(lv); 
}


bool InfCylinder::_range_not_satisfied(double theta) {
  return (theta < arc.get_thmin() || theta > arc.get_thmax());
}

vector<overlap> InfCylinder::overlap_vector(Capsule& body) {
  vector<overlap> ret; // return value
  if (!solid) {
    return ret;
  }
  Vector3d origin = get_origin();
  // check special case of two contacts here! 
  // check this case works TODO
  // current implementation does a lot of work to check for this case, catch it earlier?

  /* off
  Vector3d hdv  = body.get_headpt().xzproject()-origin;
  Vector3d tdv  = body.get_endpt().xzproject()-origin;
  if (hdv.len()-arc.r < body.R && tdv.len()-arc.r < body.R) {
    double m = body.length/2.;
    // head and tail contacts
    // should check ranges?
    if (!_range_not_satisfied(std::atan2(hdv.z, hdv.x)) 
        && !_range_not_satisfied(std::atan2(tdv.z, tdv.x)) )
    {
      //cout << "two contacts with cylinder special case" << endl;
      ret.push_back(std::make_tuple(m, hdv.len()-arc.r, hdv.unit()));
      ret.push_back(std::make_tuple(-m, tdv.len()-arc.r, tdv.unit()));
      return ret; // no more contacts permitted
    }
  }
  */

  // otherwise we have only have one possible contact 
  Lvec bodylv = body.get_lvec();
  double t = bodylv.xzproject().closest_t(get_origin());  
  Vector3d lpt = bodylv.get_pt(t); // this should be the closest pt on the body axis segment
  Vector3d rv = Vector3d(lpt.x - origin.x, 0, lpt.z - origin.z); // origin - lpt

  // check this overlap is compatible with the arc
  // TODO condition may be too strict
  if (_range_not_satisfied(std::atan2(rv.z, rv.x))) {
    return ret; // no contact
  }
  double rd = rv.len();

  if (rd-arc.r < body.R) { // intersection
    double m = (0.5-t) * body.length; // body.get_lvec() is head -> tail
    ret.push_back(std::make_tuple(m, rd - arc.r, rv.unit())); // check this
  }
  return ret; // can be empty vector
}


const Matrix3d dexzde = Matrix3d(1,0,0,0,0,0,0,0,1);

Eigen::RowVectorXf InfCylinder::surface_grad_part(Capsule& body, contactf contact, overlap over) 
{
  vector<double> pk = {0,0,0};
  Rk bodyRk = body.get_Rk();
  double m = get<0>(over);
  double r = get<1>(over);
  Vector3d rhat = get<2>(over);

  double fg = -contact(r);
  Vector3d fgrad = fg*rhat;
  for (int k=0; k<3; k++) {
      pk[k] = fg * m * dot(rhat, dexzde * bodyRk.Rkvez[k]);
  }
  Eigen::RowVectorXf grad(6);
  grad << fgrad.x, fgrad.y, fgrad.z, pk[0], pk[1], pk[2];
  return grad;
}


Eigen::RowVectorXf InfCylinder::surface_grad(Capsule& body, contactf contact) {
  //
  vector<overlap> inter = this->overlap_vector(body); // overlaps are in xz plane
  if (inter.size() == 0)
  {
    return Eigen::RowVectorXf::Zero(6);
  }

  Vector3d fgrad;
  double fg;
  vector<double> pk = {0,0,0};
  
  Rk bodyRk = body.get_Rk();

  double m, r;
  Vector3d rhat;
  for (overlap over : inter) {
    // m is computed here in 3d. We project \ehat and m remains the same
    m = get<0>(over); 
    r = get<1>(over);
    rhat = get<2>(over);
    fg = -contact(r);
    fgrad += fg * rhat;
    for (int k=0; k<3; k++) {
      pk[k] += fg * m * dot(rhat, dexzde * bodyRk.Rkvez[k]);
    }
  }

  Eigen::RowVectorXf grad;
  grad << fgrad.x, fgrad.y, fgrad.z, pk[0], pk[1], pk[2];
  return grad;
}


double InfCylinder::operator()(double x) {
  assert(is_lower_hemisphere || is_upper_hemisphere);
  double xp = x - get_origin().x;
  double sgn = (is_lower_hemisphere) ? -1 : 1;
  double r = get_r();
  return get_origin().z + sgn * sqrt(r*r - xp*xp);
}

std::string InfCylinder::__str__() {
  std::string form = "Cylindrical: origin%s radius(%12.6f) arc(%12.6f -> %12.6f)";
  return ( boost::format(form) % get_origin().__str__() % get_r() % get_thmin() % get_thmax() ).str();
}

std::ostream& operator<<(std::ostream& os, InfCylinder& cyl)
{
  os << cyl.__str__();
  return os;
}


//////////////////////////////////////////////////////////////////////////////////
// Step 


float InfSteps::_get_n(double X) 
{
  return (X/_sep) + 0.5;
}

double InfSteps::_get_X(float n)
{
  return (n-0.5)*_sep;
}

double InfSteps::_get_X(int n)
{
  return _get_X((float)n);
}

bool InfSteps::is_lower(int n)
{
  return ((n % 2) == 0);
}

bool InfSteps::is_close_prev(float n)
{
  return n - floor(n) < _close_n;
}

bool InfSteps::is_close_next(float n)
{
  return ceil(n) - n < _close_n;
}


std::shared_ptr<PartialPlaneX> InfSteps::_get_partial_plane(int n) 
{
  if (!_xplanes.count(n)) {
    double xmin = _get_X(n) + _smallr;
    double xmax = _get_X(n+1) - _smallr;
    double height = (n % 2 == 0) ? 0. : _height;
    Vector3d origin{(xmax+xmin)/2., 0, height};
    auto pp = std::make_shared<PartialPlaneX>(origin, e_z, xmin, xmax);
    _xplanes[n] = pp; // shared_ptr is moved or copied?
  }
  return _xplanes[n];
}

std::shared_ptr<PartialPlaneZ> InfSteps::_get_partial_vertical(int n) 
{
  if (!_zplanes.count(n)) {
    double x = _get_X(n);
    Vector3d xnormal = (n % 2 == 0) ? e_x : -1*e_x;
    auto pp = std::make_shared<PartialPlaneZ>(Vector3d(x,0,(h2+h1)/2), xnormal, h1, h2);
    _zplanes[n] = pp;
  }
  return _zplanes[n];
}

std::shared_ptr<InfCylinder> InfSteps::_get_prev_corner_at(int n)
{
  if (!_pcorners.count(n)) {
    double x = _get_X(n);
    std::shared_ptr<InfCylinder> pp;
    if (is_lower(n)) {
      Vector3d origin{x+_smallr, 0, _smallr};
      pp = std::make_shared<InfCylinder>(origin, _smallr, -M_PI, -M_PI/2., false);
    }
    else { 
      Vector3d origin{x+_smallr, 0, h2};
      pp = std::make_shared<InfCylinder>(origin, _smallr, M_PI/2., M_PI, true);
    }
    _pcorners[n] = pp;
  }
  return _pcorners[n];
}

std::shared_ptr<InfCylinder> InfSteps::_get_next_corner_at(int n)
{
  if (!_ncorners.count(n)) {
    double xp = _get_X(n+1);
    std::shared_ptr<InfCylinder> pp;
    if (is_lower(n)) {
      Vector3d origin{xp-_smallr, 0, _smallr};
      pp = std::make_shared<InfCylinder>(origin, _smallr, -M_PI/2., 0., false); 
    }
    else { 
      Vector3d origin{xp-_smallr, 0, h2};
      pp = std::make_shared<InfCylinder>(origin, _smallr, 0., M_PI/2., true); 
    }
    _ncorners[n] = pp;
  }
  return _ncorners[n];
}

vector<std::shared_ptr<Plane>> InfSteps::_get_forwards_planes(int n) {
  vector<std::shared_ptr<Plane>> geometry;
  geometry.push_back(_get_partial_plane(n));
  geometry.push_back(_get_partial_vertical(n+1));
  geometry.push_back(_get_partial_plane(n+1));
  return geometry;
}

vector<std::shared_ptr<Plane>> InfSteps::_get_forwards_geometry(int n) {
  vector<std::shared_ptr<Plane>> geometry;
  geometry.push_back(_get_partial_plane(n));
  geometry.push_back(_get_next_corner_at(n));
  geometry.push_back(_get_partial_vertical(n+1));
  geometry.push_back(_get_prev_corner_at(n+1));
  geometry.push_back(_get_partial_plane(n+1));
  return geometry;
}

vector<std::shared_ptr<Plane>> InfSteps::_get_backwards_geometry(int n) {
  vector<std::shared_ptr<Plane>> geometry;
  geometry.push_back(_get_partial_plane(n));
  geometry.push_back(_get_prev_corner_at(n));
  geometry.push_back(_get_partial_vertical(n));
  geometry.push_back(_get_next_corner_at(n-1));
  geometry.push_back(_get_partial_plane(n-1));
  return geometry;
}


boost::optional<Vector3d> InfSteps::intersects(Lvec& lv) 
{
  // get 
  double x1 = lv.get_origin().x;
  double x2 = lv.get_endpt().x;
  float n1 = _get_n(x1);
  float n2 = _get_n(x2);
  int n = (int)floor(n1);
  
  if (floor(n1) == floor(n2)) {
    // segment is contained on one planar element
    // order n1, n2 
    if (n1 > n2) { std::swap(n1, n2); }
    auto inter = _get_partial_plane(n)->intersects(lv);
    if (inter) { 
      return  inter;
    }
    // check the previous corner
    if (is_close_prev(n1)) {
      auto inter = _get_prev_corner_at(n)->intersects(lv);
      if (inter) {
        return inter;
      }
    }
    // check the next corner
    if (is_close_next(n2)) {
      auto inter = _get_next_corner_at(n)->intersects(lv);
      if (inter) {
        return inter;
      }
    }
    // finally 
    return boost::none;
  }
  else  
  {
    // need to get all the elements from n1 -> n2 
    bool is_forwards = (n2 > n1);

    vector<std::shared_ptr<Plane>> step;
    if (is_forwards) {
      // partial plane, next corner, vertical (n2), prev_corner_n2, partial plane n2
      step = _get_forwards_geometry(n);
    }
    else {
      step = _get_backwards_geometry(n);
    }
    for (std::shared_ptr<Plane> surface_part : step ) {
      auto inter = surface_part->intersects(lv);
      if (inter) {
        return inter;
      }
    }
    return boost::none;
  }
  // never gets here
}

vector<overlap> InfSteps::overlap_vector(Capsule& body) {
  if (!_now_touching) {
    _now_touching = std::vector<std::shared_ptr<Plane>>();
  }
  _now_touching->clear(); // this method is responsible for updating the touching list with anything we are touching
  vector<overlap> ret; 
  Lvec bodyax = body.get_lvec();
  double x1 = bodyax.get_origin().x; 
  double x2 = bodyax.get_endpt().x; 
  // sort x1, x2
  if (x1 > x2) {
    std::swap(x1, x2);
  }
  // add bacteria radius
  double ext1 = x1 - body.R;
  double ext2 = x2 + body.R;
  float n1 = _get_n(ext1);
  float n2 = _get_n(ext2);
  int n = (int)floor(n1);
  if (floor(n1) == floor(n2)) {
    // not touching sides, not touching corners (can't touch bottom corner)
    std::shared_ptr<Plane> plane = _get_partial_plane(n);
    ret = plane->overlap_vector(body); // [0,2] overlaps
    for (int i=0; i<ret.size(); i++) {
      _now_touching->push_back(plane);
    }
    // not finished. If we are on the upper surface and we are not touching
    // a plane then we may be touching a corner
    if (ret.empty()) 
    {
      // let these options be mutually exclusive 
      std::shared_ptr<InfCylinder> corner;
      if (is_close_prev(n1)) {
        corner = _get_prev_corner_at(n);
        ret = corner->overlap_vector(body);
      }
      else if (is_close_next(n2)) 
      {
        corner = _get_next_corner_at(n);
        ret = corner->overlap_vector(body);
      }
      if (!ret.empty()) { 
        _now_touching->push_back(corner); 
      }
    }
    return ret;
  }
  else
  {
    vector<std::shared_ptr<Plane>> step;
    // can touch different geometry objects each having different gradient calculations
    // 1. check if we touch the upper corner

    /*
      Get all the possible contacts and then prune `duplicates`
     */

    std::shared_ptr<Plane> corner;
    if (is_lower(n)) {
      corner = _get_prev_corner_at(n+1); 
    }
    else {
      corner = _get_next_corner_at(n);
    }

    vector<overlap> corner_touch = corner->overlap_vector(body);
    for (int i = 0; i < corner_touch.size(); i++) {
      ret.push_back(corner_touch[i]);
      _now_touching->push_back(corner);
    }

    // 2. deal with planes
    step = _get_forwards_planes(n); 
    int z_touching = 0;
    for (std::shared_ptr<Plane> surface_part : step) {
      vector<overlap> inter = surface_part->overlap_vector(body); 
      for (int i = 0; i <inter.size(); i++) {
        overlap over = inter[i];
        ret.push_back(over);
        _now_touching->push_back(surface_part);
      }
    }

    return ret; 
  }
  // never gets here 
}


int InfSteps::get_num_contacts(Capsule& body) {
  // fast implementation doesn't call overlap_vector an additional time 
  // but we must have called overlap_vector or grad before calling this!
  if (!_now_touching) {
    this->overlap_vector(body); // update the touching list
  }
  return _now_touching->size();
}

Eigen::RowVectorXf InfSteps::surface_grad(Capsule& body, contactf contact)  {
  // must call overlap_vector to update _now_touching
  vector<overlap> inter = this->overlap_vector(body); 
  Eigen::RowVectorXf grad = Eigen::RowVectorXf::Zero(6);

  // accumulate gradient contributions from all touching objects
  for (int i = 0; i < inter.size(); i++) {
    overlap over = inter[i];
    std::shared_ptr<Plane> surface_part = _now_touching->at(i);
    // calling surface_grad calls overlap vector
    // need to separate gradient calculation and overlap calculation
    Eigen::RowVectorXf grad_part = surface_part->surface_grad_part(body, contact, over);

    //cout << "touching " << surface_part->__str__() << endl;
    //cout << over << endl;
    //cout << "Add surface gradient conribution " << grad_part << endl;
    grad += grad_part;
  }

  return grad;
}



double InfSteps::operator()(double x) {
  float nvar = _get_n(x);
  int n = (int)floor(nvar);
  if (is_close_prev(nvar)) {
    return _get_prev_corner_at(n)->operator()(x);
  }
  else if (is_close_next(nvar)) {
    return _get_next_corner_at(n)->operator()(x);
  }
  else {
    if (is_lower(n)) {
      return 0.;
    }
    else {
      return _height;
    }
  }
  // never gets here
}


std::string InfSteps::report_touching() {
  if (!_now_touching) {
    return "No touching data.";
  }
  std::string report = "--Touching List Start\n";
  for (auto surface_part : _now_touching.get() ) {
    report += surface_part->__str__();
    report += "\n";
  }
  report += "--End";
  return report;
}

std::string InfSteps::__str__() {
  std::string form = "InfSteps: height(%f) separation(%f) smallr(%f)";
  return ( boost::format(form) % _height % _sep % _smallr ).str();
}

//const std::string InfSteps::__class__ = "InfSteps";
