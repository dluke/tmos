
#include <iostream>
#include <exception>
#include <boost/format.hpp>

#include "sphere.hpp"
#include "matrix3d.hpp"

using std::endl;
using std::cout;

const double pi = M_PI;


/////

bool tolequals(double x, double y) { return (std::abs(x-y) < Base::FLOAT_TOL); }

bool Sphere::contains(Vector3d pt) {
  return ((pt - frame.origin).len() < R);
}

bool Sphere::on_surface(Vector3d pt) {
  return tolequals( (pt - frame.origin).len(), this->R );
}

Vector3d Sphere::radv(Vector3d pt) {
  return pt - frame.origin;
}

Vector3d Sphere::normal(Vector3d pt) {
  return radv(pt).unit();
}

Vector3d Sphere::spherical(double theta, double phi) {
  return this->frame.spherical(this->R, theta, phi);
}

// todo write test 
//https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
boost::optional<Vector3d> Sphere::intersects(Lvec& lv) 
{
  // unit() fails when a segment of the pilus has 0 length 
  // I thought this shouldn't happen but since it does we can fix 
  Vector3d lvv = lv.line.prot_unit();
  Vector3d lsc = lv.source - this->get_origin(); // vector from sphere centre to line source
  double b = lvv.dot(lsc);
  double sqdet = b*b- lsc.dot(lsc) + R*R;
  if (sqdet < 0)
    return boost::none; // no intersection
  double d1 = -b - sqrt(sqdet);
  double d2 = -b + sqrt(sqdet);
  // check roots are on line segment
  double ld = lv.len();
  bool d1valid = ((d1 >= 0.) && (d1 <= ld));
  if (d1valid)
    return lv.get_pt(d1/ld);
  bool d2valid = ((d2 >= 0.) && (d2 <= ld));
  if (d2valid)
    return lv.get_pt(d2/ld);
  if (!d1valid && !d2valid) {
    return boost::none; // intersection of line but not of this line segment
  }
  else 
    { throw std::runtime_error("Failed to determine line vector -> sphere intersection"); } 
}

vector<overlap> Sphere::overlap_vector(Capsule& body)
{
  boost::optional<overlap> inter = body.intersects(*this);
  if (inter) {
    return {*inter};
  }
  return {};
}

std::string Sphere::__str__() {
  std::string form = "Sphere: origin%s Radius(%8.2f)";
  return ( boost::format(form) % frame.origin.__str__() % R ).str();
}

//////////////////////////////////////////////////////////////////////////////////

Vector3d Shell::normal(Vector3d pt) {
  return (get_origin() - pt).unit();
}

Vector3d Shell::distance_vector(Vector3d pt) 
{
  return (R/pt.len() - 1) * -pt;
}

vector<overlap> Shell::overlap_vector(Capsule& body)
{
  vector<overlap> ret;
  // check head and tail
  Vector3d hpt = body.get_headpt();
  Vector3d tpt = body.get_endpt();

  Vector3d hdv = distance_vector(hpt);
  Vector3d tdv = distance_vector(tpt);

  double d;
  d = hdv.len();
  if (d < body.R) {
    overlap over = std::make_tuple(-body.length/2., d, hdv/d);
    ret.push_back(over);
  }
  
  d = tdv.len();
  if (d < body.R) {
    overlap over = std::make_tuple(body.length/2., d, tdv/d);
    ret.push_back(over);
  }
  return ret;
}

/////////////////////////////////////////////////////////////////////////////////

int Capsule::HEAD = 0;
int Capsule::BODY = 1;
int Capsule::TAIL = 2;

// because of the way we defined the body frame with e_z along body axis
// This method returns body frame vector
Vector3d Capsule::alt_spherical(double theta, double phi)
{
  // theta and phi are in the default x,y,z frame
  // Need to rotate the result onto the capsule frame
  //Matrix3d M = this->frame.get_tmatrix();

  // construct spherical coordinate in body frame 
  Vector3d bodya = Vector3d(
      this->R * sin(theta) * cos(phi),
      this->R * sin(theta) * sin(phi),
      this->R * cos(theta)
      );

  return bodya;
}

double Capsule::long_length() { return this->length + pi * R; }
double Capsule::true_length() { return this->length + 2 * R; }

Vector3d mirroru(Vector3d v, const Vector3d u) {
  Vector3d vdot = v.dot(u) * u;
  return v - 2*vdot;
}
// mirror along long axis (in the body frame of course)
Vector3d Capsule::mirror(Vector3d v) {
  return mirroru(v, this->get_axis());
}


int Capsule::get_segment(double d)
{
  double sq = (pi * R)/2.;
  if (d <= sq) {
    return Capsule::HEAD;
  }
  else if (sq <= d && d < (sq+length)) {
    return Capsule::BODY;
  }
  else if (sq+length <= d && d <= long_length()) {
    // theta around the opposite pole
    return Capsule::TAIL;
  }
  else { 
    throw std::runtime_error(
        (boost::format("Surface length %f larger than the long surface length %f") 
        % d % long_length()
        ).str()
          ); 
  }
}

std::pair<Vector3d, Vector3d> Capsule::surface_pt(double d, double phi, double inlength)
{
  //IMPORTANT we designate a length inside the cell to anchor the pili
  double infactor = (R - inlength)/R;

  double sq = (pi * R)/2.;
  int segment = get_segment(d);
  //cout << "segment " << segment << endl;
  std::pair<Vector3d,Vector3d> spt_unitax;
  if (segment == Capsule::HEAD) {
    // 2 * pi * ( d / 2 * pi * R)
    Vector3d rad = infactor * alt_spherical(d / R, phi);
    Vector3d pt = Vector3d(0,0,this->length/2.) + rad; 
    spt_unitax = std::make_pair(pt, rad.unit());
  }
  else if (segment == Capsule::BODY) {
    Vector3d rad = alt_spherical(pi/2., phi);
    double m = 1. - (d-sq);
    Vector3d pt = infactor * rad + m*e_z;
    spt_unitax = std::make_pair(pt, rad.unit());
  }
  else if (segment == Capsule::TAIL) {
    // theta around the opposite pole
    double r_theta = pi/2. - (d - length - sq)/R;

    Vector3d alt = mirroru( alt_spherical(-r_theta, phi), e_z );  // choosing -r_theta, phi
    //Vector3d body_endpt = length * -1 * e_z;
    Vector3d pt =  Vector3d(0,0,-this->length/2.) + infactor * alt;
    spt_unitax = std::make_pair(pt, alt.unit());
  }
  else { 
    throw std::runtime_error(
        (boost::format("Segment %d is not any of %d,%d,%d") 
        % segment % Capsule::HEAD % Capsule::BODY % Capsule::TAIL
        ).str()
          ); 
  }
  return spt_unitax;
}

Vector3d Capsule::spt(double d, double phi)
{
  return this->surface_pt(d, phi, 0.).first;
}

// can be furthur simplified, see cylinder.cpp
boost::optional<overlap> Capsule::intersects(Sphere& sph) {
  Lvec capslv = this->get_lvec();
  Vector3d pt = sph.get_origin();
  Vector3d clpt = capslv.closest_pt(pt);
  
  // vector from centerpoint to closest pt
  Vector3d mv = (clpt - this->get_centerpt());

  double msign = (mv.dot(this->get_axis()) >= 0) ? 1 : -1;
  double m = msign * mv.len();

  // vector from the sphere to the capsule segment
  Vector3d dv = clpt - pt;
  double dist = dv.len();
  double sumR = this->R + sph.R;
  if (dist >= sumR) {
    return boost::none;
  } else {
    return std::make_tuple(
        m,
        //dist - this->R, // how can this be wrong?
        dist - sph.R, 
        dv.unit()
        );
  }
}

void Capsule::set_rotation(Vector3d p) {
  this->frame.set_rotation(p);
}

// slow? todo
Vector3d Capsule::get_axisangle() {
  Matrix3d Rp{frame.e1, frame.e2, frame.e3};
  Vector3d p = axisangle(Rp);
  return p;
}

// return a line vector of the caspule axis
Lvec Capsule::get_lvec() {
  return Lvec(get_headpt(), get_vector());
}
// normal at pt on capsule surface
Vector3d Capsule::normal(Vector3d pt) {
  return get_lvec().distance_vector(pt).unit();
}
// assume pt is on the head sphere then return normal
Vector3d Capsule::head_normal(Vector3d pt) {
  return pt - frame.origin;
}
/* */

double Capsule::sgn_distance_to_surface(Vector3d pt) {
  Lvec bodyax{Vector3d(0,0,length/2.),Vector3d(0,0,-length)};
  //cout << bodyax << endl;
  //cout << bodyax.distance_vector(pt) << endl;
  return bodyax.distance_vector(pt).len() - this->R; 
}

std::string Capsule::__str__() {
  std::string form = "Capsule: radius(%4.2f) length(%4.2f)\n->%s";
  return ( boost::format(form) % R % length % frame.__str__() ).str();
}


std::ostream& operator<<(std::ostream& os, Capsule& cap)
{
  os << cap.__str__();
  return os;
}


namespace Caps {
  Capsule headat(Vector3d origin, Vector3d orient) {
    return Capsule(origin + -1. * orient, orient);
  }
  
  Capsule tailat(Vector3d origin, Vector3d orient) {
    return Capsule(origin + orient, orient);
  }

  Capsule body_from_lv(Lvec lv) {
    return Capsule(lv.get_pt(0.5), -lv.line.unit());
  }
}


////////////////////////////////////////////////////////////////////////////////



