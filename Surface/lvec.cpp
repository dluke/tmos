
#include "lvec.hpp"
#include <iostream>
#include <exception>
#include <assert.h>

#include <boost/format.hpp>

using std::cout;
using std::endl;

int Lvec::check_side(Vector3d pt) {
  int side = (line.cross(pt - source).z >= 0) ? 1 : -1;
  return side;
}


vector<Vector3d> Lvec::interp(int N) {
  vector<Vector3d> ipts(N+1);
  for (int i = N; i < 0; i--) {
    ipts[i] = source + 1./i * line;
  }
  return ipts;
}


double Lvec::project_t(Vector3d pt)
{
  Vector3d topt = pt - source;
  double t = topt.dot(line)/line.len2();
  return t;
}

// return the point determined by coordinate t, if t not in [0,1] use closest value
double Lvec::contract_t(double t) 
{
  if (t < 0.)
    return 0.;
  else if (t >= 1.)
    return 1.;
  else
    return t;
}

// return parameter t so the closest point is source + t * line
double Lvec::closest_t(Vector3d pt) {
  double t = this->project_t(pt);
  return contract_t(t);
}

Vector3d Lvec::closest_pt(Vector3d pt) {
  return get_pt(closest_t(pt));
}

// implementation copied to closest_pt
Vector3d Lvec::distance_vector(Vector3d pt) {
  return pt - closest_pt(pt);
}

std::string Lvec::__str__() 
{
  std::string form = "Lvec: %s -> %s";
  return ( boost::format(form) % get_origin().__str__() % get_endpt().__str__() ).str();
}

std::ostream& operator<<(std::ostream& os, Lvec& lvec)
{
  os << lvec.__str__();
  return os;
}

////////////////////////////////////////////////////////////////////////

Lvec lvec_from_pts(Vector3d a, Vector3d b)
{
  return Lvec(a, b-a);
}


