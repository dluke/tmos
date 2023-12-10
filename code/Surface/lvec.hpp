

#ifndef __LVEC_HPP__
#define __LVEC_HPP__

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <Eigen/Dense>

#include "sample.hpp"
#include "frame.hpp"
#include "vector3d.hpp"
#include "matrix3d.hpp"

#include <vector>
using std::vector;

using std::abs;

// LineBase
// Base class for all things Line-like

// line segments are not shapes, but they will have 
//similar methods: intersects, translate, rotate, intersects should return a pt or None
//
// Line segments can be represented as a point plus a vector

class Lvec
{

  public:
  Vector3d source, line;
  // attachment tracking
  double inter_t;
  Lvec(Vector3d source, Vector3d line) : source(source), line(line) {}

  void translate(Vector3d r) {
    this->source += r;
  }

  void rotate(Matrix3d M) {
    this->line = M * line;
  }

  Vector3d get_origin() { return source; }
  Vector3d get_endpt() { return source + line; }
  Vector3d get_pt(double t) { return source + t * line; }

  // Just want to know if the point is on the positive (anti-clockwise) side of the line 
  int check_side(Vector3d pt);

  // interpolate the line with N+1 points (include endpoints)
  vector<Vector3d> interp(int N);

  // return the t coordinate for projecting pt onto the line defined by this segment
  double project_t(Vector3d pt);
  // return the point determined by coordinate t, if t not in [0,1] use closest value
  double contract_t(double t);
  // same as contract_t(project_t)
  double closest_t(Vector3d pt);
  // same as this->operator()(closest_t(pt)) Vector3d closest_pt(Vector3d pt);

  // closest pt on this segment to pt
  Vector3d closest_pt(Vector3d pt);
  // the shortest vector from the line to the point
  Vector3d distance_vector(Vector3d pt);

  // Length
  double len_to_inter() { return inter_t*len(); }
  double len() { return line.len(); }
  // The minmum height of line segment above the surface
  double minz() { return std::min(source.z, get_endpt().z); }

  // convenience method for projecting onto xyplane through origin
  Lvec xyproject() 
  { return Lvec(source.xyproject(), line.xyproject()); }
  Lvec xzproject() 
  { return Lvec(source.xzproject(), line.xzproject()); }

  std::string __str__();

};

std::ostream& operator<<(std::ostream& os, Lvec& lv);

Lvec lvec_from_pts(Vector3d, Vector3d);


#endif

