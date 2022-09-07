
#ifndef __FRAME_HPP__
#define __FRAME_HPP__

#include <cmath>
#include <array>
#include <string>
#include <tuple>
#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "axisangle.hpp"


// Coordinate frames are general objects. Here we only consider stationary
// frames with mutually perpendicular basis vectors.

class Frame
{
  
  public:
  Vector3d origin, e1, e2, e3;

  Frame() :
    origin(Vector3d()), e1(Vector3d(1.,0.,0.)), e2(Vector3d(0.,1.,0.)), e3(Vector3d(0.,0.,1.)) {}


  // Frame with position only
  Frame(Vector3d r) : origin(r) 
  {
    e1 = Vector3d(1.,0.,0.);
    e2 = Vector3d(0.,1.,0.);
    e3 = Vector3d(0.,0.,1.);
  }

  // Frame from rotation matrix 
  Frame(Vector3d r, Matrix3d R) : origin(r) 
  {
    e1 = R * e_x;
    e2 = R * e_y;
    e3 = R * e_z;
  }

  // Fully defined frame
  Frame(Vector3d origin, Vector3d e1, Vector3d e2) :
    origin(origin), e1(e1), e2(e2) 
  {
    e3 = e1.cross(e2);
  }

  // All vectors defined
  Frame(Vector3d origin, Vector3d e1, Vector3d e2, Vector3d e3) :
    origin(origin), e1(e1), e2(e2), e3(e3) {}

  // frame with normal vector only defined
  // so this defines the e3 vector!
  // e1 and e2 vectors are chosen for a right handed coordinate system
  // so if normal = e3 = Vector3d(0.,0.,1.) = e_z
  // then e1 = ex, e2 = ey
  Frame(Vector3d origin, Vector3d normal) : origin(origin), e3(normal) {
    // ry90 sends e_z to e_x
    e1 = perp(normal);
    e2 = e3.cross(e1); // left handed coordinate system
  }

  Vector3d get_origin() { return origin; }


  void translate(Vector3d r) {
    this->origin += r;
  }

  void rotate(Matrix3d M) {
    e1 = M * e1;
    e2 = M * e2;
    e3 = e1.cross(e2);
  }


  // set the orientation of the frame with a matrix M such that M*e_z = e3
  void set_rotation(Matrix3d M) {
    e3 = M * e_z;
    e1 = M * e_x;
    e2 = e3.cross(e1);
  }
  void set_rotation(Vector3d p) {
    this->set_rotation(aarmatrix(p));
  }

  std::array<double, 3> orthogonal_error();
  std::array<double, 3> unit_error();

  void cross_norm();
  double orthogonalise();
  double normalise();

  Matrix3d get_rmatrix();
  Matrix3d get_tmatrix();
  Vector3d to_lab_rt(Vector3d v);
  Vector3d to_lab(Vector3d v);

  // obtain the spherical polar coordinate pt in this frame for a radius R
  static Vector3d spherical(double R, double theta, double phi)
  {
    double x = sin(theta) * cos(phi);
    double y = sin(theta) * sin(phi);
    double z = cos(theta);
    return R * Vector3d(x,y,z);
  }


  std::string __str__();

};

std::ostream& operator<<(std::ostream& os, Frame& frame);

// if this frame is constructed out of cartesian basis vectors then return 
// the indices of the those vectors like {0:e_x, 1:e_y, 2:e_z}

std::tuple<int,int,int> cartesian_basis_vectors(Frame& frame);

#endif
