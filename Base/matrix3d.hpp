
#ifndef __MATRIX3D_HPP__
#define __MATRIX3D_HPP__

// having looked at the commonly used libraries they all seem 
// excessively complicated in architecture and somewhat lacking in 
// my specific requirements.
// It is difficult to naturally extend a library like Boost.geometry or CGAL
// without being a professional C++ programmer.
// Also I only need a very small and specific subset of features for this simulation.
//
// Ultimately it I feel it's preferable to personally implement Vectors, matrics, coordinate frames, 
//  line segments, sphere, etc.
// as well as the specfic surface and bacterial geometry of interest.
//
// This choice has the advantages of code simplicity and few dependencies.
// But may cause problems for implementing new features in the future.
// These free software libraries are also heavily optimised so I may lose speed.

#include <cmath>
#include <string>
#include "vector3d.hpp"

#include <Eigen/Dense>

// x1 x2 x3
// y1 y2 y3
// z1 z2 z3

class Matrix3d
{
  public:
  Matrix3d() : 
    x1(0.), x2(0.), x3(0.),
    y1(0.), y2(0.), y3(0.),
    z1(0.), z2(0.), z3(0.) {}

  Matrix3d(
      double x1, double x2, double x3,
      double y1, double y2, double y3,
      double z1, double z2, double z3)
    : x1(x1), x2(x2), x3(x3), y1(y1), y2(y2), y3(y3), z1(z1), z2(z2), z3(z3) 
  {}

  // Matrix from column vectors
  Matrix3d(Vector3d e1, Vector3d e2, Vector3d e3) 
    : Matrix3d(e1.x, e2.x, e3.x, e1.y, e2.y, e3.y, e1.z, e2.z, e3.z)
  {}

  

  // / x1 x2 x3 \ / v.x \
  // | y1 y2 y3 | | v.y |
  // \ z1 z2 z3 / \ v.z /
  Vector3d operator*(const Vector3d v) const
  {
    double x = x1 * v.x + x2 * v.y + x3 * v.z;
    double y = y1 * v.x + y2 * v.y + y3 * v.z;
    double z = z1 * v.x + z2 * v.y + z3 * v.z;
    return Vector3d(x,y,z);
  }

  Matrix3d operator-() const
  {
    return Matrix3d(-x1,-x2,-x3,-y1,-y2,-y3,-z1,-z2,-z3);
  }

  // transpose
  Matrix3d T() const
  {
    return Matrix3d(x1, y1, z1, x2, y2, z2, x3, y3, z3);
  }

  // Identity
  static Matrix3d I() 
  {
    return Matrix3d(1,0,0,0,1,0,0,0,1);
  }

  // Addition
  Matrix3d operator+(const Matrix3d m) const
  {
    return Matrix3d(
        x1+m.x1,
        x2+m.x2,
        x3+m.x3,
        y1+m.y1,
        y2+m.y2,
        y3+m.y3,
        z1+m.z1,
        z2+m.z2,
        z3+m.z3
        );
  }

  // / x1 x2 x3 \ / x1 x2 x3 \ 
  // | y1 y2 y3 | | y1 y2 y3 | 
  // \ z1 z2 z3 / \ z1 z2 z3 / 
  Matrix3d operator*(const Matrix3d m) const
  {
    return Matrix3d(
        x1*m.x1 + x2*m.y1 + x3*m.z1,
        x1*m.x2 + x2*m.y2 + x3*m.z2,
        x1*m.x3 + x2*m.y3 + x3*m.z3,

        y1*m.x1 + y2*m.y1 + y3*m.z1,
        y1*m.x2 + y2*m.y2 + y3*m.z2,
        y1*m.x3 + y2*m.y3 + y3*m.z3,

        z1*m.x1 + z2*m.y1 + z3*m.z1,
        z1*m.x2 + z2*m.y2 + z3*m.z2,
        z1*m.x3 + z2*m.y3 + z3*m.z3
        );
  }

  bool operator==(const Matrix3d m) const
  {
    return (x1 == m.x1 
        && x2 == m.x2
        && x3 == m.x3
        && y1 == m.y1
        && y2 == m.y2
        && y3 == m.y3
        && z1 == m.z1
        && z2 == m.z2
        && z3 == m.z3);
  }

  double trace(void) {
    return x1 + y2 + z3;
  }

  Eigen::MatrixXf to_eigen();

  Matrix3d svdnormalise();

  // Is using 9 variables like this a huge mistake?
  // We only need a few rotations of vectors.
  // stop worrying so much
  double x1, x2, x3, y1, y2, y3, z1, z2, z3;


  std::string __str__() const;

};

std::ostream& operator<<(std::ostream& os, Matrix3d m);

// non member functions

Matrix3d m3f_to_matrix3d(Eigen::MatrixXf m);

Matrix3d operator*(double const& c, Matrix3d m);
Matrix3d operator*(Matrix3d m, double const& c);
  
///////////////////////////////////////////////////////////
// Just some methods for constructing rotation matrices

// flip vector by exactly 90 degrees anticlockwise around e_x
const Matrix3d rx90 = Matrix3d(1.,0.,0.,0.,0.,-1.,0.,1.,0.);

const Matrix3d ry90 = Matrix3d(0.,0.,1.,0.,1.,0.,-1.,0.,0.);

const Matrix3d rx180 = Matrix3d(1.,0.,0.,0.,-1.,0.,0.,0.,-1.);

// construct matrix from a vector. This is not the only way to construct a 
// matrix from a vector but its the only one we need right now.
Matrix3d skew_matrix(Vector3d v);
// skew matrix but normalised
Matrix3d tilde(Vector3d v);

Matrix3d Rmatrix(Vector3d v, Vector3d u);

#endif

