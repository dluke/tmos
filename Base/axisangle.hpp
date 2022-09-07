
#ifndef __AXISANGLE_HPP__
#define __AXISANGLE_HPP__

#include <cmath>
#include <vector>
#include "vector3d.hpp"
#include "matrix3d.hpp"

using std::cos;
using std::sin;
using std::vector;

// get rotation matrix from axisangle
Matrix3d aarmatrix(Vector3d p);

// get axisangle from rotation matrix
Vector3d axisangle(Matrix3d M);

Vector3d perp(Vector3d v);

// compute the deriviatives of the rotation matrix R(p) 
class Rk {
  public:

  // calculate the derivatives of the rotation matrix related to axis angle vector p
  Rk(Vector3d p);
  Rk() : Rk(e_z) {}

  Matrix3d Rk1, Rk2, Rk3;

  Matrix3d Rkget(int k) { return Rkv[k]; }

  vector<Matrix3d> Rkv;
  vector<Vector3d> Rkvez;
};

extern Matrix3d skewx;
extern Matrix3d skewy;
extern Matrix3d skewz;
extern vector<Matrix3d> skewk;


#endif
