
#include "axisangle.hpp"

#include <stdexcept>
#include <iostream>

using std::cout;
using std::endl;

using std::cos;
using std::sin;

Matrix3d skewx = skew_matrix(e_x);
Matrix3d skewy = skew_matrix(e_y);
Matrix3d skewz = skew_matrix(e_z);

vector<Matrix3d> skewk = {skewx, skewy, skewz};


const double theta_tol = 0.000001;

// get rotation matrix from axisangle
Matrix3d aarmatrix(Vector3d p) {
  double theta = p.len();
  // theta can be zero then we need to exit early to avoid NAN
  if (theta == 0.) 
    return Matrix3d::I();
  Matrix3d pskew = skew_matrix(p);
  return Matrix3d::I() + (1 - cos(theta))/(theta*theta) * pskew*pskew
   + sin(theta)/theta * pskew;
}

// get an arbitrary vector perpendicular to v
Vector3d perp(Vector3d v) {
  if ((v.x==0.) && (v.y==0)) {
    if (v.z==0.) {
      throw 
        std::invalid_argument("cannot find vector perpendicular to zero vector");
    }
    // if v is close to e_z/-e_z then return e_x
    return Vector3d(1,0,0);
  }
  return Vector3d(-v.y, v.x, 0).unit();
}

// convert rotation matrix to axisangle
Vector3d axisangle(Matrix3d M) {
  // get the axisangle from a rotation matrix
  // what if M is not a normalised rotation matrix?

  double ctheta = ((M.trace() -1)/2);
  // sanity check
  if (abs(ctheta) > (1. + theta_tol) ) {
    cout << "warning ctheta = " << ctheta << endl;
  }
  ctheta = std::max(-1., std::min(ctheta, 1.)); //careful

  double theta;
  Vector3d phat;

  if (ctheta == 1) {
    return Vector3d();
  }
  else if (ctheta == -1) { //critical case 
    theta = M_PI;
    phat = (M*e_x + e_x).prot_unit(); // e_x is arbitrary choice // does this axis have the right sign? TODO
    return theta * phat;
  }
  else { 
  //https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis%E2%80%93angle
    Vector3d p = Vector3d(M.z2-M.y3, M.x3-M.z1, M.y1-M.x2);
    phat = p.prot_unit();
      //
    theta = std::acos(ctheta);
    // check sign of theta
    /*if (Mv.len() > 1. + theta_tol) {*/
      //cout << "Warning M is not length preserving, Mv.len() = " << Mv.len() << endl;
      //cout << "Mv " << Mv << endl;
    /*}*/
  }
  return theta * phat;
}


Rk::Rk(Vector3d p) { // check identical for pairs of p vector which differ in axis
  double theta = p.len();
  double thetasq = theta*theta;
  bool thetaiszero = (theta == 0.);
  if (!thetaiszero) {
    Matrix3d pskew1 = (1/theta) * skew_matrix(
        Vector3d(1-p.x*p.x/thetasq, -p.x*p.y/thetasq, -p.x*p.z/thetasq)
        );
    Matrix3d pskew2 = (1/theta) * skew_matrix(
        Vector3d(-p.x*p.y/thetasq, 1-p.y*p.y/thetasq, -p.y*p.z/thetasq)
        );
    Matrix3d pskew3 = (1/theta) * skew_matrix(
        Vector3d(-p.x*p.z/thetasq, -p.z*p.y/thetasq, 1-p.z*p.z/thetasq)
        );
    Matrix3d pskew = (1/theta) * skew_matrix(p);
    Rk1 = p.x*sin(theta)/theta * pskew*pskew
      + (1-cos(theta))*(pskew1*pskew + pskew*pskew1)
      + p.x*cos(theta)/theta * pskew
      + sin(theta)*pskew1;
    Rk2 = p.y*sin(theta)/theta * pskew*pskew
      + (1-cos(theta))*(pskew2*pskew + pskew*pskew2)
      + p.y*cos(theta)/theta * pskew
      + sin(theta)*pskew2;
    Rk3 = p.z*sin(theta)/theta * pskew*pskew
      + (1-cos(theta))*(pskew3*pskew + pskew*pskew3)
      + p.z*cos(theta)/theta * pskew
      + sin(theta)*pskew3;
  }
  else {
    // critical case theta == 0.
    Vector3d p1 = p;
    p1.x = 1.;
    Rk1 = tilde(p1);
    Vector3d p2 = p;
    p2.y = 1.;
    Rk2 = tilde(p2);
    Vector3d p3 = p;
    p3.z = 1.;
    Rk3 = tilde(p3);
  }
  Rkv.push_back(Rk1);
  Rkv.push_back(Rk2);
  Rkv.push_back(Rk3);
  Rkvez.push_back(Rk1*e_z);
  Rkvez.push_back(Rk2*e_z);
  Rkvez.push_back(Rk3*e_z);
}

