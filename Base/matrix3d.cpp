
#include "matrix3d.hpp"

#include <iostream>
#include <boost/format.hpp>
using std::cout;
using std::endl;


std::string Matrix3d::__str__() const {
  return ( boost::format(
        "\n/%12.6f,%12.6f,%12.6f\\\n"
        "|%12.6f,%12.6f,%12.6f|\n"
        "\\%12.6f,%12.6f,%12.6f/"
        ) 
      % x1 % x2 % x3
      % y1 % y2 % y3
      % z1 % z2 % z3
      ).str();
}

std::ostream& operator<<(std::ostream& os, Matrix3d m)
{
  os << m.__str__();
  return os;
}


Matrix3d Matrix3d::svdnormalise() {
  Eigen::MatrixXf m = this->to_eigen();
  Eigen::JacobiSVD<Eigen::MatrixXf> svd(m, 
      Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXf muv = svd.matrixU() * svd.matrixV();
  return m3f_to_matrix3d(muv);
}

Eigen::MatrixXf Matrix3d::to_eigen() {
  Eigen::Matrix3f m;
  m << x1, x2, x3, y1, y2, y3, z1, z2, z3;
  return m;
}


Matrix3d operator*(double const& c, Matrix3d m)
{
  return Matrix3d(
      c*m.x1,
      c*m.x2,
      c*m.x3,
      c*m.y1,
      c*m.y2,
      c*m.y3,
      c*m.z1,
      c*m.z2,
      c*m.z3
      );
}
Matrix3d operator*(Matrix3d m, double const& c) {
  return c*m;
}

Matrix3d m3f_to_matrix3d(Eigen::MatrixXf m) {
  return Matrix3d(
      m(0,0),m(0,1),m(0,2),
      m(1,0),m(1,1),m(1,2),
      m(2,0),m(2,1),m(2,2)
      );
}

Matrix3d skew_matrix(Vector3d v) {
  return Matrix3d(0,-v.z,v.y,v.z,0,-v.x,-v.y,v.x,0);
}

Matrix3d tilde(Vector3d v) {
  return (1/v.len()) * skew_matrix(v);
}

// the matrix rotating v onto u
// need to make sure that v, u are unit vectors 
// May be numerically unstable for u close to -e_z
Matrix3d Rmatrix(Vector3d v, Vector3d u)
{
  Vector3d vu = v.cross(u);
  Matrix3d skm = skew_matrix(vu);
  double cosvu = v.dot(u);
  //cout << "Rmatrix cosvu == " << cosvu << endl;
  if (cosvu == -1) {
    cout << "Rmatrix with cosvu == -1" << endl;
    return rx180; // can't do this
  }
  double sfactor = 1/(1 + cosvu);
  return Matrix3d::I() + skm + sfactor * skm*skm;
}

