
#include "frame.hpp"
#include <boost/format.hpp>
#include <iostream>

using std::cout;
using std::endl;

//https://stackoverflow.com/questions/34391968/how-to-find-the-rotation-matrix-between-two-coordinate-systems
// M.T()
// Vector Frame -> Lab
Matrix3d Frame::get_rmatrix() {
  // e1 e2 e3 
  return Matrix3d(
      e1.x, e2.x, e3.x, 
      e1.y, e2.y, e3.y, 
      e1.z, e2.z, e3.z
      );
}

// M
// Vector Lab -> Frame
Matrix3d Frame::get_tmatrix() {
  // e1
  // e2
  // e3
  return Matrix3d(
      e1.x, e1.y, e1.z, 
      e2.x, e2.y, e2.z, 
      e3.x, e3.y, e3.z
      );
}

// rotate vector from this frame into lab frame
Vector3d Frame::to_lab_rt(Vector3d v) {
  return get_rmatrix() * v;
}

// find the actual point targeted by this vector in lab frame
Vector3d Frame::to_lab(Vector3d v) {
  // assume lab frame has origin at \vec{0}
  return origin + to_lab_rt(v);
}

std::array<double, 3> Frame::orthogonal_error() {
  return std::array<double, 3>{e1.dot(e2), e1.dot(e2), e2.dot(e3)};
}
// the length of the basis vectors should be unity
std::array<double, 3> Frame::unit_error() {
  return std::array<double, 3>{e1.len(), e2.len(), e3.len()};
}

void Frame::cross_norm() {
  Vector3d c23 = e2.cross(e3);
  e1 = c23.unit();
  Vector3d c31 = e3.cross(e1);
  e2 = c31.unit();
  e3 = e1.cross(e2);
}

// re-orthogonalise matrix to account for numerical drift
//https://stackoverflow.com/questions/23080791/eigen-re-orthogonalization-of-rotation-matrix
double Frame::orthogonalise() {
  // assume error is equal in all directions
  double error = e1.dot(e2); // zero for perfect R matrix
  //cout << "orthogonalise error " << error << endl;
  Vector3d ne1 = e1 - (error/2.)*e2;
  Vector3d ne2 = e2 - (error/2.)*e1;
  Vector3d ne3 = ne1.cross(ne2);
  // normalise all vectors
  ne1.unit();
  ne2.unit();
  ne3.unit();
  this->e1 = ne1;
  this->e2 = ne2;
  this->e3 = ne3;
  return error;
}

// some facet of C++ which I don't understand prevents this method from 
// influencing the object state correctly
// It must act on a copy and not the object itself?
// Or not set the resulting values
double Frame::normalise() {
  //std::cout << "CALLING NORMALISE" << std::endl;
  e1 = this->e1.unit();
  e2 = this->e2.unit();
  e3 = this->e3.unit();
}


std::string Frame::__str__() {
  return ( boost::format("Frame: \n->origin%s \n->e1%s \n->e2%s \n->e3%s") 
      % origin.__str__() % e1.__str__() % e2.__str__() % e3.__str__() 
      ).str();
}

std::ostream& operator<<(std::ostream& os, Frame& frame)
{
  os << frame.__str__();
  return os;
}


int which_cart(Vector3d v) {
  int i;
  if (tolequals(v, e_x)) {
    i = 0;
  }
  else if (tolequals(v, e_y)) {
    i = 1;
  }
  else if (tolequals(v, e_z)) {
    i = 2;
  }
  else { i = Base::NOTCART; }
  return i;
}

std::tuple<int,int,int> cartesian_basis_vectors(Frame& frame)
{
  int i1,i2,i3;
  i1 = which_cart(frame.e1);
  i2 = which_cart(frame.e2);
  i3 = which_cart(frame.e3);
  return std::make_tuple(i1,i2,i3);
}


