
#include "vector3d.hpp"
#include <boost/format.hpp>

namespace Base
{
  double FLOAT_TOL = 1e-8;
  int NOTCART = 9;
}



// we implement copy assigment so we should probably implement all special methods
Vector3d Vector3d::operator=(const Vector3d rhs)
{
  x = rhs.x;
  y = rhs.y;
  z = rhs.z;
  return *this;
}

Vector3d Vector3d::operator+(const Vector3d v) const
{
  return Vector3d(x + v.x, y + v.y, z + v.z);
}

Vector3d Vector3d::operator-(const Vector3d v) const
{
  return Vector3d(x - v.x, y - v.y, z - v.z);
}

Vector3d Vector3d::operator-() const
{
  return Vector3d(-x, -y, -z);
}

Vector3d Vector3d::operator*(const double c) const
{
  return Vector3d(c*x, c*y, c*z);
}

Vector3d Vector3d::operator/(const double c) const
{
  return Vector3d(x/c, y/c, z/c);
}

Vector3d Vector3d::operator*=(const double c)
{
  x *= c; y *= c; z *= c;
  return *this;
}

bool Vector3d::operator==(const Vector3d v) const
{
  return (x == v.x && y == v.y && z == v.z);
}

Vector3d Vector3d::operator+=(const Vector3d v)
{
  x += v.x;  y += v.y;  z += v.z;
  return *this;
}

Vector3d Vector3d::operator-=(const Vector3d v)
{
  x -= v.x;  y -= v.y;  z -= v.z;
  return *this;
}

double Vector3d::dot(const Vector3d v)
{
  return x*v.x + y*v.y + z*v.z;
}

Vector3d Vector3d::cross(const Vector3d v)
{
  return Vector3d(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
}

double Vector3d::xycross(const Vector3d v)
{
  return x*v.y - y*v.x;
}

Vector3d Vector3d::xyproject()
{
  return Vector3d(this->x, this->y, 0);
}
Vector3d Vector3d::xzproject()
{
  return Vector3d(this->x, 0, this->z);
}

double Vector3d::len() { return sqrt(x*x + y*y + z*z); }
double Vector3d::xylen() { return sqrt(x*x + y*y); }
double Vector3d::len2() { return x*x + y*y + z*z; }
void Vector3d::scale(double s) { x *= s; y *= s;  z *= s;  }
Vector3d Vector3d::scaled(double s) { return Vector3d(s*x,s*y,s*z);  }


void Vector3d::norm()
{
  double len = this->len();
  if (len != double(0))
    this->scale(double(1)/len);
}


Vector3d Vector3d::prot_unit()
{
  double len = this->len();
  if (len == 0.)
    return Vector3d();
  return Vector3d(x/len,y/len,z/len);
}


Vector3d Vector3d::unit()
{
  double len = this->len();
  if (len == double(0)) {
    return Vector3d();
  }
  return Vector3d(x/len,y/len,z/len);
}

Vector3d Vector3d::parallel_projection(const Vector3d v)
{
  Vector3d n = v;
  n.norm();
  double proj = this->dot(n);
  return Vector3d(proj*n.x, proj*n.y, proj*n.z);
}

Vector3d Vector3d::perp_projection(const Vector3d v)
{
  Vector3d n = this->parallel_projection(v);
  return Vector3d(x-n.x, y-n.y, z-n.z);
}

Vector3d Vector3d::rotate(const double phi, const Vector3d v)
{
  double s = std::sin(phi);
  double c = std::cos(phi);
  double k = 1.0 - c;
  
  double nx = x * (c + k * v.x * v.x) + y * (k * v.x * v.y - s * v.z) + z * (k * v.x * v.z + s * v.y);
  double ny = x * (k * v.x * v.y + s * v.z) + y * (c + k * v.y * v.y) + z * (k * v.y * v.z - s * v.x);
  double nz = x * (k * v.x * v.z - s * v.y) + y * (k * v.y * v.z + s * v.x) + z * (c + k * v.z * v.z);
  
  return Vector3d(nx,ny,nz);
}

Vector3d Vector3d::xyrotate(const double phi)
{
  double s = std::sin(phi);
  double c = std::cos(phi);
  double nx = c * x  -s * y;
  double ny = s * x + c * y;
  return Vector3d(nx,ny,z);
}

Vector3d Vector3d::perp2d() 
{
  return Vector3d(-y, x, z);
}

//////////////////////////////////////////////////////////////////////////////////


double theta(Vector3d vec)
{
  return angle(e_x, vec);
}

Vector3d theta_axis(double theta)
{
  return Vector3d(cos(theta), sin(theta), 0.);
}


// signed angle makes sense in xy plane
double angle(Vector3d a, Vector3d b)
{
  double a_dot_b = dot(a.unit(),b.unit());
  // these checks ensure no errors from numerical rounding
  if (a_dot_b > 1.0) a_dot_b = 1.0;
  else if (a_dot_b < -1.0) a_dot_b = -1.0;
  double phi = std::acos(a_dot_b);
  double sign = cross(a,b).z;
  if (sign >= 0) 
    return phi;
  else
    return -phi;
}



std::string Vector3d::__str__() const {
  return ( boost::format("(%12.6f,%12.6f,%12.6f)") % x % y % z ).str();
}

// how slow?
double& Vector3d::operator[] (const int index) 
{
  if (index == 0) return x;
  if (index == 1) return y;
  return z;
}

std::ostream& operator<<(std::ostream& os, Vector3d v)
{
  os << v.__str__();
  return os;
}


Vector3d cross(const Vector3d v1, const Vector3d v2)
{
  return Vector3d(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

Vector3d operator*(const double c, const Vector3d v)
{
  return Vector3d(c*v.x, c*v.y, c*v.z);
}

double dot(const Vector3d v1, const Vector3d v2)
{
  return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

Vector3d mirror2d(Vector3d pt, Vector3d ax)
{
  
  ax = ax.unit();
  double pdotx = pt.dot(ax);
  double x = 2 * pdotx * ax.x - pt.x;
  double y = 2 * pdotx * ax.y - pt.y;
  return Vector3d(x, y, pt.z);
}



bool tolequals(Vector3d v, Vector3d u, double tol) {
  return (
      abs(v.x-u.x) < tol &&
      abs(v.y-u.y) < tol &&
      abs(v.z-u.z) < tol
      );
}

bool tolequals(Vector3d v, Vector3d u) { return tolequals(v, u, Base::FLOAT_TOL); }
