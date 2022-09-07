/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/


#ifndef __VECTOR3D_HPP__
#define __VECTOR3D_HPP__

#include <string>
#include <cassert>
#include <cmath>
#include <iostream>

using std::abs;
using std::cos;
using std::sin;
using std::sqrt;

// 
namespace Base
{
  extern double FLOAT_TOL;
  extern int NOTCART;
}

// took this Vector class from SAMoS when the code was much smaller and 
// had no dependencies

// Now it would be perhaps preferable to use Eigen vectors and matrices
// Eigen doesn't need python bindings because we can directly pass numpy arrays by copying


/*! Vector3d class
 *  Handles vectors in 3d Eucledian space
 */
class Vector3d
{
public:
  
  //! Deault constructor
  Vector3d() : x(0.0), y(0.0), z(0.0) { }
  //! Constructor for a vector object
  Vector3d(double x, double y, double z) : x(x), y(y), z(z) { }

  Vector3d operator=(const Vector3d rhs);
   
  //! Add two vectors
  Vector3d operator+(const Vector3d v) const;
  
  //! Subtract two vectors (constant version)
  Vector3d operator-(const Vector3d v) const;
  
  //! Negate a vector
  Vector3d operator-() const;
  
  //! Scale vector by a constant
  Vector3d operator*(const double c) const;
  Vector3d operator/(const double c) const;
  
  //
  Vector3d operator*=(const double c);
  
  //! Test equality
  bool operator==(const Vector3d v) const;
  
  //! Add vector to current vector
  Vector3d operator+=(const Vector3d v);
  
  //! Subtract vector from current vector
  Vector3d operator-=(const Vector3d v);
  
  //! Euclidean dot product with another vector
  double dot(const Vector3d v);
  
  //! Cross prduct with another vector
  Vector3d cross(const Vector3d v);

  double xycross(const Vector3d v);
  Vector3d xyproject();
  Vector3d xzproject();
  
  //! Vector length 
  double len(); 
  double xylen(); 
  
  //! Vector length squared
  double len2(); 
  
  //! Rescale vactor
  void scale(double s); 
  
  //! Return rescale vactor
  Vector3d scaled(double s); 
  
  //! normalise
  void norm();
  
  //! norm but protected so that zero vector returns zero vector
  Vector3d prot_unit();
  
  //! Return unit vector in the direction of this vector
  Vector3d unit();
  
  //! Get the part of the vector parallel to a given vector
  Vector3d parallel_projection(const Vector3d v);
      
  //! Get the part normal to a given vector
  Vector3d perp_projection(const Vector3d v);
  
  Vector3d rotate(const double phi, const Vector3d v);
  
  Vector3d xyrotate(const double phi);

  Vector3d perp2d();

  std::string __str__() const;

  // how slow is this and does it matter?
  double& operator[] (const int index);
  
  ///@{
  double x, y, z;
  //@}

};


// should clean up these un-inlines inline methods

//! Compute cross product between two vectors
Vector3d cross(const Vector3d v1, const Vector3d v2);

//! Scale vector by a number
Vector3d operator*(const double c, const Vector3d v);

//! Compute dot product between two vectors
double dot(const Vector3d v1, const Vector3d v2);

// Angle in right handed coordinate system
//! \param a first vector
//! \param b second vector
double angle(Vector3d a, Vector3d b);

const Vector3d e_x = Vector3d(1.,0.,0.);
const Vector3d e_y = Vector3d(0.,1.,0.);
const Vector3d e_z = Vector3d(0.,0.,1.);

double theta(Vector3d vec);

Vector3d theta_axis(double theta);

std::ostream& operator<<(std::ostream& os, Vector3d v);

Vector3d mirror2d(Vector3d pt, Vector3d axis);

bool tolequals(Vector3d v, Vector3d u, double tol);
bool tolequals(Vector3d v, Vector3d u);


///////////////////////////////////////////////////////////////////////////////////




#endif

