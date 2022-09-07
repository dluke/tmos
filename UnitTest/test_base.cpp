

#include <time.h>
#include <iostream>
#include <Eigen/Dense>

#include "sample.hpp"
#include "point.hpp"
#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "axisangle.hpp"

#define BOOST_TEST_MODULE BaseTest
#include <boost/test/unit_test.hpp>

using std::cout;
using std::endl;
using std::vector;

// test the point class
// check I implemented new operators * and / correctly
BOOST_AUTO_TEST_CASE( point )
{
  Point2d p1{1,2};
  Point2d phalf = p1/2;
  //cout << phalf << endl;
  //cout << (p1 == Point2d(0,0)) << endl;
  //cout << (phalf == Point2d(0,1)) << endl;
}


BOOST_AUTO_TEST_CASE( sample )
{
  // Generator is not randomly seeded by default
  long t = long( time(NULL) ); 
  // awful method of seeding, only for testing code
  //cout << "using time seed " << t << endl;
  init_generator(t);

  int n = 5;
  Eigen::RowVectorXd phir = phi_row(n);
  //cout << phir << endl;
  BOOST_REQUIRE( phir.size() == n );

}


BOOST_AUTO_TEST_CASE( matrix3d_Rmatrix )
{
  double tol = 0.000000001;

  Vector3d ez{0.,0.,1.};
  Vector3d ey{0.,1.,0.};
  Matrix3d My = Rmatrix(ez, ey);
  Matrix3d Mz = Rmatrix(ey, ez);
  Vector3d recover_ey = My * ez;
  Vector3d recover_ez = Mz * ey;
  BOOST_REQUIRE( tolequals(recover_ey, ey, tol) );
  BOOST_REQUIRE( tolequals(recover_ez, ez, tol) );
}

