
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <typeinfo>

#include <memory>
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "axisangle.hpp"

#include "lvec.hpp"
#include "chain.hpp"
#include "sphere.hpp"
#include "shape.hpp"
#include "plane.hpp"

#include "pili.hpp"
#include "cell.hpp"
#include "cell3d.hpp"

#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>

using doctest::Approx;
using std::get;
using std::cout;
using std::endl;
using std::vector;

// for checking floating point equality
double float_tol = 0.00001;

bool gequals(Vector3d v, Vector3d u) { return tolequals(v,u,float_tol); }
  
TEST_CASE( "pt_in_sphere" )
{
  Sphere S = Sphere(1.);
  Vector3d pt = Vector3d(0.5,0.5,0.);
  CHECK( S.contains(pt) == true);

  Vector3d pt2 = Vector3d(1.5,0.5,0.);
  CHECK( S.contains(pt2) == false);

}
TEST_CASE( "closest_to_lvec" )
{
  // from e_z to e_z + e_x
  Lvec lv = Lvec(Vector3d(0.,0.,1.), Vector3d(1.,0.,0.));
  Vector3d pt1 = Vector3d(); // (0.,0.,0.)
  cout << lv << endl;
  cout << lv.closest_pt(pt1)  << endl;
  CHECK( lv.closest_pt(pt1) == Vector3d(0.,0.,1.) );

  Vector3d pt2 = Vector3d(0.5,0.,-1.);
  CHECK( lv.closest_pt(pt2) == Vector3d(0.5,0.,1.) );

  Vector3d pt3 = Vector3d(1.5,0.,-1.); 
  CHECK( lv.closest_pt(pt3) == Vector3d(1.,0.,1.) );
}

TEST_CASE( "distance_to_lvec" )
{
  Lvec lv = Lvec(Vector3d(), Vector3d(1.,0.,0.));
  Vector3d pt1 = Vector3d(); // (0.,0.,0.)
  CHECK( lv.distance_vector(pt1) == Vector3d(0.,0.,0.) );

  Vector3d pt2 = Vector3d(-1.,0.,0.); // (-1., 0., 0.)
  CHECK( lv.distance_vector(pt2) == Vector3d(-1.,0.,0.) );

  Vector3d pt3 = Vector3d(2.,1.,0.); // (1., 1., 0.)
  CHECK( lv.distance_vector(pt3) == Vector3d(1.,1.,0.) );

  Vector3d pt4 = Vector3d(0.2,5.,0.); // (0., 5., 0.)
  CHECK( lv.distance_vector(pt4) == Vector3d(0.,5.,0.) );

;
  Lvec lv2{Vector3d(0,0,1),Vector3d(0,0,-2)};
  Vector3d pt5 = Vector3d(-0.2534,-0.3839,-0.6); 
  CHECK( lv2.distance_vector(pt5).len() == Approx(0.46) );

}


TEST_CASE( "pt_to_plane" )
{
  Plane plane = Plane();
  Vector3d pt1 = Vector3d(); // (0.,0.,0.)
  CHECK( plane.project(pt1) == Vector3d(0,0,0) );
  CHECK( plane.distance_vector(pt1) == Vector3d(0,0,0) );

  Vector3d pt2 = Vector3d(3.,4.,1.); 
  CHECK( plane.project(pt2) == Vector3d(3,4,0) );
  CHECK( plane.distance_vector(pt2) == Vector3d(0,0,1) );

  Vector3d pt3 = Vector3d(3.,4.,-1.); 
  CHECK( plane.distance_vector(pt3) == Vector3d(0,0,-1) );

  
}

/*
TEST_CASE( "lvec_plane_intersect" )
{

  std::shared_ptr<Plane> plane = std::make_shared<Plane>();
  Vector3d pt1 = Vector3d(); // (0.,0.,0.)
  Lvec lv1 = Lvec(Vector3d(0., 0., 1.), Vector3d(0., 1., 2.));
  //CHECK( lv1.intersect(plane) == NULL );
  CHECK( !lv1.intersect(*plane) );

  Lvec lv2 = Lvec(-Vector3d(0., 0., 1.), -Vector3d(0., 1., 2.));
  CHECK( !lv2.intersect(*plane) );

  Lvec lv3 = Lvec(Vector3d(1., 0., -1.), Vector3d(2., 0., 2.));
  //cout << lv3.intersect(plane) << endl;
  CHECK( *lv3.intersect(*plane) == Vector3d(2., 0., 0.) );


  Lvec lv4{e_z, Vector3d(0.972979, 0.000000, -0.034979)-e_z};
  //cout << "intersects? " << *lv4.intersect(plane) << endl;
  CHECK( lv4.intersect(*plane) );

}
*/

TEST_CASE( "frame_norm" )
{
  Frame ext{Vector3d(), 2*e_x, 2*e_y};
  ext.normalise();
  CHECK( tolequals(ext.e1, e_x, float_tol) );
  CHECK( tolequals(ext.e2, e_y, float_tol) );
  CHECK( tolequals(ext.e3, e_z, float_tol) );
}

TEST_CASE( "frame_rmatrix" )
{

  Frame zero{};
  Frame frame{Vector3d(), e_x};
  zero.rotate( frame.get_rmatrix() );
  CHECK( tolequals( zero.e1, frame.e1, float_tol) );
  CHECK( tolequals( zero.e2, frame.e2, float_tol) );
  CHECK( tolequals( zero.e3, frame.e3, float_tol) );
}

const double R = 0.5;
const double l = 2;

TEST_CASE( "pili_spherical" )
{
  //Frame fez{Vector3d(0,0,R), e_z};
  Frame fex{Vector3d(0,0,0), e_x};
  Capsule caps{Vector3d(0, 0, 0.0), Vector3d(1,0,0).unit(), R, l};
  
  CHECK( tolequals(caps.alt_spherical(0, 0), R*e_z, float_tol) );
  CHECK( tolequals(caps.alt_spherical(0, M_PI/2.), R*e_z, float_tol) );
}


TEST_CASE( "test_capsule_params" )
{
  Capsule caps{Vector3d(0, 0, 0.4), Vector3d(1,0,-1).unit(), 0.5, 1};
  //cout << caps.halfal() << endl;
  Capsule caps_2{Vector3d(0, 0, 0.4), Vector3d(1,0,-1).unit(), 0.5, 2};
  //cout << caps_2.halfal() << endl;
}

TEST_CASE( "test_capsule_overlap" )
{
  Plane plane = Plane();
  Capsule caps{Vector3d(0, 0, 0.4), Vector3d(1,0,-1).unit(), 0.5, 2};
  // needs updated
  //boost::optional<overlap> vv = caps.overlap_vector(plane);
  //CHECK( vv );
  //CHECK( abs(get<0>(*vv) - 0.1) < float_tol );
  //CHECK( get<1>(*vv) == Vector3d(0,0,1) );
  //CHECK( tolequals( get<2>(*vv), Vector3d(0,0,-0.1), float_tol) );
}

/*TEST_CASE( "test_capsule_surfacept" )*/
//{
  //Capsule caps{Vector3d(0, 0, 0), Vector3d(1,0,0).unit(), R, l};
  //double sq = (pi*R)/2.;

  //CHECK( gequals( caps.spt(0., 0.), Vector3d(0,0,1.5) ) );
  //CHECK( gequals( caps.spt(sq, 0.), Vector3d(0.5,0,1) ) );
  //CHECK( gequals( caps.spt(sq + 1., 0.), Vector3d(0.5,0,0) ) );
  //CHECK( gequals( caps.spt(sq + 2., 0.), Vector3d(0.5,0,-1) ) );
  //CHECK( gequals( caps.spt(2*sq + 2. -float_tol/10., 0.), Vector3d(0, 0, -1.5) ) );
  //CHECK( gequals( caps.spt(2*sq + 2., 0.), Vector3d(0, 0, -1.5) ) );
  ////CHECK( gequals( caps.spt(sq, pi/2.), Vector3d(0,0.5,0) ) );
//}

/*TEST_CASE( "test_alt_spherical" )*/
//{
  //Capsule caps{Vector3d(0, 0, 0), Vector3d(1,0,0).unit(), R, l};
  //double sq = (pi*R)/2.;

  //Vector3d b1 = caps.spt(0., 0.);
  //Vector3d b2 = caps.spt(sq, 0.);
  //Vector3d b3 = caps.spt(sq + 1., 0.);
  //Vector3d b4 = caps.spt(sq + 2., 0.);
  //Vector3d b5 = caps.spt(2*sq + 2., 0.);
  //Vector3d b6 = caps.spt(sq, pi/2.);

  //Matrix3d M_tolab = caps.frame.get_rmatrix();
  //CHECK( gequals( M_tolab * b1, Vector3d(0.5,0,0) ) );
  //CHECK( gequals( M_tolab * b2, Vector3d(0,0.5,0) ) );
  //CHECK( gequals( M_tolab * b3, Vector3d(-1,0.5,0) ) );
  //CHECK( gequals( M_tolab * b4, Vector3d(-2,0.5,0) ) );
  //CHECK( gequals( M_tolab * b5, Vector3d(-2.5,0,0) ) );
  //CHECK( gequals( M_tolab * b6, Vector3d(0,0,-0.5) ) );
/*}*/


/*TEST_CASE( "Lvec_Step_interaction" )*/
//{
  ////cout << "test Step intersection " << endl;
  //double h = 1;
  //double B = 4;

  //Plane *p = new Step(h, B);

  //Lvec lv = Lvec(Vector3d(-0.1, 0, 0.5), Vector3d(1., 0, 0));
  //boost::optional<Vector3d> inter = p->intersects(lv);
  //if (inter) {
  //}
  //else {
  //}
/*}*/

TEST_CASE( "chain_rotation" )
{
  Vector3d source = Vector3d(0,0,0);
  vector<Vector3d> vc{Vector3d(1.,1.,0.),Vector3d(0.,0.,-2.),Vector3d(-1,-1,0)};
  Chain ch{vc};
  Matrix3d toez = Rmatrix(Vector3d(1,1,0).unit(),e_z);
  Matrix3d fromez = Rmatrix(e_z,Vector3d(1,1,0).unit());
  ch.rotate(toez);
  ch.rotate(fromez);
  ch.compute_targets();
}

TEST_CASE( "plane_chain_intersection" )
{
  Plane plane = Plane();
  Vector3d source = Vector3d(0,0,1);
  vector<Vector3d> vc{Vector3d(1.,1.,0.),Vector3d(0.,0.,-2.),Vector3d(-1,-1,0)};
  Chain ch{vc};
  ch.translate(source);
  ch.compute_targets();

  boost::optional<Vector3d> inter = plane.intersects(ch);
  //cout << *inter << endl;
  CHECK( inter );
  CHECK( tolequals(*inter, Vector3d(1,1,0), float_tol) );
}


TEST_CASE( "capsule_sphere_intersection" )
{
  //cout << "Test capsule sphere intersection" << endl;
  Capsule caps{Vector3d(0, 0, 1.0), Vector3d(1,0,0).unit(), 0.5, 2};
  Sphere sph1{Vector3d(0,0,0), 0.55};
  boost::optional<overlap> over1 = caps.intersects(sph1);
  Sphere sph2{Vector3d(0,0,0), 0.45};
  boost::optional<overlap> over2 = caps.intersects(sph2);

  //cout << "has intersection" << endl;
  CHECK( bool(over1) );
  //cout << "No intersection" << endl;
  CHECK( !bool(over2) );
}


TEST_CASE( "sphere_line_intersection" )
{
  //cout << "Test sphere line intersection" << endl;
  Sphere sph1{Vector3d(0,0,0), 1};
  Lvec lv1 = Lvec(Vector3d(0,-2,0.5), Vector3d(0,5,0));
  boost::optional<Vector3d> inter = sph1.intersects(lv1);
  CHECK( tolequals(*inter, Vector3d(0, -0.866025, 0.500), float_tol) );
  Lvec lv2 = Lvec(Vector3d(0,-2,1.5), Vector3d(0,5,0));
  boost::optional<Vector3d> inter2 = sph1.intersects(lv2);
  CHECK( !bool(inter2) );
}


// from test_base

// TEST_CASE( "test_aarmatrix" )
// {
//   Vector3d p1{0, M_PI, 0};
//   Vector3d p2{0, -M_PI, 0};
//   //cout << aarmatrix(p1) << endl;
//   //cout << aarmatrix(p2) << endl;
//   Rk brk1{p1};
//   Rk brk2{p2};
//   cout << "from p1 " << endl;
//   cout << brk1.Rk1 << endl;
//   cout << brk1.Rk2 << endl;
//   cout << brk1.Rk3 << endl;
//   cout << "--" << endl;
//   cout << "from p2 " << endl;
//   cout << brk2.Rk1 << endl;
//   cout << brk2.Rk2 << endl;
//   cout << brk2.Rk3 << endl;
//   cout << "--" << endl;
// }
