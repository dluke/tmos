
#include <cmath>
#include <iostream>
#include <vector>

#include "frame.hpp"
#include "vector3d.hpp"
#include "point.hpp"

#include "lvec.hpp"
#include "chain.hpp"
#include "sphere.hpp"
#include "shape.hpp"
#include "plane.hpp"
#include "hexsurface.hpp"
#include "curvestep.hpp"
#include "periodic.hpp"

#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>

#define BOOST_TEST_MODULE SurfaceTest
#include <boost/test/unit_test.hpp>

// boost floating point equality convenience
#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

const double ftol = 0.000001;

using std::cout;
using std::endl;
using std::vector;

BOOST_AUTO_TEST_CASE( shape ) {
  Shape* shape;
}


BOOST_AUTO_TEST_CASE( hex_coordinates ) {
  Hexc hc{1,2};
  Cubec cube = hc.cube();
  BOOST_TEST( cube[2] == -3 );
  Hexc hx = cube.hex();
  BOOST_TEST( hx == hc );
}

BOOST_AUTO_TEST_CASE( absolute_map ) {
  //cout << "Test Hexgrid::absmap" << endl;
  Hexc hr{0,1};
  //Hexc hq{1,0};
  //Hexc hz{-1,-1};
  BOOST_TEST( -hr == Hexc(0,-1) );
  BOOST_TEST( HexGrid::absmap[-hr] == hr );
}

BOOST_AUTO_TEST_CASE( test_basis, * utf::tolerance(0.000001) ) {
  //cout << "Test Hgrid basis" << endl;
  double R = 1.;
  HexGrid Hgrid{R};
  BOOST_TEST( Hgrid.hrr.len() == 2*R );
  BOOST_TEST( Hgrid.hrq.len() == 2*R );
}

BOOST_AUTO_TEST_CASE( get_rq ) {
  //cout << "Test get Hexc from vector pt" << endl;
  double R = 1.;
  HexGrid Hgrid{R};

  Hexc hr = Hgrid.get_rq(Hgrid.hrr);
  Hexc hq = Hgrid.get_rq(Hgrid.hrq);
  Hexc hz = Hgrid.get_rq(Hgrid.hrz);
  BOOST_TEST( hr == HexGrid::hr );
  BOOST_TEST( hq == HexGrid::hq );
  BOOST_TEST( hz == HexGrid::hz );
}

BOOST_AUTO_TEST_CASE( get_xyz ) {
  //cout << "Test get Hexc from vector pt" << endl;
  double R = 1.;
  HexGrid Hgrid{R};

  Vector3d hrr = Hgrid.get_xyz(HexGrid::hr);
  Vector3d hrq = Hgrid.get_xyz(HexGrid::hq);
  Vector3d hrz = Hgrid.get_xyz(HexGrid::hz);

  BOOST_TEST( tolequals(hrr, Hgrid.hrr, ftol) );
  BOOST_TEST( tolequals(hrq, Hgrid.hrq, ftol) );
  BOOST_TEST( tolequals(hrz, Hgrid.hrz, ftol) );
}

BOOST_AUTO_TEST_CASE( hex_line ) {
  
  HexGrid Hgrid{1.};
  //Vector3d target{5.,5.,0.};
  Vector3d target{0.1,0.1,0.};
  Lvec lv{Vector3d(0,0,0), target};

  vector<Hexc> hline = Hgrid.line(lv);
  cout << "hexes in line " << lv << endl;
  for (Hexc hx : hline) {
    cout << hx << endl;
  }

}

BOOST_AUTO_TEST_CASE( chain_set ) {
  cout << "Test HexSphereGrid::chain_set" << endl;

  HexSphereGrid Hgrid{1.};
  Vector3d source{0,0,1.2};
  Vector3d axis{1,0,0};
  double a = 2;
  vector<Vector3d> vc{
    a*Vector3d(1,0,0).unit(),
    a*Vector3d(0,-1,-1).unit()
  };
  Chain ch{source, axis, a, vc};
  ch.compute_targets();
  cout << "Test Chain" << endl;
  cout << ch.__str__() << endl;

  cout << "Chain intersection" << endl;
  boost::optional<Vector3d> inter = Hgrid.intersects(ch);
  if (inter)
    cout << *inter << endl;
  else
    cout << "none" << endl;
  //cout << *inter << endl;

}

BOOST_AUTO_TEST_CASE( caps_overlaps ) {
  cout << "Test HexSphereGrid::overlap_vector(Capulse caps) " << endl;
  HexSphereGrid Hgrid{1.};
  // head pointing +y +deltaz, tail intersect -y direction
  Capsule caps{Vector3d(0, 0, 0.49), Vector3d(1,1,0.1).unit(), 0.5, 4};
  cout << "Here caps" << endl;
  cout << caps.__str__() << endl;
  cout << caps.get_endpt().__str__() << endl;

  vector<overlap> overs = Hgrid.overlap_vector(caps);
  cout << "overlaps ..." << endl;
  for (overlap over : overs) {
    cout << over << endl;
  }

}

BOOST_AUTO_TEST_CASE( axisangle_matrix )
{
  Vector3d p{-2.92616927, -1.13154542, -0.15923197};
  //Vector3d p{-3.92616927, -1.13154542, -0.15923197};
  Matrix3d M = aarmatrix(p);
  //cout << M * e_z << endl;
  //cout << M * e_x << endl;
  //cout << e_
  cout << aarmatrix(p) << endl;

  Capsule caps{Vector3d(), e_z, 0.5, 2.};
  caps.set_rotation(p);
  Vector3d pdash = caps.get_axisangle();
  cout << pdash << endl;
}



BOOST_AUTO_TEST_CASE( curve_element )
{
  cout << "-----------------------------------------" << endl;
  cout << "CurveElement Line Segment Intersection." << endl;

  double panewidth = 1.;
  double planeangle = M_PI/6.;
  //PartialPlane plane{Vector3d(), e_z, e_x, 0, 1};
  CurveElement celement{Frame(), planeangle, panewidth};
  double xperiod = celement.get_xlength();
  double delta = 0.01;
  // forwards
  Lvec lv{e_z, e_x-e_z};
  // backwards
  Lvec lvr{e_z+celement.get_xlength()*e_x, -1*e_x-e_z};
  // backwards, two intersections
  Lvec lvhr{0.5*e_z+celement.get_xlength()*e_x, -celement.get_xlength()*e_x};
  // backwards, end segment
  Lvec lvbe{celement.get_lmpt() + e_x, -2* e_x - 0.1*e_z};
  // backwards, first segment
  Lvec lvbf{celement.get_lstart() + 0.5*e_x, -2*e_x - 0.1*e_z};
  // forwards, off the end
  Lvec lvfe{celement.get_lmpt() - e_x, 2*e_x - 0.1*e_z};
  // forwards, onto the start
  Lvec lvff{celement.get_lstart() - e_x, 2*e_x - 0.1*e_z};
  Lvec offleft{celement.get_lstart() - e_x, -1*e_x};
  Lvec offright{celement.get_lmpt() + e_x, e_x};

  Lvec mrvertical{Vector3d(xperiod/2. +delta, 0, 1.0),-2*e_z};
  Lvec mlvertical{Vector3d(xperiod/2. -delta, 0, 1.0),-2*e_z - delta*e_x }; // reverse direction
  
  //for (auto pplane : celement.panes) {
    //cout << pplane.get_origin() << endl;
  //}
  //cout << celement.get_lmpt() << endl;

  cout << *celement.intersects(lvr) << endl;
  cout << *celement.intersects(lvhr) << endl;
  cout << *celement.intersects(lvbe) << endl;
  cout << *celement.intersects(lvbf) << endl;
  cout << *celement.intersects(lvfe) << endl;
  cout << *celement.intersects(lvff) << endl;
  cout << mrvertical << endl;
  cout << *celement.intersects(mrvertical) << endl;
  cout << mlvertical << endl;
  cout << *celement.intersects(mlvertical) << endl;

}


BOOST_AUTO_TEST_CASE( segplane_construct )
{
  cout << "-----------------------------------------" << endl;
  cout << "SegPlane construction" << endl;
  double panewidth = 1.;
  double planeangle = M_PI/6.;
  
  SegPlane sp{Vector3d(), e_z, planeangle, panewidth};
  // convenience reference to repeating distance and surface height
  cout << "zmax " << sp.zmax << endl;
  cout << "xperiod " << sp.xperiod << endl;
  cout << "size " << sp._elements.size() << endl;
  cout << "element 0" << endl;
  cout << *sp.get_element_at(0) << endl;
  cout << "element -1" << endl;
  cout << *sp.get_element_at(-1) << endl;
}

BOOST_AUTO_TEST_CASE( segplane_overlap )
{
  cout << "-----------------------------------------" << endl;
  cout << "Body overlap" << endl;
  double R = 0.5;
  double length = 2;
  double delta = 0.01;
  double panewidth = 1.;
  double planeangle = M_PI/6.;
  
  SegPlane sp{Vector3d(), e_z, planeangle, panewidth};
  // convenience reference to repeating distance and surface height
  double zmax = sp.zmax;
  double xperiod = sp.xperiod;
  //CurveElement celement = sp._elements[0].get();

  Capsule caps1{Vector3d(0, 0, zmax+R-delta), e_x, R, length};
  Capsule caps2{Vector3d(xperiod/2., 0, 1.+R-delta), -1*e_z, R, length};
 //out << caps2 << endl;
  //cout << caps2.get_lvec().__str__() << endl;
  Capsule caps3{Vector3d(xperiod/2., 0, R-delta), e_y, R, length};
  cout << "1 touch at repeating plane" << endl;
  cout << sp.overlap_vector(caps1).size() << endl;
  cout << "2 touches at minimum" << endl;
  cout << sp.overlap_vector(caps2).size() << endl;
  cout << "4 touches at minimum" << endl;
  cout << sp.overlap_vector(caps3).size() << endl;

  //cout << "2 touches on surface" << endl;

}


BOOST_AUTO_TEST_CASE( sineplane_pt )
{
  cout << "-----------------------------------------" << endl;
  Frame frame;
  double A = 1.;
  double B = 1.;
  SinePlane sineplane{frame, A, B, 1};

  Vector3d pt1{-M_PI/2., 0, -A+0.5};
  Vector3d pt2{-M_PI/2., 99,-A+0.5};
  double x = M_PI/4.;
  double delta = 0.1;
  Vector3d pt3{x, 0, sineplane.form(x)+delta};
  Vector3d v1 = sineplane.distance_vector(pt1);
  Vector3d v2 = sineplane.distance_vector(pt2);
  Vector3d v3 = sineplane.distance_vector(pt3);
  BOOST_TEST( sineplane.iperiod == 0 );
  BOOST_TEST( sineplane.iconst == 1 );
  BOOST_TEST( sineplane.inorm == 2 );


  double t1 = sineplane.closest_t(pt1);
  double t2 = sineplane.closest_t(pt2);
  double t3 = sineplane.closest_t(pt3);

  cout << sineplane.closest_t(pt1) << endl;
  cout << sineplane.closest_t(pt2) << endl;
  cout << sineplane.closest_t(pt3) << endl;
  cout << "--" << endl;

  cout << v1 << endl;
  cout << v2 << endl;
  cout << v3 << endl;

  cout << "(1) == (2)" << endl;
  cout << (Vector3d(t1, 0, sineplane.form(t1)) - pt1).len() << endl;
  cout << (Vector3d(pt1.x, 0, sineplane.form(pt1.x)) - pt1).len() << endl;


}

BOOST_AUTO_TEST_CASE( sineplane_lv )
{
  cout << "-----------------------------------------" << endl;
  Frame frame;
  double A = 1.;
  double B = 1.;
  SinePlane sineplane{frame, A, B, 1};

  double tiny = 1e-4;
  double small = 0.1;
  double delta = 0.5;

  // (0,0,0)
  Lvec lv1{Vector3d(-delta,0,0),Vector3d(2*delta,0,0)};
  // none
  Lvec lv2{Vector3d(-delta,0,A+small),Vector3d(2*delta,0,0)};
  // none
  Lvec lv3{Vector3d(-delta-M_PI/2.,0,A-small),Vector3d(2*delta,0,0)};
  // (0,0,0)
  Lvec lv4{Vector3d(0,0,delta),Vector3d(0,0,-2*delta)};
  // none
  Lvec lv5{Vector3d(M_PI/2.,0,delta),Vector3d(-tiny,0,-2*delta)};
  // (pi/2., 0, 1)
  Lvec lv6{Vector3d(M_PI/2.,0,A+delta),Vector3d(0,0,-2*delta)};
  // left down of (pi/2, 0, 1)
  Lvec lv7{Vector3d(M_PI/2.-delta,0,A-small),Vector3d(2*delta,0,0)};
  // right down of (pi/2, 0, 1)
  Lvec lv8{Vector3d(M_PI/2.+delta,0,A-small),Vector3d(-2*delta,0,0)};
  // right up of (-pi/2., 0, -1) (rvalue)
  Lvec lv9{
    Vector3d(-M_PI/2. - M_PI/4., 0, sineplane.form(-M_PI/2. - M_PI/4.) + small),
    Vector3d(M_PI,0,-1.)};
  // ...
  Lvec lv10 = lvec_from_pts(
      Vector3d( -0.488200,   -0.439635,    0.100335),  
      Vector3d(   -0.503501,   -0.489333,    0.100335));

  vector<Lvec> lvl;
  lvl.push_back(lv1);
  lvl.push_back(lv2);
  lvl.push_back(lv3);
  lvl.push_back(lv4);
  //lvl.push_back(lv5);
  lvl.push_back(lv6);
  lvl.push_back(lv7);
  lvl.push_back(lv8);
  lvl.push_back(lv9);
  lvl.push_back(lv10);

  int ct = 1;
  for (Lvec lv : lvl) {
    cout << "CHECKING INTERSECTION " << ct << endl;
    cout << lv << endl;
    auto inter = sineplane.intersects(lv);
    if (inter)
      cout << "Intersection " << *inter << endl;
    else
      cout << "None" << endl;
    cout << "--" << endl;
    ct++;
  }

}

BOOST_AUTO_TEST_CASE( sineplane_caps )
{
  cout << "-----------------------------------------" << endl;
  Frame frame;
  double A = 1.;
  double B = 1.;
  SinePlane sineplane{frame, A, B, 1};
  double R = 0.5;
  double length = 2.;

  double delta = 0.01;
  double small = 0.1;
  // head contact
  Vector3d n1 = sineplane.normal_form(0.);
  cout << "n1" << n1 << endl;
  Capsule cap1 = Caps::headat(Vector3d() + (R-delta)*n1, -e_z);
  // tail contact
  Capsule cap2 = Caps::tailat(Vector3d(M_PI/2, 0, A+R-delta), e_z);
  // body contact
  Capsule cap3 = Capsule(Vector3d(M_PI/2, 0, A+R-delta), e_x);
  // two contacts (constant direction)
  Capsule cap4 = Capsule(Vector3d(-M_PI/2, 0, -A+R-delta), e_y);
  // two contacts (opposite ends)
  Capsule cap5 = Capsule(Vector3d(-M_PI/2, 0, -2*small), e_x);
  // one contact
  Lvec lv6 = lvec_from_pts(
    Vector3d(   -0.428971,    1.256417,    0.263319),
    Vector3d(    0.892777,    0.545908,    1.585498));
  Capsule cap6 = Caps::body_from_lv(lv6);
  Lvec lv7 = lvec_from_pts(
    Vector3d( 1.685611,     0.639798,    -0.498219),
    Vector3d(-3.031867,     1.550677,     0.667068));
  Capsule cap7 = Caps::body_from_lv(lv7);
  //Vector3d(-7.250284403562545776e-02, 
      //-6.523967385292053223e-01,
      //6.243178844451904297e-01)
      //-6.609007696187683134e-01 3.562129624703446185e-01 -6.605471202612425152e-01

  Capsule cap8 = Capsule(
    Vector3d(-1.905893802642822266e+00, 2.231079578399658203e+00, -4.142510890960693359e-01),
    Vector3d( -3.812505605884797255e-02, -9.990990809297168873e-01, 1.864152852919498393e-02)
    );

  Capsule cap9 = Capsule(
      Vector3d(-1.873687386512756348e+00, 2.224899530410766602e+00, -4.288987815380096436e-01),
     Vector3d( -1.874481811050249247e-02, -9.997851966930920709e-01, 8.842639162569887468e-03)
     );

  vector<Capsule> lcap;
  lcap.push_back(cap1);
  lcap.push_back(cap2);
  lcap.push_back(cap3);
  lcap.push_back(cap4);
  lcap.push_back(cap5);
  lcap.push_back(cap6);
  //lcap.push_back(cap7);
  lcap.push_back(cap8);
  lcap.push_back(cap9);

  int ct = 1;
  for (Capsule cap: lcap) {
    cout << "CHECKING OVERLAPS " << ct << endl;
    cout << cap << endl;
    auto overs = sineplane.overlap_vector(cap);
    if (overs.empty())
      cout << "None" << endl;
    else {
      for (overlap over : overs) {
        cout << over << endl;
      }
    }
    cout << "--" << endl;
    ct++;
  }



}
