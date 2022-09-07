
#include <memory>
#include <cmath>
#include <iostream>
#include <vector>

#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "axisangle.hpp"

#include "wlc.hpp"
#include "pili.hpp"
#include "cell.hpp"
#include "cell3d.hpp"

#include "sphere.hpp"
#include "shape.hpp"
#include "plane.hpp"
#include "step.hpp"

#define BOOST_TEST_MODULE SphereTest
#include <boost/test/unit_test.hpp>

using std::cout;
using std::endl;
using std::vector;


//# constants
double R = 0.5;
double l = 2.;

// helper
// create a cell with an exactly determined pili configuration
Cell3d* default_cell3d() {
  Capsule caps{Vector3d(0, 0, 0), Vector3d(1,0,0).unit(), R, l};
  Plane* plane = new Plane();
  Cell3d* cell = new Cell3d(0, 4, caps, plane);
  cell->common_init();
  return cell;
}

BOOST_AUTO_TEST_CASE( cell_copy ) 
{
  ACell* cell = default_cell3d();
  ACell* clone = cell->clone();
  clone->idx = 1;
  BOOST_REQUIRE( cell->idx == 0 );
  BOOST_REQUIRE( clone->idx == 1 );
  clone->pili[0]->elongate();
  BOOST_REQUIRE( cell->pili[0]->leq == Pili::d_free);
  BOOST_REQUIRE( clone->pili[0]->leq == Pili::d_free + Pili::d_free );

}

BOOST_AUTO_TEST_CASE( cell3d_init )
{
  // init all types of 3d cell and check the typestr of their pili

  Capsule caps{Vector3d(0, 0, 0.0), Vector3d(1,0,0).unit(), R, l};
  Plane* plane = new Plane();
  std::shared_ptr<WLCgenerator> wlc = std::make_shared<WLCgenerator>(1., 1., 3);
   
  Cell3d* cell = new Cell3d(0, 4, caps, plane);
  cell->common_init();
  BOOST_REQUIRE( cell->pili[0]->get_typestr() == "rod" );
  BOOST_REQUIRE( cell->pili.size() == 4 );
  CellWLC3d* wlccell = new CellWLC3d(0, 4, caps, plane, wlc);
  wlccell->common_init();
  BOOST_REQUIRE( wlccell->pili[0]->get_typestr() == "wlc" );
  BOOST_REQUIRE( wlccell->pili.size() == 4 );

  // the above two step initialisation seems to work so why didn't calling common_init
  // in the respective constructs work -- todo check this later
}

BOOST_AUTO_TEST_CASE( cell_update_pos_ax )
{
  double tol = 0.000001;
  Cell3d cell = *default_cell3d();
  Vector3d x = cell.get_body().get_origin();
  Vector3d ax = cell.get_body().get_axis();
  Vector3d norv{1,1,1};
  cell.set_origin(norv);
  BOOST_REQUIRE( tolequals(cell.get_body().get_origin(), norv, tol) );
  Vector3d naxis{1,1,0};
  naxis.norm();
  cell.set_rotation( Rmatrix(e_z, naxis) );
  cell.set_rotation( Rmatrix(e_z, naxis) );
  BOOST_REQUIRE( tolequals(cell.get_body().get_axis(), naxis, tol) );
}

BOOST_AUTO_TEST_CASE( cell_surface_force )
{
  Cell3d::eps = 100.;
  Cell3d cell = *default_cell3d();

  //cout << cell.sforce(1. , 0.5) << endl;
  BOOST_REQUIRE( cell.sforce(1. , 0.5) == 0);
  BOOST_REQUIRE( cell.sforce(1. , 0.49) == 0);
  //cout << cell.sforce(1. , 0.495) << endl;
  //cout << cell.sforce(1. , 0.485) << endl;

}

BOOST_AUTO_TEST_CASE( cell_state_vector )
{
  Cell3d cell = *default_cell3d();
  cout << "Test cell state vector" << endl;
  Capsule& body = cell.get_body();
  Vector3d p = body.get_axisangle();
  Matrix3d Rp = aarmatrix(p);

  BOOST_REQUIRE( body.get_phi() == 0. );
  body.set_phi( M_PI/2. );
  BOOST_REQUIRE( body.get_phi() == M_PI/2. );
  
  //cout << "axis_angle p" << endl;
  body.set_rotation(Vector3d(0,0,M_PI));
  //cout << body.get_axisangle() << endl;
  BOOST_REQUIRE( body.get_axisangle().len() == 0. );

  // set phi again, axis angle is 0.
  body.set_phi( M_PI/2. );
  BOOST_REQUIRE( body.get_phi() == M_PI/2. );
}



