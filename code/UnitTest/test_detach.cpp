
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

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

#include "plane.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::shared_ptr;

//# constants
double R = 0.5;
double l = 2.;

// global
int nmax = 10;

std::shared_ptr<WLCgeneratorBase> new_kp() {
  double T = 273+30;
  const double kb = 1.38064852 * (10e-23 * 10e18);
  double Lp = 5.;
  double ka = Lp * kb*T;
  double a = 1.0;
  auto KP = std::make_shared<KPgenerator>(ka,a,nmax,T);
  return KP;
}

TEST_CASE("Test WLC model") {

  KPgenerator KP = *(std::dynamic_pointer_cast<KPgenerator>(new_kp()));

  double angle_1 = KP.one_theta();
  double angle_2 = KP.one_theta();
  CHECK((angle_1 > 0. && angle_1 < M_PI/2));
  CHECK(angle_1 != angle_2);

  Eigen::VectorXd theta = *(KP.sample_theta());
  CHECK( theta.size() == nmax );
  CHECK(( theta[0] > 0. && theta[0] < M_PI/2 ));
  CHECK(( theta[nmax-1] > 0. && theta[nmax-1] < M_PI/2 ));

  auto phi = *(KP.sample_phi());
  CHECK( phi.size() == nmax );
  CHECK(( phi[0] > 0. && phi[0] < 2*M_PI ));
  CHECK(( phi[nmax-1] > 0. && phi[nmax-1] < 2*M_PI ));

}

TEST_CASE("detach from debug simulation") {

  Pili::inside_length = 0.04;
  Pili::d_bound = 0.004;
  Pili::d_free = 0.004;

  // auto wlc = std::make_shared<WLCgenerator>(-1, 1.0, 19);
  // wlc needs to be filled with data
  // ...
  std::shared_ptr<WLCgeneratorBase> KP = new_kp();

  auto *plane = new Plane();
  Capsule caps{Vector3d(0, 0, 0), Vector3d(1,0,0).unit(), R, l};
  auto cell = new CellWLC3d(0, 0, caps, plane, KP);
  cell->common_init();

  SUBCASE("") {
    CHECK( cell->pilivar != 0 );
    CHECK( cell->npili == 0 );
    CHECK( cell->pili_min_length > 0. );
    CHECK( cell->pili_min_length < Pili::inside_length - Cell3d::attract_width );
  }

  SUBCASE("spawn pili") {
    auto abits = cell->pili_body_distrib();
    CHECK( !isnan(abits.anchor.x) );
    CHECK( !isnan(abits.axis.x) );
    CHECK( abits.axis.z > 0 );
    for (int i = 0; i < 1000; i++) {
      auto abits = cell->pili_body_distrib();
      cout << abits.axis << endl;
    }
    for (int i = 0; i < 1000; i++) {
      std::shared_ptr<Pili> pilus = cell->spawn_pilus();
      cout << pilus->axisEq << endl;
    }
  }

  SUBCASE("pili counter") {
    // cout << pilus->__str__() << endl;
    std::shared_ptr<Pili> pilus = cell->spawn_pilus();
    cell->add_pilus(pilus);
    CHECK(cell->pili_counter == 1);
    CHECK(cell->pili.size() == 1);
    cell->dissolve_pilus(pilus);
    CHECK(cell->pili_counter == 1);
    CHECK(cell->pili.size() == 0);
  }

  SUBCASE("") {
    // construct pilus with specific geometry
  }

  //   std::shared_ptr<PiliWLC3d> pilus (new PiliWLC3d(pili_counter,
  //         pl, pl, vector3d(), vector3d(1,0,0), 0, 0, wlc));
  // }


}

