

#include <iostream>
#include "event.hpp"
#include "cell.hpp"
#include "plane.hpp"
#include "sphere.hpp"
#include "cell3d.hpp"
#include "kmc.hpp"

#define BOOST_TEST_MODULE Event
#include <boost/test/unit_test.hpp>


using std::cout;
using std::endl;
using std::vector;


Cell3d* default_cell3d() {
  Capsule caps{Vector3d(0, 0, 0.0), Vector3d(1,0,0).unit(), 0.5, 2};
  std::shared_ptr<Plane> plane = std::make_shared<Plane>();
  Cell3d* cell = new Cell3d(0, 4, caps, plane);
  cell->common_init();
  return cell;
}

BOOST_AUTO_TEST_CASE( contstr )
{

  Cell3d* cell = default_cell3d();
  ACell* ccell = cell->clone();
  cell->update_anchors();
  CellEvent ev{0., 1, "string", cell};
  cout << ev.get_data()->anchors[0] << endl;
  cout << ccell->anchors[0] << endl;
}

