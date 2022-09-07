

// We will check that the string output of each object is correct and
// correctly formatted by eye rather than write automated unit tests

#include <iostream>
#include "vector3d.hpp"
#include "frame.hpp"
#include "matrix3d.hpp"
#include "shape.hpp"
#include "plane.hpp"
#include "sphere.hpp"

#include "pili.hpp"
#include "cell.hpp"
#include "cell3d.hpp"
//#include "mdynamics.hpp"

using std::cout;
using std::endl;

int main() {

  cout << "Vector3d() class" << endl;
  cout << Vector3d() << endl;
  cout << "Vector3d() using __str__" << endl;
  cout << Vector3d().__str__() << endl;

  Capsule cap = Capsule(Vector3d(), e_x, 0.5, 2.);
  cout << "Capsule.__str__()" << endl;
  cout << cap.__str__() << endl;

  Plane pl = Plane(Vector3d(), e_z);
  cout << "Plane.__str__()" << endl;
  cout << pl.__str__() << endl;

  //Cell sicell = Cell(0, Vector3d(), Vector3d(1,0,0));
  Cell3d cell = Cell3d(0, 4, cap, pl);
  cout << "Cell3d.__str__()" << endl;
  cout << cell.__str__() << endl;


}
