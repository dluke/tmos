

#include <string>
#include "shape.hpp"

#include <boost/format.hpp>

using std::get;

std::ostream& operator<<(std::ostream& os, const overlap& over)
{
  std::string form = "overlap: m %f\n r %f\n direction %s";
  return os << ( boost::format(form) 
      % get<0>(over)
      % get<1>(over)
      % get<2>(over).__str__()
        ).str();
}



// Chain intersection depends on Line Vector intersection method so this
// implementation can be rather generic

boost::optional<Vector3d> Shape::intersects(Chain& ch) {
  for (Lvec lv : ch.line_components() )
  {
    boost::optional<Vector3d> inter = this->intersects(lv);
    if (inter) {
      return inter;
    }
  }
  return boost::none;
}


std::string Shape::__str__() {
  return "Abstract Base Shape";
}
