
#include "curvestep.hpp"

#include <boost/format.hpp>

#include <iostream>

using std::cout;
using std::endl;


boost::optional<Vector3d> Polygon::intersects(Lvec& lv)
{
  // find the intersection with home plane of polygon if it exists
  boost::optional<Vector3d> inter = AnyPlane::intersects(lv);
  if (inter && this->contains(*inter)) {
    return inter;
  }
  return boost::none;
}

vector<overlap> Polygon::overlap_vector(Capsule& body) 
{
}

  


//https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
template int sgn<double>(double val); // explicitely create sgn(double val) function

int _crosssign_ez(Vector3d v1, Vector3d v2) {
  // return the sign of (v1.cross(v2)).dot(e_z)
  return sgn(v1.x*v2.y - v1.y*v2.x);
}

bool ZPolygon::contains(Vector3d pt) 
{
  // iterate through the 'spokes' v1 - pt, v2 - pt ...
  Vector3d spa, spb;
  spb = this->verts.back() - pt; // previous spoke
  for (Vector3d va : this->verts) {
    spa = va - pt; // next anticlockwise spoke
    if (_crosssign_ez(spb, spa) < 0) {
      return false;
    }
    spb = spa;
  }
  return true;
}


//////////////////////////////////////////////////////////////////////////////////
// CurveElement


boost::optional<Vector3d> CurveElement::intersects(Lvec& lv) 
{
  // assume that atleast one end of lv can be projected onto this surface by e_z
  assert(this->x_in_range(lv.get_origin().x) || this->x_in_range(lv.get_endpt().x));

  // above surface check is done in higher level container
  Vector3d pta = lv.get_origin();
  Vector3d ptb = lv.get_endpt();
  //if (pta.z > _zmax && ptb.z > _zmax) {
    //return boost::none; // above surface, exit early 
  /*}*/

  double x1, x2;
  x1 = pta.x; x2 = ptb.x;
  int direction = (x2 >= x1) ? 1 : -1; // direction of pilus relative to e_x
  if (x1 > x2) { std::swap(x1, x2); } // ascending order
  // 
  // find bounds
  vector<double>::iterator lower,upper;
  // returns forward iterator
  lower = lower_bound(panex.begin(), panex.end(), x1);
  upper = upper_bound(panex.begin(), panex.end(), x2);

  //cout << "lower " << lower-panex.begin() << endl;
  //cout << "upper " << upper-panex.begin() << endl;

  // construct iterator which acts in the direction of pilus
  vector<PartialPlane>::iterator it, target;
  it = panes.begin(); target = panes.begin();
  it += lower-panex.begin()-1;
  it = std::max(panes.begin(), it);
  target = it + (upper-lower) + 1;
  target = std::min(panes.end(), target);

  //cout << "upper - lower " << upper-lower << endl;
  //cout << "pre target " << target-panes.begin() << endl;

  if (direction == 1) {
    // pass
  }
  else {
    std::swap(it, target);
    it--;
    target--;
  }

  //cout << "direction " << direction << endl;
  //cout << "it idx " << it-panes.begin() << endl;
  //cout << "target idx " << target-panes.begin() << endl;

  int idx;
  for (it; it != target; it+=direction)
  {
    idx = it-panes.begin(); 
    PartialPlane& pane = *it;
    boost::optional<Vector3d> inter = pane.intersects(lv);
    //cout << "checking intersection with pane " << idx << endl;
    if (inter) {
      //cout << "intersection is found with pane " << idx << endl;
      return inter;
    }
  }
  return boost::none;
}

// TODO improve efficiency
vector<overlap> CurveElement::overlap_vector(Capsule& body) 
{
  vector<overlap> overs;
  Lvec bodyax = body.get_lvec();
  // project bodyax onto xy plane at y=0
  Lvec xzbody = bodyax.xzproject();
  double R2 = body.R*body.R;

  // determine relevant panes of this segment

  // iterate over vertices
  // TODO it's not 100% clear that this correctly handles all edge cases
  double t, tt, m, d;
  for (int i = 0; i < this->panes.size(); ++i) {
    int cv = this->curvature[i];
    Vector3d vert = get_vert(i);
    PartialPlane pane = panes[i];

    // either vertex. project head/tail onto surface
    Vector3d edge = get_edge(i);
    // construct the line segment for this edge
    Lvec edgelv{vert, edge}; 
    Vector3d hpt = xzbody.get_origin(); 
    Vector3d tpt = xzbody.get_endpt(); 
    // distance vectors from head to surface segment
    double t = edgelv.closest_t(hpt);
    double tt = edgelv.closest_t(tpt);
    Vector3d vh = hpt - edgelv.get_pt(t);
    Vector3d vt = tpt - edgelv.get_pt(tt);
    if (vh.len2() < R2) {
      d = sqrt(vh.len2());
      overlap over = std::make_tuple(-body.length/2., d, vh/d);
      //cout << "contact with edge " << i << " and hpt " << hpt << endl;
      overs.push_back(over);
    }
    if (vt.len2() < R2) {
      d = sqrt(vt.len2());
      overlap over = std::make_tuple(body.length/2., d, vt/d);
      //cout << "contact with edge " << i << " and tpt " << tpt << endl;
      overs.push_back(over);
    }

    if (cv == 1) {
      // + vertex. also project vertex onto body line segment
      t = xzbody.closest_t(vert);
      // special case for end points t = 0. || t = 1. 
      // because end points are handled above
      if (t == 0. && pane.on_plane(hpt)) {
        continue;
      }
      if (t == 1. && pane.on_plane(tpt)) {
        continue;
      }

      // distance vector
      Vector3d vv  =  xzbody.get_pt(t) - vert;
      if (vv.len2() < R2) {
        // intersection
        // use t to obtain the point in 3d space
        m = body.get_m(t);
        d = sqrt(vv.len2());
        overlap over = std::make_tuple(m, d, vv/d);
        //cout << "contact with positive vertex " << i << endl;
        overs.push_back(over);
      }
    }
  }
  // return should be at most size 4 for this geometry, for suitable curvature
  return overs;
}


std::string CurveElement::__str__() {
  std::string form = "CurveElement: origin%s normal%s\n";
  std::string st = ( boost::format(form) % get_origin().__str__() % frame.e3.__str__() ).str();
  st += "Vertices\n";
  for (Vector3d vert : this->verts)  {
    st += vert.__str__() + "\n";
  }
  return st;
}

std::ostream& operator<<(std::ostream& os, CurveElement& celement)
{
  os << celement.__str__();
  return os;
}


//////////////////////////////////////////////////////////////////////////////////
// SegPlane


boost::optional<Vector3d> SegPlane::intersects(Lvec& lv) 
{
  // check if segment is above surface 
  Vector3d pta = lv.get_origin();
  Vector3d ptb = lv.get_endpt();
  if (pta.z > zmax && ptb.z > zmax) {
    return boost::none; // above surface, exit early 
  }
  // transform to integer coordinates
  int nstart = get_n(pta.x);
  int nend = get_n(ptb.x);
  // check intersections
  boost::optional<Vector3d> inter = _elements[nstart]->intersects(lv);
  if (inter) {
    return inter;
  }
  // handle edge case where lv crosses a repeating plane
  if (nstart != nend) 
  {
    boost::optional<Vector3d> inter = _elements[nend]->intersects(lv);
    if (inter)
      return inter;
  }
  // 
  return boost::none;
}

vector<overlap> SegPlane::overlap_vector(Capsule& body)
{
  // for now assert that xperiod < body.true_length()
  //
  // decide which repeating element to compare with
  // due to construction, edge case only happens for positive curvature 
  Lvec bodyax = body.get_lvec();
  double x1 = bodyax.get_origin().x;
  double x2 = bodyax.get_endpt().x;
  if (x1 > x2) {
    std::swap(x1,x2);
  }
  double xlextent = x1 - body.R;
  double xrextent = x2 + body.R;
  int nl = get_n(xlextent);
  int nr = get_n(xrextent);
  //cout << "nl  nr " << nl << " " << nr << endl;
  vector<overlap> l_overs = _elements[nl]->overlap_vector(body);
  //cout << "l touches " << l_overs.size() << endl;
  if (nl == nr) {
    return l_overs;
  }
  else {
    // covering two repeating elements
    vector<overlap> r_overs = _elements[nr]->overlap_vector(body);
    //cout << "r touches " << r_overs.size() << endl;
    l_overs.insert( l_overs.end(), r_overs.begin(), r_overs.end() );
  }
  return l_overs;
}


std::string SegPlane::__str__() {
  return "SegPlane object"; 
}
