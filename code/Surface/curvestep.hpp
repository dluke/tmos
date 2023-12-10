
#ifndef __CURVESTEP_HPP__
#define __CURVESTEP_HPP__

#include <vector>
#include <iterator>
#include <utility>
#include <map>

#include <boost/optional.hpp>

#include "vector3d.hpp"
#include "shape.hpp"
#include "plane.hpp"
#include "chain.hpp"
#include "sphere.hpp"

using std::cos;
using std::sin;

using std::vector;
using std::pair;

/*
   infinite sequence of curved steps formed from series of planes
 */


class Polygon: public AnyPlane
{
  /*
     Simple polygons
   */

  public:
  vector<Vector3d> verts;

  Polygon(Vector3d normal, vector<Vector3d> verts) : verts(verts) 
  {
    this->frame = Frame(Vector3d(0,0,0), normal);
    // always construct polygon with vertices in anticlockwise order
  }


  // geometry
  virtual bool contains(Vector3d pt) = 0;

  // Inherited from AnyPlane

  // REMINDER
  //typedef std::tuple<double, double, Vector3d> overlap;
  // distance along axis, m
  // overlap distance
  // unit vector in direction away from the surface

  virtual boost::optional<Vector3d> intersects(Lvec& lv) override;
  virtual vector<overlap> overlap_vector(Capsule& body) override;

};

class ZPolygon: public Polygon
{
  /*
     Polygon with normal in e_z direction
     put the long edge last 
  */
  public:
  using Polygon::Polygon;

  ZPolygon(vector<Vector3d> verts) : Polygon(e_z, verts) {}
  
  virtual bool contains(Vector3d pt) override;
  

};

class CurveElement: public AnyPlane
{
  // repeating segment for SegPlane, constructed from PartialPlanes
  double _planeangle, _panewidth;
  double _alpha, _th1; // the angles are [th1, planeangle, th1, -th1, ...]
  double _zmax;

  public:
  vector<PartialPlane > panes;
  vector<double> panex; // x coordinate of all the pane origins
  // store curvatures and vertices in xz plane for quick access
  const vector<int> curvature = {1,1,-1,-1,-1,1};
  vector<Vector3d> verts;
  vector<Vector3d> edges;
  int nseg = 6;

  CurveElement(Frame frame, double planeangle, double panewidth)
    : _planeangle(planeangle), _panewidth(panewidth)
  {
    this->frame = frame;
    _alpha = (3*M_PI - 2*planeangle)/3.;
    _th1 = M_PI/2. - _alpha/2.;
    _zmax = panewidth * (2*sin(_th1)  + sin(planeangle));
    this->_construct();
  }

  double get_alpha(void) { return this->_alpha; }
  double get_zmax(void) { return this->_zmax; }
  double get_xlength(void) { return get_lmpt().x - get_lstart().x; }
  // get the signs of the 
  vector<int> get_pane_csgn(int idx) { 
    return vector<int>{curvature[idx], curvature[idx+1]};
  }
  Vector3d get_pane_normal(int idx) { return panes[idx].frame.e3; }
  Vector3d get_vert(int idx) {
    if (idx == verts.size()) return get_lmpt(); 
    return verts[idx]; 
  }
  Vector3d get_edge(int idx) { return edges[idx]; }

  void _construct(void) {
    //
    panes.reserve(6); // 6 pane segment
    panex.reserve(6);
    verts.reserve(6);
    edges.reserve(6);
    this->_pane_construct(get_origin()+_zmax*e_z, -_th1);
    this->_pane_construct(panes.back().get_lmpt(), -_planeangle);
    this->_pane_construct(panes.back().get_lmpt(), -_th1);
    this->_pane_construct(panes.back().get_lmpt(), _th1);
    this->_pane_construct(panes.back().get_lmpt(), _planeangle);
    this->_pane_construct(panes.back().get_lmpt(), _th1);
    // let the end vertex belong to the next repeating Element
  }

  // convenience method for _construct
  void _pane_construct(Vector3d origin, double angle) {
    // infinite direction is e_y, repeating direction is e_x
    Vector3d e1 = Vector3d(cos(angle), 0, sin(angle));
    Vector3d normal = e1.cross(e_y);
    panes.push_back(PartialPlane(origin, normal, e1, 0, _panewidth));
    panex.push_back(origin.x);
    verts.push_back(origin);
    edges.push_back(_panewidth * e1);
  }

  Vector3d get_lmpt(void) { return panes.back().get_lmpt(); }
  Vector3d get_lstart(void) { return panes.front().get_origin(); }

  // check whether this repeating element contains coordinate x
  bool x_in_range(double x) { return (x >= get_lstart().x || x <= get_lmpt().x); }

  virtual boost::optional<Vector3d> intersects(Lvec& lv) override;
  virtual vector<overlap> overlap_vector(Capsule& body) override;

  virtual std::string __str__();

};

std::ostream& operator<<(std::ostream& os, CurveElement& celement);


class SegPlane: public AnyPlane
{

  double _planeangle, _panewidth;
  double _alpha;
  // curved plane constructed from a linearly repeating pattern of CurveElement
  // shared_ptr helps with python bindngs

  public:
  std::map<int,std::shared_ptr<CurveElement> > _elements;

  const int xnwidth = 10; // number of repeating elements in each direction
  double xperiod, zmax;

  SegPlane(Vector3d origin, Vector3d normal, double planeangle, double panewidth)
    : AnyPlane(origin, normal), _planeangle(planeangle), _panewidth(panewidth)
  {
    // create template element
    CurveElement tmpelement{Frame(), _planeangle, _panewidth};
    // get period
    this->xperiod = tmpelement.get_lmpt().x - tmpelement.get_origin().x; 
    this->zmax = tmpelement.get_zmax();

    // coordinate system is ready. construct the repeating elements.
    _plane_construct();
  }

  void _plane_construct() {
    // populate a mapping to hold the repeating elements
    for (int i = -xnwidth; i < xnwidth; ++i) {
      Frame fr{Vector3d(get_X(i),0,0), e_x, e_y, e_z};
      auto element = std::make_shared<CurveElement>(fr, _planeangle, _panewidth);
      _elements[i] = element;
    }
  } 

  // default constructor
  SegPlane() : SegPlane(Vector3d(), e_z, M_PI/6., 1.) {} 
  SegPlane(double planeangle, double panewidth) 
    : SegPlane(Vector3d(), e_z, planeangle, panewidth) {} 


  // corrdinate system mapping e1 coordinate to it's repeating element
  int get_n(double X) { return std::floor(X/xperiod); }
  double get_X(int n) { return n*xperiod; } 

  std::shared_ptr<CurveElement> get_element_at(int n) { return _elements[n]; }

  // geometry
  virtual boost::optional<Vector3d> intersects(Lvec& lv) override;
  virtual vector<overlap> overlap_vector(Capsule& body) override;


  virtual std::string __str__();


};

class CurveSteps: public AnyPlane
{

  double _height, _sep, _planeangle, _panewidth;

  double _alpha;

  public:
  CurveSteps(Vector3d origin, Vector3d normal, double height, double sep,
      double planeangle, double panewidth)
    : AnyPlane(origin, normal), _height(height), _sep(sep),
    _planeangle(planeangle), _panewidth(panewidth)
  {
    //double _alpha = (3*M_PI - 2*planeangle)/3.;
  }

  //_construct_chain_element() 

};




#endif
