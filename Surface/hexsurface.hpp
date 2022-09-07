

#ifndef __HEXSURFACE_HPP__
#define __HEXSURFACE_HPP__

#include <cmath>
#include <string>
#include <functional>
#include <unordered_map>

#include <Eigen/Dense>

#include "point.hpp"
#include "vector3d.hpp"
#include "shape.hpp"
#include "plane.hpp"
#include "lvec.hpp"
#include "sphere.hpp"
#include "cell.hpp"

#include <boost/optional.hpp>
#include <boost/functional/hash.hpp>

using std::sqrt;
using std::sin;
using std::cos;
using std::abs;
using std::round;
using std::vector;

/* Define coordinate objects. Axial and Cube coordinates */
//https://www.redblobgames.com/grids/hexagons/#rounding

//class representing a simple 2d hexagon and define a contains method
class Hex
{
  Vector3d origin;
  double R;
};

// need forward dec?
class Cubec;

class Hexc: public Point2d
{
  public:
  using Point2d::Point2d;

  int len() const;
  Hexc unit() const;
  Cubec cube() const;

  Hexc operator-() const;

  // define hash function
  struct hashf {
    std::size_t operator()(Hexc const& hx) const {
      std::size_t seed = 0;
      boost::hash_combine(seed, hx[0]);
      boost::hash_combine(seed, hx[1]);
      return seed;
    }
  };

  std::string __str__() const;

};

class Cubec: public Point3d
{
  public:
  using Point3d::Point3d;
  Hexc hex() const;
};


/* define a set using slow O(N) insertion and operator== using non-member 
   function.
   unique_append
   --
*/
/* START vector set dec*/
void unique_append(vector<Hexc>& vhx, Hexc newhx);
void unique_extend(vector<Hexc>& vhx, vector<Hexc> extra);

/* END vector set dec*/



// Infinite geometric Line
class uHexLine
{
  public:

  struct Lequation {
    Lequation(int vari, int constv) : vari(vari), constv(constv) {}
    int vari;
    int constv;
  };

  Hexc source, line;

  uHexLine(Hexc source, Hexc line) : source(source), line(line) 
  { assert( line.len() == 1 ); }

  // 
  bool parallel( uHexLine const other ); 
  // check point on line
  bool coincide(Hexc const hp);
  // check line line coincidence
  bool operator==( uHexLine const other );

  // lines in hex space can be defined by x = c, y = c or z = c
  // Lequation defines an equation [x,y,z] = const
  Lequation get_line_eq() const;
  
  // return the intersection of two lines or boost::none if they are parallel
  // Assume that the lines are not coincident
  boost::optional<Hexc> intersect(uHexLine const hl);

};

// Half Space
class Hspace: public uHexLine
{
  struct Hequation: public uHexLine::Lequation {
    Hequation(int vari, int constv, std::function<bool(int,int)> comp) 
      : Lequation(vari, constv), g_or_l(comp) 
    {}
    std::function<bool(int, int)> g_or_l;
    bool operator()(int k) { return g_or_l(k, constv); }
    bool operator()(Cubec cube) { return g_or_l(cube[vari], constv); }
  };

  Hspace(Hexc source, Hexc line, std::function<bool(int,int)> g_or_l) 
    : uHexLine(source, line), g_or_l(g_or_l) {}
  std::function<bool(int,int)> g_or_l;

  Hequation get_space_eq();

};

// Line Segment
class HexLineSeg: public uHexLine
{
  public:
  HexLineSeg(Hexc source, Hexc line) : uHexLine(source, line) {}

  int len() { return line.len(); }
};


// 
const double r3 = sqrt(3.);

// define axial coordinate system
//https://www.redblobgames.com/grids/hexagons/
class HexGrid {
  public:

  static const Hexc origin;
  
// hr + hq + hz = 0
  static const Hexc hr;
  static const Hexc hq;
  static const Hexc hz;

  // define positive and negative directions
  static std::unordered_map<Hexc,Hexc,Hexc::hashf> absmap; 

  // define the basis according to a mapping of princple directions
  static std::unordered_map<int,Hexc> basis;
  static std::unordered_map<Hexc,int,Hexc::hashf> basis_ord;

//
  Vector3d hrr;
  Vector3d hrq;
  Vector3d hrz;

  // radius of a circle inscribed inside a hex
  double R;
  HexGrid(double R) : R(R) {
    // define the real space vectors
    hrr = 2*R * Vector3d(r3/2., 1/2., 0.);
    //hrr = 2*R * Vector3d(r3/2., 1/2., 0.);
    hrq = 2*R * Vector3d(0., 1., 0.);
    hrz = -hrr-hrq; // used?
  }

  // integer distance in hexgrid space
  static int hex_distance(Hexc hpa, Hexc hpb);
  // convert hex coordinate into real space coordinate
  Vector3d get_xyz(Hexc hp);
  // convert real space coordinate into hex coordinate
  Hexc get_rq(Vector3d pt);

  // Get the set of hexes under a line
  vector<Hexc> line(Lvec lv);

  // get a list of hexes in a range around target size N
  vector<Hexc> coordinate_range(Hexc target, int N);

  // Get the set of hexes in a polygon
  vector<Hexc> poly(vector<uHexLine> hlines);

};

// Specfic methods for determining how bacterium geometry and packed sphere surface interact
// remove HexGrid has parent and encapsulate instead
class HexSphereGrid: public Plane
{
  public:
  HexSphereGrid(double R, Vector3d origin, Vector3d normal)
    : Plane(origin, normal)
  {
    baseplane = Plane(origin, normal);
    this->hgrid = new HexGrid(R);
    have_overlaps = 0;
  }
  HexSphereGrid(double R) : HexSphereGrid(R, Vector3d(0,0,0), Vector3d(0,0,1)) {}

  // store a reference to the base plane object. It is useful for algorithms.
  Plane baseplane;
  HexGrid* hgrid;
  // use for intermediate calculations
  vector<Sphere*>  closesph;
  vector<overlap> curr_overlaps;
  bool have_overlaps;

  Sphere* get_sphere(Hexc hx);
  vector<Sphere*> get_closesph(Capsule body);
  void store_closesph(Capsule body);
  void reset() override {
    // delete any Spher
    for (Sphere* sph : this->closesph) {
      delete sph;
    }
    closesph.clear();
  }

  Vector3d get_hrr(void) { return hgrid->hrr; }
  Vector3d get_hrq(void) { return hgrid->hrq; }
  Vector3d get_hrz(void) { return hgrid->hrz; }


  // Get set of hexes which the Chain may be close to
  //vector<Hexc> chain_set(Chain ch);

  // Get the set of hexes under the capsule body 
  vector<Hexc> body_set(Capsule caps);

  virtual boost::optional<Vector3d> intersects(Lvec& lv)
  { throw std::runtime_error("Not Implemented"); }
  virtual boost::optional<Vector3d> intersects(Chain& ch);

  // bacterium body -> surface collision
  virtual vector<overlap> overlap_vector(Capsule& body);
  Eigen::RowVectorXf surface_grad(Capsule& body, contactf contact) override; 

  /* Explicitly state which methods of superclass we have not bothered to implemenet */
  virtual Vector3d normal(Vector3d pt)  
  { throw std::runtime_error("Not Implemented"); }
  virtual double distance(Vector3d pt) 
  { throw std::runtime_error("Not Implemented"); }
  virtual Vector3d project(Vector3d pt)
  { throw std::runtime_error("Not Implemented"); }
  virtual Vector3d distance_vector(Vector3d pt) 
  { throw std::runtime_error("Not Implemented"); }

  virtual int get_num_contacts(Capsule& body) override;
  virtual bool contains(Vector3d pt) override;

  // getter
  double get_R() { return this->hgrid->R; }

  // Hexgrid methods interface here because I can't subclass HexGrid 
  Vector3d get_xyz(Hexc hp) { return this->hgrid->get_xyz(hp); }
  Hexc get_rq(Vector3d pt) { return this->hgrid->get_rq(pt); }
  vector<Hexc> line(Lvec lv) { return this->hgrid->line(lv); }
  vector<Hexc> coordinate_range(Hexc target, int N) 
  { return this->hgrid->coordinate_range(target, N); }
  vector<Hexc> poly(vector<uHexLine> hlines) 
  { return this->hgrid->poly(hlines); }

  ~HexSphereGrid(void) {
    delete this->hgrid;
    this->reset();
  }
  
  protected:
  double _surface_energy(ACell& cell, contactf energyf) override;

};

#endif
