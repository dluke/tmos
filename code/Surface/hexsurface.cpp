

#include "hexsurface.hpp"

#include <iostream>
#include <tuple>

using std::get;
using std::min;
using std::max;

//////////////////////////////////////////////////////////////////////////////////
// coordinates

// Axial
int Hexc::len() const { return HexGrid::hex_distance(HexGrid::origin, *this); }
Hexc Hexc::unit() const { int ll = len(); return Hexc(operator[](0)/ll, operator[](1)/ll); }
Cubec Hexc::cube() const { 
  return Cubec((*this)[0], (*this)[1], -(*this)[0]-(*this)[1]); 
}
Hexc Hexc::operator-() const {
  return Hexc(-operator[](0), -operator[](1));
}

// Cube
Hexc Cubec::hex() const { return Hexc((*this)[0], (*this)[1]); }

// Hexc
std::string Hexc::__str__() const {
  return ( boost::format("(%d,%d)") % (*this)[0] % (*this)[1] ).str();
}

std::ostream& operator<<(std::ostream& os, Hexc& hx)
{
  os << hx.__str__();
  return os;
}

// utility
Hexc round_cube_to_axial(Vector3d cube) {
  double rx = round(cube.x);
  double ry = round(cube.y);
  double rz = round(cube.z);
  
  double x_diff = abs(rx - cube.x);
  double y_diff = abs(ry - cube.y);
  double z_diff = abs(rz - cube.z);

  if ( (x_diff > y_diff) && (x_diff > z_diff) )
    rx = -ry-rz;
  else if (y_diff > z_diff)
    ry = -rx-rz;
  else
    rz = -rx-ry;

  return Hexc(rx, ry);
}


//////////////////////////////////////////////////////////////////////////////////

void unique_append(vector<Hexc>& vhx, Hexc newhx) {
  bool isnew = true;
  for (Hexc hx : vhx) {
    if (hx == newhx) {
      isnew = false;
    }
  }
  if (isnew) {
    vhx.push_back(newhx);
  }
}
void unique_extend(vector<Hexc>& vhx, vector<Hexc> extra) {
  for (Hexc hx : extra) {
    unique_append(vhx, hx);
  }
}

//////////////////////////////////////////////////////////////////////////////////
// uHexLine

bool uHexLine::parallel( uHexLine const other ) {
  Hexc aline = HexGrid::absmap[line];
  Hexc oline = HexGrid::absmap[other.line];
  return (aline == oline);
}

bool uHexLine::coincide(Hexc const hp)
{
  Hexc aline = HexGrid::absmap[line];
  // same line if q == other.q
  if (aline == HexGrid::hr)
    return source[1] == hp[1];
  // same line if r == other.r
  else if (aline == HexGrid::hq)
    return source[0] == hp[0];
  // 
  else if (aline == HexGrid::hz)
    return ((source[0] - hp[0]) == (source[1] - hp[1]));
  else
    assert(false);
}

bool uHexLine::operator==( uHexLine const other )
{
  if (!parallel(other))
    return false;
  return this->coincide(other.source);
}

uHexLine::Lequation uHexLine::get_line_eq() const
{
  int hx = HexGrid::basis_ord[HexGrid::absmap[this->line]];
  return {hx, source.cube()[hx]};
}


boost::optional<Hexc> uHexLine::intersect(uHexLine const hl) 
{
  // temporarily assert that lines cannot be the same
  assert(!(*this == hl));
  // parrallel check
  if (parallel(hl))
    return boost::none;
  // intersection
  Lequation Lq = get_line_eq();
  Lequation hLq = hl.get_line_eq();
  int v1 = Lq.vari; int v2 = hLq.vari;
  if (v1 > v2) // fast sort two ints
    std::swap(v1, v2);
  vector<int> have{v1, v2};
  vector<int> all{0,1,2};
  vector<int> u(1);
  std::set_difference(all.begin(), all.end(), have.begin(), have.end(), u.begin());
  std::array<int, 3> cube;
  cube[Lq.vari] = Lq.constv;
  cube[hLq.vari] = hLq.constv;
  cube[u[0]] = -Lq.constv-hLq.constv;
  return Cubec(cube).hex();
}

/////////////////////////////////////////////////////////////////////////////////
// Half Space
Hspace::Hequation Hspace::get_space_eq()
{
  Lequation Lq = get_line_eq();
  return Hequation{Lq.vari, Lq.constv, this->g_or_l};
}



//////////////////////////////////////////////////////////////////////////////////
// HexGrid


const Hexc HexGrid::origin = Hexc(0,0);
const Hexc HexGrid::hr = Hexc(1,0);
const Hexc HexGrid::hq = Hexc(0,1);
const Hexc HexGrid::hz = Hexc(-1,-1); 

//static std::unordered_map<Hexc,Hexc,Hexc::hashf> absmap; 
std::unordered_map<Hexc,Hexc,Hexc::hashf> HexGrid::absmap =
{ {hr,hr},{-hr,hr},{hq,hq},{-hq,hq},{hz,hz},{-hz,hz} };

std::unordered_map<int,Hexc> HexGrid::basis({ {0,hr},{1,hq},{2,hz} });
std::unordered_map<Hexc,int,Hexc::hashf> HexGrid::basis_ord({ {hr,0},{hq,1},{hz,2} });

// manhattan distance on cube coordinates
//https://www.redblobgames.com/grids/hexagons/
int HexGrid::hex_distance(Hexc hpa, Hexc hpb) {
  return (abs(hpa[0] - hpb[0]) + abs(hpa[1] - hpb[1])
      + abs(-hpa[0]-hpa[1] + hpb[0]+hpb[1]))/2.;
}


Vector3d HexGrid::get_xyz(Hexc hp) {
  return hp[0] * hrr + hp[1] * hrq;
}

Hexc HexGrid::get_rq(Vector3d pt) {
  double var_r = 1./(2*R) * 2./r3 * pt.x;
  double var_q = 1./(2*R) * (-pt.x/r3 + pt.y);
  return round_cube_to_axial(Vector3d(var_r, var_q, -var_r-var_q));
}

// Get the set of hexes under a line
vector<Hexc> HexGrid::line(Lvec lv) {
  Hexc lsource = get_rq(lv.source);
  Hexc lend = get_rq(lv.get_endpt());
  int hxd = HexGrid::hex_distance(lsource, lend);
  vector<Hexc> ret;
  ret.push_back(lsource);
  if (hxd == 0) {
    assert(lsource == lend);
    return ret; // lsource = lend
  }
  if (hxd == 1) {
    vector<Hexc> adj {lsource, lend};
    return adj;
  }
  vector<Vector3d> linepts = lv.interp(hxd);
  for (vector<Vector3d>::iterator it = linepts.begin()+1; it != linepts.end()-1; it++)
  {
    // not guaranteed unique
    ret.push_back(get_rq(*it));
  }
  // finally add the endpoint
  ret.push_back(lend);
  return ret;
}


vector<Hexc> HexGrid::coordinate_range(Hexc target, int N)
{
  vector<Hexc> ret;
  Cubec ctarget = target.cube();
  int rmin = ctarget[0] - N;
  int rmax = ctarget[0] + N;
  int qmin = ctarget[1] - N;
  int qmax = ctarget[1] + N;
  int zmin = ctarget[2] - N;
  int zmax = ctarget[2] + N;
  for (int r = rmin; r < rmax+1; ++r) {
    for (int q = max(qmin, -r-zmax); q < min(qmax, -r-zmin)+1; ++q) {
      ret.push_back(Hexc(r,q));
    }
  }
  return ret;
}

// Get the set of hexes in a polygon
vector<Hexc> HexGrid::poly(vector<uHexLine> hlines) {
}


//////////////////////////////////////////////////////////////////////////////////
// HexSphereGrid

Sphere* HexSphereGrid::get_sphere(Hexc hx) {
  return new Sphere(this->get_xyz(hx), this->get_R());
}

// used for attachment of pili
boost::optional<Vector3d> HexSphereGrid::intersects(Chain& ch) {
  // we can immediately discount chain segments which are above the surface
  double minheight = get_origin().z + this->get_R();
  // declare the output variable
  boost::optional<Vector3d> inter; 
  for (Lvec lv : ch.line_components()) {
    if ( (abs(lv.get_origin().z) < minheight) || (abs(lv.get_endpt().z) < minheight) ) {
      // line segment not above the surface 
      for (Hexc hx : this->line(lv.xyproject())) { // iterate over close spheres
        // check intersection against these spheres
        Sphere*  sph = this->get_sphere(hx);
        inter = sph->intersects(lv);
        if (inter && inter->z > 0.) { // need to set orientation of surface, todo
          // If our intersection is below the baseplane then discard
          return inter;
        }
        delete sph;
      }
      // check intersection with the planar part of the surface
      inter = this->baseplane.intersects(lv);
      if (inter)
        return inter;
    }
  }
  return boost::none;
}


//# re-write this method

vector<Hexc> HexSphereGrid::body_set(Capsule body) {
  // need HexGrid::poly implementation to do this efficiently
  // Lets quickly hack out a method which just uses HexGrid::coordinate_range

  // tmp hack
  int reasonable_range = (int)std::ceil( body.R/this->get_R() );
  // get range for midpoint and both ends
  Hexc hhex = this->get_rq(body.get_headpt());
  Hexc chex = this->get_rq(body.get_centerpt());
  Hexc endhex = this->get_rq(body.get_endpt());

  // initialise return array with the first range
  vector<Hexc> ret = this->coordinate_range(hhex, reasonable_range);
  unique_extend(ret, this->coordinate_range(chex, reasonable_range));
  unique_extend(ret, this->coordinate_range(endhex, reasonable_range));
  return ret;
}


vector<Sphere*> HexSphereGrid::get_closesph(Capsule body)
{
  if (this->closesph.empty()) {
    this->store_closesph(body);
  }
  return this->closesph;
}

void HexSphereGrid::store_closesph(Capsule body) {
  vector<Hexc> hxset = this->body_set(body);
  for (Hexc hx : hxset) {
    Sphere* sph = this->get_sphere(hx);
    this->closesph.push_back(sph);
  }
}


vector<overlap> HexSphereGrid::overlap_vector(Capsule& body) {
  vector<overlap> ret;
  for (Sphere* sph : this->get_closesph(body)) {
    boost::optional<overlap> inter = body.intersects(*sph);
    if (inter) {
      ret.push_back(*inter);
    }
  }
  curr_overlaps = ret;
  have_overlaps = 1;

  // calling reset here defeats the purpose of storing results
  this->reset();
  return ret;
}

int HexSphereGrid::get_num_contacts(Capsule& body) {
  return curr_overlaps.size();
}

bool HexSphereGrid::contains(Vector3d pt) {
  return (baseplane.contains(pt) || this->get_sphere(this->get_rq(pt))->contains(pt));
}


double HexSphereGrid::_surface_energy(ACell& cell, contactf energyf) 
{
  Capsule& body = cell.get_body();
  double energy = 0;
  vector<overlap> inter = this->overlap_vector(body);
  for (overlap over : inter) {
    double r = get<1>(over);
    energy += energyf(r);
  }
  // add energy for baseplane 
  // can set a cutoff on surface sphere radius to skip this calculation
  double base_energy = this->baseplane.surface_energy(cell, energyf);
  return base_energy + energy;
}

Eigen::RowVectorXf HexSphereGrid::surface_grad(Capsule& body,
    contactf contact) {
  Vector3d fgrad;
  double fg; 
  vector<double> pk = {0,0,0};

  Rk bodyRk = body.get_Rk();

  vector<overlap> inter = this->overlap_vector(body);
  for (overlap over : inter) {
    double m = get<0>(over);
    double r = get<1>(over);
    Vector3d rhat = get<2>(over);

    // the surface repulsive function, take (-) to get direction up the gradient
    fg = -contact(r);
    // really? todo
    fgrad += fg * rhat;
    Vector3d dLjde = fg * m * rhat;
    for (int k=0; k<3; k++) {
      pk[k] += dot(dLjde, bodyRk.Rkvez[k]);
    }
  }
  Eigen::RowVectorXf grad(6);
  grad << fgrad.x, fgrad.y, fgrad.z, pk[0], pk[1], pk[2];

  // check collision with baseplane
  Eigen::RowVectorXf basegrad = this->baseplane.surface_grad(body, contact);

  return grad + basegrad;
}

