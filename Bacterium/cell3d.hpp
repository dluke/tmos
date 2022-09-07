
#ifndef __CELL3D_HPP__
#define __CELL3D_HPP__ 

#include <unordered_set>
#include <utility>
#include <iterator>
#include <algorithm>
#include <vector>
#include <string>
#include <memory>
#include <tuple>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <Eigen/Dense>

#include "sample.hpp"
#include "vector3d.hpp"
#include "plane.hpp"
#include "sphere.hpp"
#include "potentials.hpp"
#include "pili.hpp"
#include "cell.hpp"

#include "piligrad.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::string;


class Cell3d: public ACell
{
  public:

  Plane* surface;
  Capsule body;

  // make private
  boost::function<double(double)> sforce;
  boost::function<double(double)> senergy;
  // Only reason to store this would be for passing to python as a wrapped type?
  std::shared_ptr<Spotential> sp;

  void _setup_sforce() {
    this->sp = std::make_shared<Spotential>(eps, get_body().R - attract_width, get_body().R);
    if (repulsive_only) {
      this->sforce = boost::bind(&Spotential::wca_force, this->sp, _1);
      this->senergy = boost::bind(&Spotential::wca_energy, this->sp, _1);
    }
    else {
      this->sforce = boost::bind(&Spotential::sforce, this->sp, _1);
      this->senergy = boost::bind(&Spotential::senergy, this->sp, _1);
    }
  }

  // getters
  Capsule& get_body();
  Plane* get_surface() { return surface; } 

/*  [> debugging surface init <]*/
  Cell3d(int idx, int npili)
    : ACell(idx, npili) {}

  // Empty cell object which we can then draw to a vtk file
  Cell3d(Capsule body) : ACell(0), body(body) {}


  Cell3d(int idx, int npili, Capsule body, Plane* surface) 
    : ACell(idx, npili), body(body), surface(surface)
  {
    ACell::pili_min_length = Pili::inside_length - 2*Cell3d::attract_width;
    _setup_sforce();
    init_anchors();
  }

  // Copy constructor
  Cell3d(const Cell3d& other) : ACell(other) {
    // check how often this is called
    _setup_sforce();
    // python constructed object. Just copy and keep reference.
    surface = other.surface; 
    body = other.body;
  }
  // virtual copy but populate the anchors and axisEqs vectors with data that we don't want to recompute
  // at the analysis stage
  std::unique_ptr<ACell> clone() override { 
    auto thisc = std::make_unique<Cell3d>(*this);
    return thisc;
  }
  Cell3d& operator=( Cell3d tmp);
  friend void swap(Cell3d& lhs, Cell3d& rhs);
  
  //~Cell3d(void) { delete surface; }
  
  void common_init() override {
    ACell::common_init(); 
  }

  void reset_surface() { surface->reset(); } 
  int get_num_contacts() { return this->get_surface()->get_num_contacts(get_body()); }

  Vector3d get_lab_anchor(Pili& pilus) {
    return body.frame.to_lab(pilus.anchor);
  }
  Vector3d get_lab_axisEq(Pili& pilus) {
    return body.frame.to_lab_rt(pilus.axisEq);
  }
  vector<Vector3d> get_bound_anchors();
  
  // want to update the bound pili anchors after updating the cell position
  void update_anchors() override {
    bool boundonly=true;
    for (std::shared_ptr<Pili> pilus : pili) {
      // always update anchors unless boundonly=true and isbound = false
      if ( !(boundonly && !pilus->isbound) ) {
        pilus->lab_axisEq = get_lab_axisEq(*pilus);
        pilus->lab_anchor = get_lab_anchor(*pilus);
      }
    }
  }

  Vector3d real_pos() override;
  Vector3d get_headpt() override;
  Vector3d get_trail() override;
  Vector3d real_centre() override;
  Vector3d get_vector() override;
  Vector3d get_axis() override;
  double get_length() override { return this->body.length; }

  anchorbits pili_body_distrib(); 
  // for testing with python
  Vector3d pili_body_position(); 
  
  // create pili3d object
  std::shared_ptr<Pili> spawn_pilus() override;

  CellEvent* attach(Pili& pilus, double tme) override;
  CellEvent* detach(Pili& pilus, double tme) override;
  
  vector<fpt> pili_force() override;
  vector<fpt> surface_force();
  fpt total_ft();

  Eigen::RowVectorXf surface_grad();
  Eigen::RowVectorXf pili_grad();
  Eigen::RowVectorXf grad();
  void set_state(Eigen::RowVectorXf state);
  Eigen::RowVectorXf get_state(void);
  double state_energy(Eigen::RowVectorXf state);
  Eigen::RowVectorXf state_gradient(Eigen::RowVectorXf state);
  
  void update(Vector3d dxy, Matrix3d deltaR) override;
  // for MD
  void update_pos(Vector3d dxy) override; 
  void update_axis(Matrix3d deltaR) override;
  // for line search
  void set_origin(Vector3d orv) { this->body.set_origin(orv); }
  void set_rotation(Vector3d p);

  void update_pili();
  double orthogonalise() { this->body.frame.orthogonalise(); } 

  // output utilities

  std::string report_touching() { return get_surface()->report_touching(); }
  std::string __str__() override;

  double energy_surface();
  double energy();
  
  // tune the strenth of the surface interaction
  // surface interaction stength is a static cell property
  static double eps; 
  static bool repulsive_only; 
  static double attract_width; 
  static bool cost_anchor_intersection;
  static string distrib_type;

};

// This subclass manages a reference to WLCgenerator
class CellWLC3d: public Cell3d
{
  public:

  std::shared_ptr<WLCgeneratorBase> wlc;
  std::shared_ptr<WLCgeneratorBase> get_wlc() { return wlc; }

  using Cell3d::add_pilus;
  using Cell3d::common_init;



  CellWLC3d(int idx, int npili, Capsule body, Plane* surface,
      std::shared_ptr<WLCgeneratorBase> wlc) 
    : Cell3d(idx, npili, body, surface), wlc(wlc)
  {}

  std::shared_ptr<Pili> spawn_pilus() override;

  // Copy constructor
  CellWLC3d(const CellWLC3d& other) : Cell3d(other) {
    // python constructed object. Just copy and keep reference.
    wlc = other.wlc;
  }
  /* clone using copy constructor */
  virtual std::unique_ptr<ACell> clone() override 
  { return std::make_unique<CellWLC3d>(*this); } 
  CellWLC3d& operator=( CellWLC3d tmp);
  friend void swap(Cell3d& lhs, Cell3d& rhs);


};

#endif

