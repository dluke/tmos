

#ifndef __CELL_HPP__
#define __CELL_HPP__


#include <iostream>
#include <utility>
#include <unordered_set>
#include <iterator>
#include <exception>
#include <vector>
#include <cmath>
#include <string>
#include <memory>

#include <boost/optional.hpp>
#include "sample.hpp"
#include "event.hpp"
#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "sphere.hpp"
#include "potentials.hpp"

#include "pili.hpp"

#include <boost/format.hpp>

using std::begin;
using std::end;
using std::unordered_set;
using std::tuple;
using std::string;

class ACell;

class CellEvent: public Event<ACell>
{
  public:
  CellEvent(std::unique_ptr<ACell> data) : Event<ACell>(0., std::move(data)) {}
  CellEvent(double tme, int pidx, string process, 
      std::unique_ptr<ACell> data, string trigger)
    : Event<ACell>(tme, std::move(data)), pidx(pidx), 
    process(process), trigger(trigger) 
    {}

  int pidx;  
  string process;
  // The reason for the event i.e. attach, release, spawn, istaut, isslack 
  string trigger;

  string __str__();

};


// reduced to a named pair of vectors
struct anchorbits {
  anchorbits(Vector3d anchor, Vector3d axis)
    : anchor(anchor), axis(axis) {}
  Vector3d anchor, axis;
};

//////////////////////////////////////////////////////////////////////////
class ACell // base class
{
  protected: 
  ACell(int idx, int npili) : idx(idx), npili(npili) 
  { /* for two step construction */
      pili_counter = 0;
  }
  ACell(int idx) : ACell(idx, 1)  {}
  
  // copy constructor
  protected: ACell(const ACell& other) {
    idx = other.idx;
    npili = other.npili;
    pili_counter = other.pili_counter;
    for (std::shared_ptr<Pili> p : other.pili) {
      pili.push_back(p->clone());
    }
    anchors = other.anchors;
    axisEqs = other.axisEqs;
  }
  //ACell& operator=( ACell tmp);

  public:

  friend void swap(ACell& lhs, ACell& rhs);
  virtual std::unique_ptr<ACell> clone() = 0;

  // used by surface module
  virtual vector<Vector3d> get_bound_anchors() = 0;

  // IMPLEMENTED
  virtual CellEvent* create_event();
  virtual void common_init();
  virtual vector<vector<double>> get_rates();
  virtual vector<double> get_cell_rates();
  virtual CellEvent* elongate(Pili& pilus, double tme);
  virtual CellEvent* shrink(Pili& pilus, double tme);
  virtual CellEvent* switch_se(Pili& pilus, double tme);
  virtual std::shared_ptr<Pili> spawn_pilus() = 0;

  // ABSTRACT
  virtual Vector3d real_pos() = 0;
  virtual Vector3d get_headpt() = 0;
  virtual Vector3d get_trail() = 0;
  virtual Vector3d real_centre() = 0;
  virtual Vector3d get_vector() = 0;
  virtual Vector3d get_axis() = 0;
  virtual double get_length() = 0;

  virtual void set_origin(Vector3d orv) = 0;
  //virtual void set_rotation(Matrix3d M) = 0;
  virtual void set_rotation(Vector3d p) = 0;

  //virtual void add_pilus(int i, double pl) = 0;
  virtual void add_pilus(std::shared_ptr<Pili> pilus);
  virtual void prep_dissolve_pilus(std::shared_ptr<Pili> pilus);
  virtual void dissolve_pilus(std::shared_ptr<Pili> pilus);
  virtual CellEvent* attach(Pili& pilus, double tme) = 0;
  virtual CellEvent* detach(Pili& pilus, double tme) = 0;

  virtual void update_pos(Vector3d dxy) = 0;
  virtual void update_pili() = 0;

  // These functions will be implemented by Cell3d subclass
  // But since all methods need an implementation, their default behaviour
  // should be to raise an exception.
  // Hopefully I will be able redesign the code at some point and get rid of this
  //
  // 2D methods
  virtual vector<fpt> pili_force() { throw std::runtime_error("Not Implemented"); } 
  virtual void update_axis(double dangle) { throw std::runtime_error("Not Implemented"); };
  virtual void update(Vector3d dxy, double dangle) { throw std::runtime_error("Not Implemented"); };
  virtual double pili_distrib() { throw std::runtime_error("Not Implemented"); };
  // 3D methods
  virtual vector<fpt> surface_force() { throw std::runtime_error("Not Implemented"); };
  virtual fpt total_ft() { throw std::runtime_error("Not Implemented"); };
  virtual fltq pili_stretch() { throw std::runtime_error("Not Implemented"); } // dep this eventually
  virtual void update_axis(Matrix3d deltaR) { throw std::runtime_error("Not Implemented"); };
  virtual void update(Vector3d dxy, Matrix3d deltaR) { throw std::runtime_error("Not Implemented"); };
  virtual anchorbits pili_body_distrib() 
  { throw std::runtime_error("Not Implemented"); };
  virtual double orthogonalise() { throw std::runtime_error("Not Implemented"); } 
  virtual double track_overlap() { throw std::runtime_error("Not Implemented"); } 
  virtual Capsule& get_body() { throw std::runtime_error("Not Implemented"); } 
  std::shared_ptr<Pili> get_pilus(int pidx);
  int get_pilus_vidx(std::shared_ptr<Pili> p);
  //
  
  // IMPLEMENTED BUT DON'T USE IN 2D 
  virtual Vector3d get_lab_anchor(Pili& pilus) { return pilus.anchor; }
  virtual Vector3d get_lab_axisEq(Pili& pilus) { return pilus.axisEq; }

  // MEMBER VARIABLES

  int idx, npili, pili_counter;
  std::shared_ptr<Pili> last_pilus_ref;
  vector<std::shared_ptr<Pili> > pili;

  // store the position of the anchors of the bound pili in lab frame for force calculation.
  vector<Vector3d> anchors, axisEqs;

  // END MEMBER VARIABLES

  // call once
  void init_anchors() {
    anchors.resize(pili.size());
    axisEqs.resize(pili.size());
  }
  virtual void update_anchors() = 0;

  virtual std::string __str__() = 0;

  // getters, setters
  void set_npili(int n) { this->npili = n; }
  
  // output utilities
  int nbound() { return std::count_if(pili.begin(), pili.end(), [](auto pilus){ return pilus->isbound; }); }
  int ntaut() { return std::count_if(pili.begin(), pili.end(), [](auto pilus){ return pilus->istaut(); }); }
  // return the number of pili with length > threshold
  int num_pili(double threshold) { 
    int count = 0;
    for (auto pilus: pili) {
      if (pilus->leq > threshold) {
        count += 1;
      }
    }
    return count;
  }
  double pbrf();
  
  double pl_avg();
  double l_total();
  virtual double energy_pili();
  virtual double energy();

  static int cycles_limit, spawn_extension_state;
  static bool pilus_replacement_on;
  static double k_spawn, pilivar, maxpl, polarisation, pili_min_length;
  static bool running_start;

};

std::ostream& operator<<(std::ostream& os, ACell& cell);



#endif
