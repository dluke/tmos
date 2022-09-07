

#ifndef __KMC_HPP__
#define __KMC_HPP__

#include <cmath>
#include <random>
#include <iterator>
#include <algorithm>
#include <string>
#include <tuple>
#include <memory>

// r2() random number 
#include "sample.hpp"
#include "event.hpp"
#include "vector3d.hpp"
#include "pili.hpp"
#include "cell.hpp"
#include "cell3d.hpp"

/*
Object for managing Kinetic Monte Carlo simulation 

Keeps CellEvent object for relaying data to python
*/

class Clock
{
  public:
  Clock(double deltat) : deltat(deltat) {}

  bool now(double tme) {
    if (tme > prev_t + deltat) {
      prev_t = tme;
      return true;
    }
    return false;
  }
  
  private:
  double deltat;
  double prev_t = 0.;
};

typedef tuple<double, string, bool, int> kmctup;

// we update the simulation state by calling Kmc::kmcstep()
class Kmc
{
  public:
  Kmc(double dt) : clock(Clock(dt)) {}

  kmctup kmcstep(ACell& cell, double tme);
  void attachstep(ACell& cell, double tme);
  void postmd(ACell& cell, double tme);

  // use shared ptr becuase unique_ptr doesn't make much sense in Python (?)
  std::vector<std::shared_ptr<CellEvent> > events;
  
  private:
  Clock clock;
  
  void clear_events(void) { events.clear(); }
  void set_this_event(CellEvent* e) { events.push_back( std::shared_ptr<CellEvent>(e) ); }
  
  // std::unique_ptr<CellEvent> ev;

};


#endif
