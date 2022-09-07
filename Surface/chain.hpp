
#ifndef __CHAIN_HPP__
#define __CHAIN_HPP__


#include <string>
#include <vector>
#include <Eigen/Dense>

#include "sample.hpp"
#include "frame.hpp"
#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "lvec.hpp"

using std::vector;

class Chain
{
  public:
  // the axis is that of the anchor 
  Vector3d axis;
  //
  Vector3d source;
  bool have_targets;
  // endpoints of each vector
  vector<Vector3d> targets;
  // sequence of vectors forming a chain
  vector<Vector3d> lines; 
  double a; // segment length
  bool blank;
  // intersection tracking
  double inter_t;
  int inter_n;


  Chain() { blank = true; }
  operator bool() const { return !this->blank; }

  Chain(vector<Vector3d> lines) : lines(lines) 
  {
    a = 1.;
    axis = lines[0].unit();
    source = Vector3d();
    have_targets = false;
    blank = false;
  }

  Chain(Vector3d source, Vector3d axis, double a, vector<Vector3d> lines) 
    : source(source), axis(axis), a(a), lines(lines)
  { blank = false; }

  Chain(Vector3d source, Vector3d axis, double a, vector<Vector3d>lines, vector<Vector3d> targets) 
    : source(source), axis(axis), a(a), lines(lines), targets(targets)
  { blank = false; have_targets = true; }

  // Chain with one segment and no bends
  Chain(Vector3d axis, double a) : a(a), axis(axis) { 
    lines.push_back(a * axis.unit());
    source = Vector3d();
    blank = false;
    have_targets = false;
  }

  // Take angles in constructor and rotate each of them to get a sequence 
  // of correctly oriented line vectors
  Chain(Vector3d axis, double a, Eigen::VectorXd thetas, Eigen::VectorXd phi) 
    : a(a), axis(axis)
  {
    lines.reserve(thetas.size()+1);
    //
    source = Vector3d();
    blank = false;
    have_targets = false;
    lines.push_back(a * axis.unit());

    // construct
    Matrix3d rstate = Rmatrix(e_z, axis);

    for (int i = 0; i < thetas.size(); i++)
    {
      // get the next line with random phi
      Vector3d v1 = Frame::spherical(a, thetas[i], phi[i] );
      Matrix3d rtmp = Rmatrix(e_z, v1);
      // multiply rotation matrices to keep track of the total rotation
      rstate = rtmp * rstate;
      Vector3d l1  = rstate * e_z;
      l1.norm();
      lines.push_back( a * l1 );
    }
  }

  Vector3d get_origin() { return source; }

  void set_source(Vector3d r);
  // set the length of the last line
  void set_lasta(double la);

  void translate(Vector3d r);
  void rotate(Matrix3d M);
  void compute_targets();

  double shorten(int ln, double la);
  // Only call after Surface->intersection
  double len_to_inter(void);
  //
  double get_length(void);
  
  Lvec get_lvec(int n);
  vector<Lvec> line_components();

  // copy 
  // isn't this basically the default copy constructor?
  Chain(const Chain& ch)
    : source(ch.source), axis(ch.axis), a(ch.a), lines(ch.lines)
  {
    if (ch.have_targets) {
      targets = ch.targets;
    }
    blank = ch.blank;
    have_targets = ch.have_targets;
  }

  std::string __str__();

};


#endif

