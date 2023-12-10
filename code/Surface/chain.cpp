

#include "chain.hpp"
#include <iostream>
#include <exception>
#include <assert.h>

#include <boost/format.hpp>

using std::cout;
using std::endl;

void Chain::set_source(Vector3d r) {
  this->source = r;
}

void Chain::set_lasta(double la) {
  double lb = lines.back().len();
  lines.back() = lines.back() * (la/lb);
}

void Chain::translate(Vector3d r) {
  this->source += r;
  if (this->have_targets) {
    for (int i = 0; i < targets.size(); i++)
    {
      targets[i] += r;
    }
  }
}

void Chain::rotate(Matrix3d M) {
  Vector3d psource = source;
  // for efficiency follow the build -> rotate -> translate -> compute_targets
  assert(source == Vector3d()); 
  // axis = M * axis;
  //this->translate(-source);
  for (int i = 0; i != lines.size(); i++)
  {
    lines[i] = M * lines[i];
  }
  //this->translate(psource);
}

void Chain::compute_targets() {
  if (blank) { throw std::runtime_error("Tried to compute targets of a blank Chain object"); }

  targets.resize(lines.size());
  targets[0] = source + lines[0];
  for (int i = 1; i < lines.size(); i++)
  {
    targets[i] = targets[i-1] + lines[i];
  }
  this->have_targets = true;
}


Lvec Chain::get_lvec(int n) { 
  assert(this->have_targets);
  if (n == 0) {
    return Lvec(source, lines[0]);
  }
  return Lvec(targets[n-1], lines[n]); 
}


// a iterator type structure for returning components might be better, todo
vector<Lvec> Chain::line_components()
{
  vector<Lvec> components;
  components.reserve(lines.size());
  for (int i = 0; i != lines.size(); i++)
  {
    components.push_back(this->get_lvec(i));
  }
  return components;
}

// slow
double Chain::get_length(void) 
{
  double ll = 0;
  for (auto l : lines) 
    ll += l.len();
  return ll;
}

// return the new length of the chain
double Chain::shorten(int ln, double lt) {
  // ln is the index of the intersecting segment
  double la = inter_t * lines[inter_n].len();
  double new_length = 0;
  if (ln > 0) {
      lines.resize(ln + 1);
      new_length +=  ln * this->a;
  }
  set_lasta(la);
  new_length += la;
  return new_length;
}

double Chain::len_to_inter(void)
{
  double lti = inter_n * this->a + inter_t * get_lvec(inter_n).len();
  return lti;
}


std::string Chain::__str__() 
{
  if (!blank) {
    std::string head = ( boost::format("Chain: source(%s)\naxis%s\n") % source % axis).str();
    head += "---lines---\n";
    for (auto& ll : lines) {
      head += ll.__str__() + "\n";
    }
    if (have_targets) {
      head += "---targets---\n";
      for (auto &vv : targets) {
        head += vv.__str__() + "\n";
      }
    }
    return head;
  }
  else {
    return "Chain: blank -- \n";
  }
}
