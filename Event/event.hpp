

#ifndef __EVENT_HPP__
#define __EVENT_HPP__


// MOVED TO Base/

// Want event handling system which
// (i) returns a full set of information about the object at an instant of time
// (ii) contains any additional information related to the specific event
// (iii) Cell3dEvent should recursively gather information from the subobjects
//  of the cell which also have an event tracker (Pili)
// (iv) expose a python method to get() events

// Abstract base class for Events 
template <class T>
class Event 
{
  public: 
  Event(double time, std::unique_ptr<T> data) 
    : time(time), data(data) {}

  virtual std::unique_ptr<T> get_data() { return data; }
  virtual double get_time() { return time; }

  std::unique_ptr<T> data;
  double time;
};





#endif
