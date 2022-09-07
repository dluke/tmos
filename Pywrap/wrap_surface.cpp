
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "vector3d.hpp"
#include "frame.hpp"
#include "shape.hpp"
#include "plane.hpp"
#include "cylinder.hpp"
#include "potentials.hpp"
#include "sphere.hpp"
#include "hexsurface.hpp"
#include "curvestep.hpp"
#include "periodic.hpp"

// tmp
//#include "cell.hpp"

namespace py = pybind11;

void def_surface(py::module &m) {

  auto msurface = m.def_submodule("surface");

  msurface.def("wca_force_abs", &wca_force_abs, "");
  msurface.def("wca_energy_abs", &wca_energy_abs, "");
  py::class_<Spotential> spotential(msurface, "Spotential");
  spotential
    .def(py::init<double,double,double>())
    .def("sforce", &Spotential::sforce, "")
    .def("senergy", &Spotential::senergy, "")
    ;

  py::class_<Shape> shape(msurface, "Shape");
  shape
    .def("get_frame", &Shape::get_frame, "")
    .def("get_origin", &Shape::get_origin, "")
    .def_readwrite("frame", &Shape::frame, "")
    ;

  //py::class_<Plane> plane(msurface, "Plane", shape);
  py::class_<Plane> plane(msurface, "Plane", shape);
  plane
    .def(py::init<>())
    .def(py::init<Vector3d&,Vector3d&>())
    .def("normal", &Plane::normal, "")
    .def("__str__", &Plane::__str__, "")

    .def_readwrite("frame", &Plane::frame, "")
      .def("get_origin", &Plane::get_origin, "")
      ;
    
  py::class_<AnyPlane> anyplane(msurface, "AnyPlane", plane);

  py::class_<Sphere> sphere(msurface, "Sphere", anyplane);
  sphere
    .def(py::init<>())
    .def(py::init<double&>())
    .def(py::init<Vector3d&,double&>())
    ;

  py::class_<Shell> shell(msurface, "Shell", sphere);
  shell
    .def(py::init<double&>())
    ;
  
  py::class_<NullPlane> nullplane(msurface, "NullPlane", plane);
  nullplane
    .def(py::init<>())
    ;

  py::class_<InfSteps> infsteps(msurface, "InfSteps", anyplane);
  infsteps
    .def(py::init<double&,double&,double&>())
    .def("report_touching", &InfSteps::report_touching, "")
    ;


  py::class_<HexGrid, std::shared_ptr<HexGrid> > hexgrid(msurface, "HexGrid", hexgrid);
  hexgrid
    .def(py::init<double&>())
    .def_readwrite("hrr", &HexGrid::hrr, "")
    .def_readwrite("hrq", &HexGrid::hrr, "")
    .def_readwrite("hrz", &HexGrid::hrr, "")
    ;

  py::class_<Hexc> hexc(msurface, "Hexc");
  hexc 
    .def(py::init<int&,int&>())
    .def("__getitem__", [](Hexc &h, int i) { 
        return h[i]; }, py::is_operator())
      ;

  py::class_<HexSphereGrid, Plane> hsphgrid(msurface, "HexSphereGrid");
  hsphgrid
    .def(py::init<double&>())
    .def(py::init<double&,Vector3d,Vector3d>())
    .def("get_R", &HexSphereGrid::get_R, "")
    .def("get_rq", &HexSphereGrid::get_rq, "")
    .def("get_xyz", &HexSphereGrid::get_xyz, "")
    .def("body_set", &HexSphereGrid::body_set, "")
    .def("coordinate_range", &HexSphereGrid::coordinate_range, "")
    .def("get_hrr", &HexSphereGrid::get_hrr, "")
    .def("get_hrq", &HexSphereGrid::get_hrq, "")
    .def("get_hrz", &HexSphereGrid::get_hrz, "")
    ;

  py::class_<SinePlane> sineplane(msurface, "SinePlane", anyplane);
  sineplane
    .def(py::init<>())
    .def(py::init<double,double>())
    .def(py::init<Frame,double,double>())
    .def(py::init<Frame,double,double,int>())
    .def("form", &SinePlane::form, "")
    .def("sp_form", &SinePlane::sp_form, "")
    .def("normal_form", &SinePlane::normal_form, "")
    //
    .def_readwrite("A", &SinePlane::A, "")
    .def_readwrite("B", &SinePlane::B, "")
    .def_readwrite("xperiod", &SinePlane::xperiod, "")
    .def_readwrite("invB", &SinePlane::invB, "")
    ;

  py::class_<SegPlane> segplane(msurface, "SegPlane", anyplane);
  segplane 
    .def(py::init<>())
    .def(py::init<double,double>())
    .def(py::init<Vector3d,Vector3d,double,double>())
    //
    .def_readwrite("zmax", &SegPlane::zmax, "")
    .def_readwrite("xperiod", &SegPlane::xperiod, "")
    ;


  py::class_<Capsule> (msurface, "Capsule", shape)
    .def(py::init<>())
    .def(py::init<Vector3d,Vector3d,double,double>())
    .def(py::init<Vector3d,Vector3d,Vector3d,double,double>())
    .def(py::init<Frame,double,double>())
    .def("alt_spherical", &Capsule::alt_spherical, "")
    .def("get_headpt", &Capsule::get_headpt, "")
    .def("get_endpt", &Capsule::get_endpt, "")
    .def("get_centerpt", &Capsule::get_centerpt, "")
    .def("get_axis", &Capsule::get_axis, "")
    .def("__str__", &Capsule::__str__, "")
    .def_readwrite("length", &Capsule::length, "")
    .def_readwrite("R", &Capsule::R, "")


    //for explicit control over state axisangle vector
    //.def("set_pstate", &Capsule::set_pstate, "")
    //.def("get_pstate", &Capsule::get_pstate, "")
    .def("init_Rk", &Capsule::init_Rk, "")
    ;

  py::class_<Chain, std::shared_ptr<Chain> > (msurface, "Chain")
    .def(py::init<>())
    .def_readwrite("axis", &Chain::axis, "")
    .def_readwrite("source", &Chain::source, "")
    .def_readwrite("a", &Chain::a, "")
    .def_readwrite("lines", &Chain::lines, "")
    .def_readwrite("targets", &Chain::targets, "")
    .def("compute_targets", &Chain::compute_targets, "")
    .def("__str__", &Chain::__str__, "")
    ;
}

