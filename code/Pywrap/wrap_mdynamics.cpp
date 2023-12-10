

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "mdynamics.hpp"
namespace py = pybind11;


void def_mdynamics(py::module &m) 
{
  auto mmdynamics = m.def_submodule("mdynamics");

  py::class_<Summary> summary(mmdynamics, "Summary");
  summary
    .def_readwrite("nsteps", &Summary::nsteps, "")
    .def_readwrite("rms", &Summary::rms, "")
    .def_readwrite("prior_state", &Summary::prior_state, "")
    .def_readwrite("accepted_state", &Summary::accepted_state, "")
    .def_readwrite("prior_pilus", &Summary::prior_pilus,  "")
    .def_readwrite("accepted_pilus", &Summary::accepted_pilus,  "")
    ;

  py::class_<MDintegrate> mdi(mmdynamics, "MDintegrate");
  mdi 
    .def("reset", &MDintegrate::reset,"")
    .def("step", &MDintegrate::step,"")
    .def("get_x", &MDintegrate::get_x,"")
    .def("get_rms", &MDintegrate::get_rms,"")
    .def("get_summary", &MDintegrate::get_summary,"")
    .def_readwrite("target", &MDintegrate::target,"")
    .def_readwrite("maxsteps", &MDintegrate::maxsteps,"")
    ;

  py::class_<Equilibrate3d> aeq3d(mmdynamics, "Equilibrate3d", mdi);
  aeq3d
    .def(py::init<double, int, double, int, double, double, double, double, double, double>())
    .def("one_step", &Equilibrate3d::one_step,"")
    .def_readwrite("dt", &Equilibrate3d::dt,"")
    .def_readwrite("alpha", &Equilibrate3d::alpha,"")
    //
    .def("get_dx", &Equilibrate3d::get_dx,"")

    ;
}
