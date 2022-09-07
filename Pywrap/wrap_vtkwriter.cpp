

#include "vtkwriter.hpp"
#include "vtksurface.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void def_vtkwriter(py::module &m) {
  auto mvtkwriter = m.def_submodule("vtkwriter");

  mvtkwriter.def("write_cell3d", &write_cell3d, "");
  mvtkwriter.def("write_pili3d", &write_pili3d, "");
  mvtkwriter.def("write_cell3d_from_body", &write_cell3d_from_body, "");
  mvtkwriter.def("write_pili_vectors", &write_pili_vectors, "");

  mvtkwriter.def("write_infsteps", &write_infsteps, "");
  mvtkwriter.def("write_hgrid", &write_hgrid, "", 
      py::arg("Hgrid"), py::arg("range") = 5, py::arg("res") = 10);
  mvtkwriter.def("write_segplane", &write_segplane, "");
  mvtkwriter.def("write_sineplane", &write_sineplane, "");
}


