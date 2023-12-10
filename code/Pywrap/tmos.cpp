
//wrap the whole project
//
//
//
//
#include <pybind11/pybind11.h>

namespace py = pybind11;

void def_base(py::module &);
void def_surface(py::module &);
void def_pili(py::module &);
void def_mdynamics(py::module &);
void def_vtkwriter(py::module &);

PYBIND11_MODULE(tmos, m) {

  // implementations of these functions are provided by the linker so long as the 
  // relevant files are part of the compilation
  def_base(m);
  def_surface(m);
  def_pili(m);
  def_mdynamics(m);
  def_vtkwriter(m);

}


