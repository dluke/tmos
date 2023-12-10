
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "sample.hpp"
#include "vector3d.hpp"
#include "frame.hpp"

void def_base(py::module &m) {

  auto mbase = m.def_submodule("base");

  mbase.def("init_generator", &init_generator, "");

  py::class_<Vector3d>(mbase, "Vector3d")
    .def(py::init<>())
    .def(py::init<double, double, double>())
    //.def(py::self = py::self)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def("__neg__", [](Vector3d &v) {
        return Vector3d(-v.x, -v.y, -v.z);
    }, py::is_operator())
    .def(py::self * double())
    .def(double() * py::self)
    .def(py::self == py::self)
    .def(py::self += py::self)
    .def(py::self -= py::self)
    .def("cross", &Vector3d::cross, "cross product")
    .def("len", &Vector3d::len, "")
    .def("len2", &Vector3d::len2, "")
    .def("scale", &Vector3d::scale, "")
    .def("norm", &Vector3d::norm, "")
    .def("unit", &Vector3d::unit, "")
    .def("rotate", &Vector3d::rotate, "")
    .def("xyrotate", &Vector3d::xyrotate, "")
    .def("perp2e", &Vector3d::perp2d, "")
    .def("dot", &dot, "dot product")
    .def("angle", &angle, "")
    .def("theta", &theta, "")
    .def("theta_axis", &theta_axis, "")
    .def_readwrite("x", &Vector3d::x)
    .def_readwrite("y", &Vector3d::y)
    .def_readwrite("z", &Vector3d::z)
    .def("__str__", &Vector3d::__str__, "")
    ;
    //.def("list", [](const Vector3d &s) { return py::list(s.x, s.y, s.z); } )

  py::class_<Matrix3d>(mbase, "Matrix3d")
    .def(py::init<>())
    .def("__str__", &Matrix3d::__str__, "")
    ;

  py::class_<Frame> frame(mbase, "Frame");
  frame
    .def(py::init<>())
    .def(py::init<Vector3d&,Vector3d&>())
    .def(py::init<Vector3d&,Vector3d&,Vector3d&>())
    .def(py::init<Vector3d&,Vector3d&,Vector3d&,Vector3d&>())
    .def("get_rmatrix", &Frame::get_rmatrix, "")
    .def("get_tmatrix", &Frame::get_tmatrix, "")
    .def("get_origin", &Frame::get_origin, "")
    .def("to_lab_rt", &Frame::to_lab_rt, "")
    .def("to_lab", &Frame::to_lab, "")
    .def("orthogonalise", &Frame::orthogonalise, "")
    .def("normalise", &Frame::normalise, "")
    .def("orthogonal_error", &Frame::orthogonal_error, "")
    .def("unit_error", &Frame::unit_error, "")
    .def("__str__", &Frame::__str__, "")
    .def_readwrite("origin", &Frame::origin)
    .def_readwrite("e1", &Frame::e1)
    .def_readwrite("e2", &Frame::e2)
    .def_readwrite("e3", &Frame::e3)
    ;

}

