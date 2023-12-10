

//using eigen and pybind11 pass a numpy array 

#include <iostream>
#include <Eigen/Dense>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

int pass_array(py::EigenDRef<Eigen::MatrixXd> m)
{
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;
}


int main()
{
  // do nothing
}


PYBIND11_MODULE(teigen, m) {

  m.def("pass_array", &pass_array, "");

}

