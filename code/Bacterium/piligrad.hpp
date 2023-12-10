
#ifndef __PILIGRAD_HPP__
#define __PILIGRAD_HPP__ 

#include <cmath>
#include <vector>

#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "axisangle.hpp"

using std::vector;

// drop extern?
extern vector<double> eref;

vector<Matrix3d> dR_phide_k(Vector3d ehat, double phi);
vector<Matrix3d> dR_phidpk(vector<Matrix3d> Rk, Vector3d ehat, double phi);
//vector<double> dekdpk(vector<Matrix3d> Rk); 
vector<Matrix3d> dR_thetadp_k(vector<Matrix3d> Rk, double thetap, Vector3d e2);
Matrix3d dR_phid_phi(double phi, Vector3d ehat);

#endif
