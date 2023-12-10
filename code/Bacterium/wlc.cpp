

#include "wlc.hpp"

#include <iostream>
#include <exception>

using std::cout;
using std::endl;


std::unique_ptr<Eigen::VectorXd> WLCgeneratorBase::sample_phi(void) {
  Eigen::VectorXd* phi = new Eigen::VectorXd(nmax);
  // how to do this better ?
  auto& p = *phi;
  for (int i = 0; i < nmax; i++) {
    p[i] = r2() * 2 * M_PI;
  }
  return std::unique_ptr<Eigen::VectorXd>(phi);
}

// ----------------------------------------------------------------------------------

double KPgenerator::one_theta(void) {
  double x = r2();
  double Ei = -log( exp(Lp/a) - Z*beta*x)/beta;
  double angle = acos(-a*Ei/ka);
  return angle;
}

std::unique_ptr<Eigen::VectorXd> KPgenerator::sample_theta(void) {
  // TODO: we dont need to sample a full chain here
  Eigen::VectorXd* theta = new Eigen::VectorXd(nmax);
  for (int i = 0; i < nmax; i++) {
    double angle = this->one_theta();
    (*theta)[i] = angle;
  }
  return std::unique_ptr<Eigen::VectorXd>(theta);
}

// ----------------------------------------------------------------------------------

void WLCgenerator::append(py::EigenDRef<Eigen::MatrixXd> cdata) 
{
  //cout << cdata << endl;
  this->vcdata.push_back(cdata);
}

void WLCgenerator::setnth(int n, py::EigenDRef<Eigen::MatrixXd> cdata) 
{
  this->vcdata[n] = cdata;
}

py::EigenDRef<Eigen::MatrixXd> WLCgenerator::getnth(int n)
{
  return this->vcdata[n];
}

std::unique_ptr<Chain> WLCgenerator::sample(int n, Vector3d axis) 
{
  int nbend = n - 1;
  if (nbend == 0) {
    return std::make_unique<Chain>(axis, get_a());
  }
  else if (nbend > 0) {
    // generate 
    py::EigenDRef<Eigen::MatrixXd> cdata = vcdata[nbend-1];
    std::uniform_int_distribution<int> uni{0,cdata.rows()-1};
    // slice a row
    int r = uni(mtgen);
    Eigen::VectorXd chain_theta = cdata.row(r);
    
    return std::make_unique<Chain>(axis, get_a(), chain_theta, *sample_phi());
  }
}

std::unique_ptr<Eigen::VectorXd> WLCgenerator::sample_theta(void) {
  py::EigenDRef<Eigen::MatrixXd> cdata = vcdata.back();
  std::uniform_int_distribution<int> uni{0,cdata.rows()-1};
  int r = uni(mtgen);
  Eigen::VectorXd* row = new Eigen::VectorXd(nmax);
  // copy
  row->operator=( cdata.row(r) );
  // wrap
  return std::unique_ptr<Eigen::VectorXd>(row);
}
