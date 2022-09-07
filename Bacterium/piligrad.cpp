
#include "piligrad.hpp"

using std::cos;
using std::sin;


vector<double> eref = {e_z.x, e_z.y, e_z.z};

vector<Matrix3d> dR_phide_k(Vector3d ehat, double phi) {
  vector<Matrix3d> pdRphi(3);
  Matrix3d skewe = skew_matrix(ehat);
  for (int k = 0; k < 3; k++) {
    Matrix3d esk = skewk[k];
    pdRphi[k] = (1-cos(phi)) * (esk*skewe + skewe*esk) + sin(phi) * esk;
  }
  return pdRphi;
}

vector<Matrix3d> dR_phidpk(vector<Matrix3d> Rk, Vector3d ehat, double phi)
{
  vector<Matrix3d> dRpk(3);
  vector<Matrix3d> dRdek = dR_phide_k(ehat, phi);
  for (int k=0; k<3; k++) {
    Matrix3d cumm{};
    for (int j=0; j<3; j++) {
      cumm = cumm + (dRdek[j] + Rk[k] * eref[j]);
    }
    dRpk[k] = cumm;
  }
  return dRpk;
}

/*vector<double> dekdpk(vector<Matrix3d> Rk) // reference axis is e_z*/
//{
  //vector<double> ret(3);
  //for (int k = 0; k < 3; k++) {
    //Matrix3d tr = Rk[k];
    //ret[k] = tr.x3 + tr.y3 + tr.z3;
  //}
  //return ret;
/*}*/

vector<Matrix3d> dR_thetadp_k(vector<Matrix3d> Rk, double thetap, Vector3d e2)
{
  vector<Matrix3d> pdR(3);
  Matrix3d ske2 = skew_matrix(e2);
  for (int k=0; k<3; k++) {
    Matrix3d Rkey = skew_matrix(Rk[k]*e_y);
    pdR[k] = (1-cos(thetap))*(Rkey*ske2 + ske2*Rkey)
      + sin(thetap) * Rkey;
  }
  return pdR;
}

Matrix3d dR_phid_phi(double phi, Vector3d ehat)
{
  Matrix3d ske = skew_matrix(ehat);
  return sin(phi) *ske*ske + cos(phi)*ske;
}

