
#include "periodic.hpp"

#include <exception>
#include <iostream>
  
#include <boost/algorithm/clamp.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>


using std::get;
using namespace boost::algorithm;

bool debug = true;
bool debug_iter = false;

nlopt::opt opt(nlopt::LD_TNEWTON, 2);
nlopt::opt yopt(nlopt::LN_NELDERMEAD, 2);
//nlopt::opt opt(nlopt::LN_COBYLA , 2);
//nlopt::opt opt(nlopt::LN_SBPLX, 2);

double ROOT_TOL = 1e-5;
double MIN_CONTACT_SEP = 1e-3;
double REL_TOL = 1e-3;
double MIN_SEP_M = 0.5;


// setup the target precision of minimisation routines
// quarter precision
int newton_bits = std::numeric_limits<double>::digits/2;
int bits = std::numeric_limits<double>::digits/4;

// tmp -- hardcoded to avoid reverse communication between modules 
const double cellR = 0.5;
const double cellL = 2.;


//int MAXITER = 100;
//int LV_MAXITER = 40;

// tmp
int LV_MAXITER = 100;
int MAXITER = 1000;

typedef std::tuple<double, double> fdf;
typedef boost::function<fdf(double)> minable;

double GRAD_MIN = 1e-8;
double f_capleq::operator()(const std::vector<double> &mx, std::vector<double> &grad)
{
  //clcl_data *clcl = (clcl_data *) data;
  ++data.itercount;
  SinePlane* sp;
  double c1,c2,l1,l2;
  c1 = data.c1;
  c2 = data.c2;
  l1 = data.l1;
  l2 = data.l2;
  sp = data.sp;
  double m,x,m1,m2,fx,dfx;
  m = mx[0];
  x = mx[1];
  m1 = m*l1 + c1;
  m2 = m*l2 + c2;
  fx = sp->form(x);
  dfx = sp->dform(x);
  if (!grad.empty()) {
    grad[0] = 2 * l1 * (m1 - x) + 2 * l2 * (m2 - sp->form(x) );
    grad[1] =  2*(x - m1) + 2*sp->dform(x) * ( sp->form(x) - m2 );
  }
  if (debug_iter) {
    cout << "evaluating f_capleq at " << m << " " << x << endl;
    cout << "gradient at " << grad[0] << " " << grad[1] << endl;
  }
  return (x - m1)*(x - m1) + (fx - m2)*(fx - m2);
}

SinePlane::SinePlane(Frame frame, double A, double B, int nsign) 
    : A(A), B(B), nsign(nsign), invB(1./B), xperiod(2 *M_PI *B)
{
  this->frame = frame;
  std::tuple<int,int,int> ixyz = cartesian_basis_vectors(frame);
  std::tie(iperiod, iconst, inorm) = ixyz;
  if (iperiod == NOTCART || iconst == NOTCART || inorm == NOTCART)
  {
    throw std::runtime_error("SinePlane should "
      "be constructed with a frame made of cartesian basis vectors.");
  }
  // construct functional struct
  // messing with scopes here
  // relative tolerance on optimisation parameters
  opt.set_xtol_rel(REL_TOL);
}


Vector3d SinePlane::sp_form(double t) 
{
  // origin is always going to be 0 in the constant direction
  Vector3d out = get_origin(); 
  out[iperiod] = t;
  out[inorm] = this->form(t);
  return out;
}

Vector3d SinePlane::normal_form(double t)
{
  Vector3d out{};
  out[iperiod] = -cos(t/B)/B;
  out[inorm] = 1;
  return out.unit();
}

Vector3d SinePlane::normal(Vector3d v)
{
  return normal_form(v[iperiod]);
}

bool SinePlane::contains(Vector3d pt) {
  if (nsign == 1) 
    return (form(pt[iperiod]) < pt[inorm]);
  else
    return (form(pt[iperiod]) > pt[inorm]);
}


boost::format brent_error_text(
    "Brent's method failed to converge after %i iterations");
boost::format warning_error_text(
    "Warning: %s");
boost::format newton_error_text(
    "NR method failed to converge after %i iterations");


bool check_iter_brent(bool satisfied, int maxiter) {
  if (!satisfied) { 
    throw std::runtime_error(
        boost::str(brent_error_text % maxiter)
        );
  }
  return satisfied;
}

bool check_iter_newton(bool satisfied, int maxiter) {
  if (!satisfied) { 
    throw std::runtime_error(
        boost::str(newton_error_text % maxiter)
        );
  }
  return satisfied;
}

// NOT USED 
bool check_iter(boost::uintmax_t it, int maxiter, boost::format error_text,
    bool warn=false) 
{
  bool failed = (it >= maxiter);
  if (failed) { 
    if (warn) {
      cout << "Warning: ";
      cout << boost::str(error_text % maxiter) << endl;
    }
    else {
      throw std::runtime_error(
          boost::str(error_text % maxiter)
          );
    }
  }
  return failed;
}

double SinePlane::_closest_t(double xc, double yc) {
  // construct distance^2 function
  if (debug_iter) {
    cout << "attempt brents method for pt " << "( " << xc << " , " << yc << " )" << endl;
  }
  boost::function<double(double)> fd2 = [this,xc,yc](double x) { 
    double fx = this->form(x); 
    double d2 = (xc-x)*(xc-x) + (yc-fx)*(yc-fx);
    if (debug_iter) {
      cout << "d2 = " << d2 << endl;
    }
    return d2;
  };

  double guess = xc;
  boost::uintmax_t it = MAXITER;
  std::pair<double,double> result = boost::math::tools::brent_find_minima(fd2, 
      guess-cellR, guess+cellR, bits, it);
  bool satisfied = (it < MAXITER);
  if (debug && !satisfied) {
    cout << "Failed brents method for pt " << "( " << xc << " , " << yc << " )" << endl;
  }
  check_iter_brent(satisfied, MAXITER);
  //cout << "number of iterations " << it << endl;
  return result.first;
}

double SinePlane::closest_t(Vector3d pt) 
{
  return this->_closest_t(pt[iperiod], pt[inorm]);
}


Vector3d SinePlane::distance_vector(Vector3d pt) 
{
  double t = closest_t(pt);
  Vector3d ret = pt;
  ret[iperiod] = t;
  ret[inorm] = form(t);
  return pt - ret;
}

void debug_newton(Lvec lv, boost::uintmax_t maxiter)
{
  cout << "NR Failed to converge for Lvec after " << maxiter << " iterations" << endl;
  cout << lv << endl;
}


bool newton_satisfied(double tvalue, double target) {
  return (abs(tvalue-target) < ROOT_TOL);
}

boost::optional<Vector3d> SinePlane::intersects(Lvec& lv)
{
  double c1, c2, l1, l2;
  bool satisfied = false;
  Vector3d src = lv.get_origin();
  Vector3d tpt = lv.get_endpt();
  c1 = src[iperiod];
  c2 = src[inorm];
  l1 = lv.line[iperiod];
  l2 = lv.line[inorm];
  // line segment equation root finding (using square of distance)
  minable lseqf = [this,l1,l2,c1,c2](double m) {
    double m1,m2,fx,dfx,eqev,sqeqev,deqev;
    m1 = m*l1 + c1;
    m2 = m*l2 + c2;
    fx = this->form(m1);
    dfx = this->dform(m1);
    eqev = m2 - fx;
    sqeqev = eqev*eqev;
    deqev = 2*eqev*(l2 - dfx * l1);
    if (debug_iter) {
      cout << "m eqev deqev "<< endl;
      cout << m << " " << eqev << " " << deqev << endl;
    }
    return std::make_tuple(sqeqev, deqev); 
  };
  
  // forward search for root, if result != 0; no root found
  double teqev, tderiv, fvalue, rvalue;
  std::tie(teqev, tderiv) = lseqf(0.);
  if (-tderiv < 0) {
    // skip ahead
  }
  else {
    boost::uintmax_t it = LV_MAXITER;
    double nrforward  = boost::math::tools::newton_raphson_iterate(lseqf, 0., 
        0., 1., newton_bits, it);
    satisfied = ( (it < LV_MAXITER) || newton_satisfied(nrforward, 1.) );
    if (!satisfied) debug_newton(lv, LV_MAXITER);
    
    check_iter_newton(satisfied, LV_MAXITER);
    fvalue = (nrforward*l2 + c2) - form(nrforward*l1 + c1);
    if (debug_iter) cout << "fvalue " << fvalue << endl;
    if (abs(fvalue) < ROOT_TOL) {
      return lv.get_pt(nrforward);
    } 
  }

  // check from the other side
  std::tie(teqev, tderiv) = lseqf(1.);
  if (-tderiv > 0) {
    //skip ahead
  }
  else {
    boost::uintmax_t it = LV_MAXITER;
    double nrbackward = boost::math::tools::newton_raphson_iterate(lseqf, 1., 
        0., 1., newton_bits, it);

    satisfied = ( (it < LV_MAXITER) || newton_satisfied(nrbackward, 0.) );
    if (!satisfied) debug_newton(lv, LV_MAXITER);
    check_iter_newton(satisfied, LV_MAXITER);

    rvalue = (nrbackward*l2 + c2) - form(nrbackward*l1 + c1);
    if (debug_iter) cout << "rvalue " << rvalue << endl;
    if (abs(rvalue) < ROOT_TOL) {
      return lv.get_pt(nrbackward);
    } 
  }
  
  // so long as we set up the system right we should have found the root by now
  return boost::none;
}

// tmp
double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}



vector<overlap> SinePlane::overlap_vector(Capsule& body) 
{ 
  // debugging
  _now_contact.clear();
  _contact_m.clear();
  // return
  vector<overlap> ret;
  overlap over; 

  double R2 = body.R*body.R;
  Lvec lv = body.get_lvec();
  Vector3d src = lv.get_origin();
  Vector3d tpt = lv.get_endpt();

  if ( (1 - abs(body.get_axis()[iconst])) < 1e3) {
    // the body essentially lines up in the constant direction and experiences a gradient less surface
    double hd2, td2, hd, td;
    Vector3d hdv, tdv;
    hdv = distance_vector(src);
    tdv = distance_vector(tpt);
    hd2 = hdv.len2();
    td2 = tdv.len2();
    if (hd2 < R2) {
      hd = sqrt(hd2);
      ret.push_back(std::make_tuple(body.get_m(0), hd, (1/hd) * hdv));
    }
    if (td2 < R2) {
      td = sqrt(td2); 
      ret.push_back(std::make_tuple(body.get_m(1), td, (1/td) * tdv));
    }
    return ret;
  }

  double c1,c2,l1,l2,hd,td,hd2,td2;
  c1 = src[iperiod];
  c2 = src[inorm];
  l1 = lv.line[iperiod];
  l2 = lv.line[inorm];
  // determine function from functional
  clcl_data data;
  data.sp = this;
  data.c1 = c1;
  data.c2 = c2;
  data.l1 = l1;
  data.l2 = l2;
  capleq.data = data;

  std::vector<double> lb(2);
  std::vector<double> ub(2);
  // 0 = m, 1 = x
  lb[0] = 0; ub[0] = 1;
  // set the x range to +/- R of the bodylv 
  lb[1] = c1; ub[1] = c1+l1;
  if (lb[1] > ub[1]) {
    double tx;
    tx = ub[1];
    ub[1] = lb[1];
    lb[1] = tx;
  }
  lb[1] -= body.R; ub[1] += body.R;
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  if (debug_iter) {
    cout << "x bounds " << lb[1] << " " << ub[1] << endl;
  }

  // set bounds // constraints
  // so ... need to pass a static function to set_min_objective?

  opt.set_min_objective(f_capleq::wrap, &capleq); 
  opt.set_maxeval(MAXITER);

  std::vector<double> mx_forward(2);
  mx_forward[0] = 0;
  mx_forward[1] = c1;

  std::vector<double> mx_reverse(2);
  mx_reverse[0] = 1;
  mx_reverse[1] = c1 + l1;
  
  double fmind2, rmind2;

  // do forward work
  try{
    capleq.data.itercount = 0;
    nlopt::result result = opt.optimize(mx_forward, fmind2);
  }
  catch(std::exception &e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
    throw e;
  }
  if (debug_iter) cout << "forward search niter " << capleq.data.itercount << endl;
  if (capleq.data.itercount == MAXITER) {
    throw std::runtime_error("EXIT: maxiterations in overlap_vector");
  }


  bool hsolved, tsolved;
  hsolved = tsolved = false;
  Vector3d dvf, dvr, lpt;
  double mf, mr;

  if (fmind2 < R2) {
    // this minimum represents and overlap
    mf = body.get_m(mx_forward[0]);
    lpt = lv.get_pt(mx_forward[0]);
    lpt[iconst] = 0;
    dvf = lpt - sp_form(mx_forward[1]);
    over = std::make_tuple(mf, sqrt(fmind2), (1/sqrt(fmind2))*dvf);
    ret.push_back(over);
    // debugging
    _now_contact.push_back(sp_form(mx_forward[1]));
    _contact_m.push_back(mf);
    hsolved = true;
  }

  // adjust bounds
  lb[0] += std::min(lb[0] + MIN_SEP_M, 1.);
  opt.set_lower_bounds(lb);

  // do reverse work
  try{
      capleq.data.itercount = 0;
      nlopt::result result = opt.optimize(mx_reverse, rmind2);
  }
  catch(std::exception &e) {
      std::cout << "nlopt failed: " << e.what() << std::endl;
      throw e;
  }
  if (debug_iter) cout << "reverse search niter " << capleq.data.itercount << endl;
  if (capleq.data.itercount == MAXITER) {
    throw std::runtime_error("EXIT: maxiterations in overlap_vector"); 
  }

  // if hsolved, check if we found the same pt

  if (rmind2 < R2) {
    // check overlaps are different
    mr = body.get_m(mx_reverse[0]);
    if (hsolved && (abs(mr - mf) < MIN_SEP_M)) 
    { 
      if (debug_iter) cout << "found the same overlap, exit early" << endl;
      return ret; 
    } else {
      // continue to add the next contact
      lpt = lv.get_pt(mx_reverse[0]);
      lpt[iconst] = 0;
      dvr = lpt - sp_form(mx_reverse[1]);
      over = std::make_tuple(mr, sqrt(rmind2), (1/sqrt(rmind2))*dvr);
      ret.push_back(over);
      // debugging
      _now_contact.push_back(sp_form(mx_reverse[1]));
      _contact_m.push_back(mr);
    }
  }

  return ret;
}



/*vector<overlap> SinePlane::overlap_vector(Capsule& body) */
//{
  //[> debugging <]
  //if (debug_iter) {
    //cout << "Attempt distance calculations for body. " << endl;
    //cout << body << endl;
  //}
  //_now_contact.clear();
  //_contact_m.clear();
  //[> <]

  //// temporary implementation without multidimensional minimisation library
  //vector<overlap> ret;
  //bool satisfied = false;

  //// we always solve in both directiions, then check if minima are the same
  //double R2 = body.R*body.R;
  //double c1,c2,l1,l2,hd,td,hd2,td2;
  //Lvec lv = body.get_lvec();
  //Vector3d src = lv.get_origin();
  //Vector3d tpt = lv.get_endpt();

  //// should first check distance vectors of head and tail? 
  
  //// NOT QUITE. This has implications for breaking contacts

  //bool hsolved, tsolved, hearly_solved, tearly_solved;
  //hsolved = tsolved = false;
  //hearly_solved = tearly_solved = false;

  
  //Vector3d hdv, tdv;
  //hdv = distance_vector(src);
  //tdv = distance_vector(tpt);
  //hd2 = hdv.len2();

  //if (hd2 < R2) {
    //hearly_solved = true;
    //[>hd = sqrt(hd2);<]
    ////ret.push_back(std::make_tuple(body.get_m(0), hd, (1/hd)*hdv));
    ////if (debug_iter) {
      ////cout << "found front overlap for m = -1" << ret.back() << endl;
    //[>}<]
  //}
  //td2 = tdv.len2();
  //if (td2 < R2) {
    //tearly_solved = true;
    //[>td = sqrt(td2);<]
    ////ret.push_back(std::make_tuple(body.get_m(1), td, (1/td)*tdv));
    ////if (debug_iter) {
      ////cout << "found rear overlap for m = 1" << ret.back() << endl;
    //[>}<]
  //}

  //c1 = src[iperiod];
  //c2 = src[inorm];
  //l1 = lv.line[iperiod];
  //l2 = lv.line[inorm];
  
  //// SLOW function
  //boost::function<double(double)> capleq = [this,l1,l2,c1,c2](double m) {
    //// x of closest pt on surface to pt on line
    //if (debug_iter) cout << "Minimise top dimension " << endl;
    //double lx,ly;
    //lx = m*l1 + c1; ly = m*l2 + c2;
    //double x = this->_closest_t(lx, ly);
    //double fx = this->form(x);
    //double capd2 = (x-lx)*(x-lx) + (fx-ly)*(fx-ly);
    //if (debug_iter) cout << "m = " << m << " capd2 = " << capd2 << endl;
    //return capd2;
  //};
  //// reverse transform
  //boost::function<double(double)> rcapleq = [capleq](double m_dash) {
    //double m = 1 - m_dash;
    //return capleq(m);
  //};

  //// look in both halves
  //double nrforward, nrbackward, df2, dr2;

  //double xf, xr;
  //Vector3d vf, vr, dvf, dvr;

  //double forward_limit = 0.5;
  //double reverse_limit = 0.5;

  //// SLOW SOLVE HEAD
  //boost::uintmax_t it;
  //if (!hsolved) {
    //it = MAXITER;
    //if (debug_iter) cout << "Attempt solve head-half intersection" << endl;
    //std::tie(nrforward, df2) = boost::math::tools::brent_find_minima(capleq,  
        //0., forward_limit, bits, it);
    //nrforward = clamp(nrforward, 0., forward_limit);
    //if (debug_iter) {
      //cout << "nrforward " << nrforward << endl;
      ////cout << body.get_m(nrforward) << endl;
    //}
    //satisfied = (it < LV_MAXITER);
    //if (debug && !satisfied ) {
      //cout << "Failed forward brent for body Lvec" << endl;
      //cout << lv << endl;
    //}
    //check_iter_brent(satisfied, LV_MAXITER);
    //if (df2 < R2) { hsolved = true; }
  //}

  ////REMINDER: overlap = {body_m (-ve points to head), distance to bodylv, unit vector}

  //overlap over; 
  //double mt, mh;
  //if (hsolved || hearly_solved) {
    //if (df2 < hd2) {
      //// another unnecessary call to iterative function
      //// ...
      //xf = this->closest_t(lv.get_pt(nrforward));
      //vf = sp_form(xf);
      //dvf = lv.get_pt(nrforward) - vf;
      //dvf[iconst] = 0; // !! PROJECTION !!

      //mh = body.get_m(nrforward);
      //over = std::make_tuple(mh, sqrt(df2), dvf.unit());
    //}
    //else {
      //hd = sqrt(hd2);
      //mh = body.get_m(0);
      //over = std::make_tuple(mh, hd, (1/hd)*hdv);
    //}

    //if (debug_iter) {
      //cout << "found front overlap " << over << endl;
    //}
    //ret.push_back(over);
    //_now_contact.push_back(body.get_pt(mh) - get<1>(over) * get<2>(over));
    //_contact_m.push_back(mh);

    //// setup for the reverse search
    ////reverse_limit = 1- body.get_t(m);
  //}


  //// SLOW  SOLVE TAIL
  //if (!tsolved) {
    //it = MAXITER;
    //double nrbackward_dash;
    //if (debug_iter) cout << "Attempt solve tail-half intersection" << endl;
    ////off
    ////std::tie(nrbackward_dash, dr2) = boost::math::tools::brent_find_minima(rcapleq,  
        ////0., reverse_limit, bits, it);
    ////nrbackward = 1 - nrbackward_dash; // reverse reverse transform
    ////on
    //std::tie(nrbackward, dr2) = boost::math::tools::brent_find_minima(capleq,  
        //reverse_limit, 1., bits, it);
    //nrbackward = clamp(nrbackward, reverse_limit, 1.);
    ////nrbackward = 1 - nrbackward_dash; // reverse reverse transform
    //if (debug_iter) {
      //cout << "nrbackward " << nrbackward << endl;
      ////cout << body.get_m(nrbackward) << endl;
    //}
    //satisfied = (it < LV_MAXITER);
    //if (debug && !satisfied ) {
      //cout << "Failed reverse brent for body Lvec" << endl;
      //cout << lv << endl;
    //}
    //check_iter_brent(satisfied, LV_MAXITER);
    //if (dr2 < R2) { tsolved = true;  }
    
  //}


  //if (tsolved || tearly_solved) {
    //if (dr2 < td2) {
      //xr = this->closest_t(lv.get_pt(nrbackward));
      //vr = sp_form(xr);
      //dvr = lv.get_pt(nrbackward) - vr;
      //dvr[iconst] = 0; // !! PROJECTION !!
      //mt = body.get_m(nrbackward);
      //over = std::make_tuple(mt, sqrt(dr2), dvr.unit());
    //}
    //else {
      //td = sqrt(td2);
      //mt = body.get_m(1);
      //over = std::make_tuple(mt, td, (1/td)*tdv);
    //}
    //if (debug_iter) {
      //cout << "found rear overlap " << over << endl;
    //}

    //// TODO check exit_early condition
    //bool exit_early =  (
        //(hsolved || hearly_solved) 
        //&& (abs(mh-mt) < MIN_CONTACT_SEP)
        //);
    //if (exit_early) { return ret; }

    //ret.push_back(over);
    //_now_contact.push_back(body.get_pt(mt) - get<1>(over) * get<2>(over));
    //_contact_m.push_back(mt);
  //}

  //return ret;
//}


std::string SinePlane::report_touching()
{
  if (_now_contact.empty()) {
    return "No touches.";
  }
  std::string report = "--Touching List Start: m v_contact\n";
  for (int i = 0; i < _now_contact.size(); ++i) {
    report += std::to_string(_contact_m[i]) + " ";
    report += _now_contact[i].__str__();
    report += "\n";
  }
  report += "--End";
  return report;
}


