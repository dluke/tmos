
#include <memory>

#include "frame.hpp"
#include "vector3d.hpp"
#include "shape.hpp"
#include "plane.hpp"
#include "sphere.hpp"
#include "pili.hpp"
#include "cell.hpp"
#include "cell3d.hpp"
#include "kmc.hpp"
#include "wlc.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
namespace py = pybind11;


void def_pili(py::module &m) {
  auto mpili = m.def_submodule("pili");

  py::class_<Pili, std::shared_ptr<Pili> > pili(mpili, "Pili");
  pili.def(py::init<int,double,double,Vector3d&,Vector3d&,int,int>())
    .def("clone", &Pili::clone, "")
    .def("istaut", &Pili::istaut, "")
    .def("tetaEq", &Pili::tetaEq, "")
    .def("update", &Pili::update, "")
    .def("force_l", &Pili::force_l, "")
    //
    .def("tau_eovers", &Pili::tau_eovers, "")
    .def("calc_k_dt", &Pili::calc_k_dt, "")
    .def("energy", &Pili::energy, "")
    .def("get_anchor_angle", &Pili::get_anchor_angle, "")
    .def("get_bound_retract_fraction", &Pili::get_bound_retract_fraction, "")
    .def("get_bound_retract_rate", &Pili::get_bound_retract_rate, "")
    .def("get_rates", &Pili::get_rates, "")
    .def("elongate", &Pili::elongate, "")
    .def("shrink", &Pili::shrink, "")
    .def("shrink_by", &Pili::shrink_by, "")
    .def("switch_se", &Pili::switch_se, "")
    // .def("attach", (int (Pili::*) (double, Vector3d)) &Pili::attach, "")
    //.def("attach", (int (PiliRod3d::*) (double, Capsule, Shape*)) &PiliRod3d::attach, "")
    .def("attach", &Pili::attach, "")
    .def("shorten_to_detach", &Pili::shorten_to_detach, "")
    .def("switch_se", &Pili::switch_se, "")
    // .def("attach", (int (Pili::*) (double, Vector3d)) &Pili::attach, "")
    //.def("attach", (int (PiliRod3d::*) (double, Capsule, Shape*)) &PiliRod3d::attach, "")
    .def("attach", &Pili::attach, "")
    .def("_attach_upkeep", &Pili::_attach_upkeep, "")
    .def("shorten_to_detach", &Pili::shorten_to_detach, "")
    .def("detach", &Pili::detach, "")
    .def_readwrite("idx", &Pili::idx)
    .def_readwrite("axisEq", &Pili::axisEq)
    .def_readwrite("anchor", &Pili::anchor)
    .def_readwrite("lab_axisEq", &Pili::lab_axisEq)
    .def_readwrite("lab_anchor", &Pili::lab_anchor)
    .def_readwrite("attachment", &Pili::attachment)
    .def_readwrite("pl", &Pili::pl)
    .def_readwrite("leq", &Pili::leq)
    .def_readwrite("isbound", &Pili::isbound)
    .def_readwrite("lastistaut", &Pili::lastistaut)
    .def_readwrite("last_shortening", &Pili::last_shortening)
    .def_readwrite("ext_motor", &Pili::ext_motor)
    .def_readwrite("ret_motor", &Pili::ret_motor)
    .def_readwrite("cycles", &Pili::cycles)
    .def_readwrite("new_cycle", &Pili::new_cycle)
    .def_readwrite("is_fully_retracted", &Pili::is_fully_retracted)
    .def_readwrite("bound_time", &Pili::bound_time)
    .def_readwrite("free_time", &Pili::free_time)
    .def_readwrite("lifetime", &Pili::lifetime)

    // model parameters ( class variables )
    .def_readwrite_static("ks", &Pili::ks)
    .def_readwrite_static("force_min_length", &Pili::force_min_length)
    .def_readwrite_static("max_length", &Pili::max_length)
    .def_readwrite_static("inside_length", &Pili::inside_length)
    .def_readwrite_static("d_free", &Pili::d_free)
    .def_readwrite_static("d_bound", &Pili::d_bound)
    .def_readwrite_static("k_adh", &Pili::k_adh)
    .def_readwrite_static("adh_accept", &Pili::adh_accept)
    .def_readwrite_static("f_stall", &Pili::f_stall)
    .def_readwrite_static("kf_sh", &Pili::kf_sh)
    .def_readwrite_static("kf_ex", &Pili::kf_ex)
    .def_readwrite_static("kb_sh", &Pili::kb_sh)
    .def_readwrite_static("kb_ex", &Pili::kb_ex)
    .def_readwrite_static("force_threshold", &Pili::force_threshold)
    .def_readwrite_static("anchor_angle_threshold", &Pili::anchor_angle_threshold)
    .def_readwrite_static("anchor_angle_smoothing_fraction", &Pili::anchor_angle_smoothing_fraction)
    //.def_readwrite_static("rse_factor", &Pili::rse_factor)
    .def_readwrite_static("f_release", &Pili::f_release)
    .def_readwrite_static("k_ext_off", &Pili::k_ext_off)
    .def_readwrite_static("k_ret_off", &Pili::k_ret_off)
    .def_readwrite_static("k_ext_on", &Pili::k_ext_on)
    .def_readwrite_static("k_ret_on", &Pili::k_ret_on)
    .def_readwrite_static("k_dissolve", &Pili::k_dissolve)
    .def_readwrite_static("k_resample", &Pili::k_resample)
    .def_readwrite_static("detach_grace_time", &Pili::detach_grace_time)
    .def_readwrite_static("detach_grace_length", &Pili::detach_grace_length)
    .def_readwrite_static("dwell_time", &Pili::dwell_time)
    .def_readwrite_static("tau_e", &Pili::tau_e)
    .def_readwrite_static("free_tau_e", &Pili::free_tau_e)
    .def_readwrite_static("allow_bound_ext_motor", &Pili::allow_bound_ext_motor)
    .def_readwrite_static("allow_bound_extension", &Pili::allow_bound_extension)
    .def_readwrite_static("force_bound_retraction", &Pili::force_bound_retraction)
    .def_readwrite_static("force_detachment", &Pili::force_detachment)
    .def_readwrite_static("simplify_free_pili_etor", &Pili::simplify_free_pili_etor)
    .def_readwrite_static("rcondition", &Pili::rcondition)
    .def_readwrite_static("enforce_stalling", &Pili::enforce_stalling)
    .def_readwrite_static("free_tau_eovers", &Pili::free_tau_eovers)
    .def_readwrite_static("low_force_proportion", &Pili::low_force_proportion)
    .def_readwrite_static("high_force_proportion", &Pili::high_force_proportion)
    .def_readwrite_static("force_max_length", &Pili::force_max_length)

    .def("get_typestr", &Pili::get_typestr, "")
    .def("__str__", &Pili::__str__, "")
    ;

  py::class_<Pili3d, std::shared_ptr<Pili3d> > pili3d(mpili, "Pili3d", pili); // use name 'pili'
  pili3d.def(py::init<int,double,double,Vector3d&,Vector3d&,int,int>())
    .def("get_last_inter_len", &Pili3d::get_last_inter_len)
    ;

  py::class_<PiliRod3d, std::shared_ptr<PiliRod3d> >(mpili, "PiliRod3d", pili3d)
    .def(py::init<int,double,double,Vector3d&,Vector3d&,int,int>())
    ;
  
  py::class_<PiliWLC3d, std::shared_ptr<PiliWLC3d> >(mpili, "PiliWLC3d", pili3d)
    .def(py::init<int,double,double,Vector3d&,Vector3d&,int,int,
        std::shared_ptr<WLCgeneratorBase> >())
    .def(py::init<int,double,double,Vector3d&,Vector3d&,int,int>())
    .def_readwrite("n", &PiliWLC3d::n, "")
    .def_readwrite("lasta", &PiliWLC3d::lasta, "")
    .def("transformed_chain_instance", &PiliWLC3d::transformed_chain_instance, "")
    ;

  //
  py::class_<WLCgeneratorBase, std::shared_ptr<WLCgeneratorBase> > wlcbase(mpili, "WLCgeneratorBase");
  wlcbase
    .def(py::init<double,double,int>())
    ;
  py::class_<KPgenerator, std::shared_ptr<KPgenerator> > kpgen(mpili, "KPgenerator", wlcbase);
  kpgen 
    .def(py::init<double,double,int,double>())
    ;

  py::class_<WLCgenerator, std::shared_ptr<WLCgenerator> > wlcg(mpili, "WLCgenerator", wlcbase);
  wlcg
    .def(py::init<double,double,int>())
    .def("append", &WLCgenerator::append, 
        py::return_value_policy::reference,
        py::keep_alive<0, 1>(),
        py::keep_alive<1, 2>()
        )
    .def("setnth", &WLCgenerator::setnth, 
        py::return_value_policy::reference,
        py::keep_alive<0, 1>(), 
        py::keep_alive<1, 3>() 
        )
    .def("getnth", &WLCgenerator::getnth, py::return_value_policy::reference)

    .def("sample", &WLCgenerator::sample, "")
    ;


  py::class_<ACell> acell(mpili, "ACell"); // expose name 'acell'
  acell
    .def_readwrite_static("pilivar", &ACell::pilivar)
    .def_readwrite_static("k_spawn", &ACell::k_spawn)
    .def_readwrite_static("spawn_extension_state", &ACell::spawn_extension_state)
    .def_readwrite_static("polarisation", &ACell::polarisation)
    .def_readwrite_static("maxpl", &ACell::maxpl)
    .def_readwrite_static("pili_min_length", &ACell::pili_min_length)
    .def_readwrite_static("running_start", &ACell::running_start)
    .def_readwrite_static("pilus_replacement_on", &ACell::pilus_replacement_on)

    .def_readwrite("idx", &ACell::idx)
    .def_readwrite("npili", &ACell::npili)
    .def_readwrite("pili", &ACell::pili) // C++ vector converts to python list
    .def_readwrite("pili_counter", &ACell::pili_counter)

    .def("create_event", &ACell::create_event, "")
    .def("common_init", &ACell::common_init, "")
    .def("set_npili", &ACell::set_npili, "")
    .def("shrink ", &ACell::shrink, "")
    .def("detach", &ACell::detach, "")
    .def("spawn_pilus", &ACell::spawn_pilus, "")
    .def("add_pilus", &ACell::add_pilus, "")
    .def("dissolve_pilus", &ACell::dissolve_pilus, "")

    // can define pure virtual methods here and have them work for subclasses
    .def("real_pos", &ACell::real_pos, "")
    .def("get_headpt", &ACell::get_headpt, "")
    .def("get_trail", &ACell::get_trail, "")
    .def("real_centre", &ACell::real_centre, "")
    .def("get_axis", &ACell::get_axis, "") 
    .def("get_pilus", &ACell::get_pilus, "") 
    .def("get_pilus_vidx", &ACell::get_pilus_vidx, "") 
    
    .def("nbound", &ACell::nbound, "")
    .def("ntaut", &ACell::ntaut, "")
    .def("num_pili", &ACell::num_pili, "")
    .def("pbrf", &ACell::pbrf, "")
    .def("pl_avg", &ACell::pl_avg, "")
    .def("l_total", &ACell::l_total, "")
    .def("energy", &ACell::energy, "")
    .def("energy_pili", &ACell::energy_pili, "")
    .def("clone", &ACell::clone, "")
    ;

  // Cell3d subclass
  py::class_<Cell3d> cell3d(mpili, "Cell3d", acell);
  cell3d
    .def(py::init<Capsule&>())
    .def(py::init<int,int,Capsule&,Plane*>(),
        py::keep_alive<1,5>()
        )
    //.def(py::init<int,int,Capsule&,Step&>()) // tmp
    .def("__str__", &Cell3d::__str__, "")
    .def("pili_body_distrib", &Cell3d::pili_body_distrib, "")
    .def("pili_body_position", &Cell3d::pili_body_position, "")

    // detailed construction
    //.def("add_pilus", (void (Cell3d::*)(int, double, Vector3d, Vector3d))&Cell3d::add_pilus, "")
    .def("add_pilus", &Cell3d::add_pilus, "")

    .def_readwrite_static("eps", &Cell3d::eps, "")
    .def_readwrite_static("repulsive_only", &Cell3d::repulsive_only, "")
    .def_readwrite_static("attract_width", &Cell3d::attract_width, "")
    .def_readwrite_static("cost_anchor_intersection", &Cell3d::cost_anchor_intersection, "")
    .def_readwrite_static("distrib_type", &Cell3d::distrib_type, "")

    .def_readwrite("body", &Cell3d::body, "")
    .def_readwrite("anchors", &Cell3d::anchors, "")
    .def_readwrite("axisEqs", &Cell3d::axisEqs, "")
    .def_readwrite("surface", &Cell3d::surface, "")

    .def("get_lab_anchor", &Cell3d::get_lab_anchor, "")
    .def("get_lab_axisEq", &Cell3d::get_lab_axisEq, "")
    .def("update_anchors", &Cell3d::update_anchors, "") // for output only
    .def("reset_surface", &Cell3d::reset_surface, "")
    .def("get_num_contacts", &Cell3d::get_num_contacts, "")
    .def("energy_surface", &Cell3d::energy_surface, "")
    .def("set_rotation", &Cell3d::set_rotation, "")
    // get the gradient as a list?
    .def("pili_grad", &Cell3d::pili_grad, "")
    .def("surface_grad", &Cell3d::surface_grad, "")
    .def("grad", &Cell3d::grad, "")
    .def("set_state", &Cell3d::set_state, "")
    .def("get_state", &Cell3d::get_state, "")
    .def("state_energy", &Cell3d::state_energy, "")
    .def("state_gradient", &Cell3d::state_gradient, "")
    .def("report_touching", &Cell3d::report_touching, "")
    ;
  // free constructors
  //mpili.def("make_cell3d_on_step", &make_cell3d_on_step, "");

  py::class_<CellWLC3d> cellwlc3d(mpili, "CellWLC3d", cell3d);
  cellwlc3d

    .def(py::init<int,int,Capsule&,
        Plane*,
        std::shared_ptr<WLCgeneratorBase> >(),
        py::keep_alive<1,5>(),
        py::keep_alive<1,6>()
        )
    .def("get_wlc", &CellWLC3d::get_wlc, "")
    .def("__str__", &CellWLC3d::__str__, "")
    ;

  py::class_<CellEvent, std::shared_ptr<CellEvent> > cellevent(mpili, "CellEvent");
  cellevent
    .def("get_time", &CellEvent::get_time, "")
    .def("get_data", &CellEvent::get_data, "", 
        py::return_value_policy::reference)
    .def_readwrite("pidx", &CellEvent::pidx, "")
    .def_readwrite("process", &CellEvent::process, "")
    .def_readwrite("trigger", &CellEvent::trigger, "")
    .def("__str__", &CellEvent::__str__, "") 
    ;

  py::class_<Kmc> kmc(mpili, "Kmc");
  kmc
    .def(py::init<double>())
    // .def("get_event", &Kmc::get_event, 
        // py::return_value_policy::reference)
    .def("kmcstep", &Kmc::kmcstep, "")
    .def("attachstep", &Kmc::attachstep, "")
    .def("postmd", &Kmc::postmd, "")
    // C++ vector converts to python list
    // default uses return_value_policy::reference_internal
    .def_readwrite("events", &Kmc::events,
        py::return_value_policy::reference)
    ;
  
}

