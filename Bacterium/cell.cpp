
#include "cell.hpp"

#include <boost/format.hpp>


// Start ACell
double ACell::polarisation = 1.;

double ACell::pilivar = 4.;
double ACell::maxpl = -1;
double ACell::pili_min_length = 0.02; // static variable set in constructor
double ACell::k_spawn = 0.58; // 35 pili min^-1
int ACell::spawn_extension_state = 1;  
// Can be an optional member variable for ACell
bool ACell::running_start = false;
bool ACell::pilus_replacement_on = false;


CellEvent* ACell::create_event() {
  return new CellEvent( std::move(this->clone()) );
}


void ACell::common_init() {
  // create n pili
  for (int pili_idx = 0; pili_idx < npili; ++pili_idx)
  {
    std::shared_ptr<Pili> pilus = this->spawn_pilus();
    this->add_pilus( pilus );
  }
}


void ACell::add_pilus(std::shared_ptr<Pili> pilus) {
  this->pili.push_back(pilus); // store the pointer
  pili_counter += 1;
  this->init_anchors();
}


void ACell::prep_dissolve_pilus(std::shared_ptr<Pili> pilus) {
  if (pilus->is_fully_retracted) {
    this->dissolve_pilus(pilus);
  }
  else {
    pilus->dissolve_on_retraction = true;
  }

}

void ACell::dissolve_pilus(std::shared_ptr<Pili> pilus) {
  pili.erase(std::remove_if(pili.begin(),pili.end(),
        [pilus](std::shared_ptr<Pili> p){return pilus->idx == p->idx;}));
  this->init_anchors();
}

vector<vector<double> > ACell::get_rates()
{
  vector<vector<double>> rates;
  for (auto& pilus : pili)
  {
    rates.push_back(pilus->get_rates());
  }
  return rates;
}

vector<double> ACell::get_cell_rates()
{
  vector<double> cell_rates{k_spawn};
  return cell_rates;
}

CellEvent* ACell::elongate(Pili& pilus, double tme) {
  double leq = pilus.elongate();

  bool pistaut = pilus.istaut();
  if (pilus.isbound && (pistaut != pilus.lastistaut))
  {
    std::string trigger = (pistaut) ? "nowtaut" : "nowslack";
    pilus.lastistaut = pistaut;
    CellEvent* e = new CellEvent(tme, pilus.idx, string("elongate"), 
        std::move(this->clone()), trigger);
    return e;
  }

  if (Pili::force_max_length && leq > Pili::max_length)
  {
    pilus.detach(tme);
    CellEvent* e = new CellEvent(tme, pilus.idx, string("elongate"), 
        std::move(this->clone()),
        string("dissolve"));
    return e;
  }
  return NULL;
}

CellEvent* ACell::shrink(Pili& pilus, double tme) {
  if (Pili::enforce_stalling && abs(pilus.force_l()) >= Pili::f_stall)
  {
    assert(false);
    // do nothing and return early
    return NULL;
  }
  // else shrink the pilus
  double leq = pilus.shrink();

  // check whether bound pili become taut this step
  bool pistaut = pilus.istaut();
  if (pilus.isbound && (pistaut != pilus.lastistaut))
  {
    std::string trigger = (pistaut) ? "nowtaut" : "nowslack";
    CellEvent* e = new CellEvent(tme, pilus.idx, string("shrink"), 
        std::move(this->clone()), trigger);
    pilus.lastistaut = pistaut;
    return e;
  }

  // this is the actual way we dissolve pili (at the end of kmc step)
  // if set !pilus.isbound condition here we get a huge number of bound pili
  // EDIT: add condition to prevent dissolving bound pili
  // if (leq < (Pili::inside_length) && pilus.cycles > 0)
  if (leq < (Pili::inside_length) && pilus.cycles > 0 && !pilus.isbound)
  {
    pilus.detach(tme);
    CellEvent* e = new CellEvent(tme, pilus.idx, string("shrink"), this->clone(),
        string("dissolve"));
    return e;
  }
  return NULL;
}

CellEvent* ACell::switch_se(Pili& pilus, double tme) {
  pilus.switch_se();
  if (pilus.isbound) {
    string sw{"switching"};
    CellEvent* e = new CellEvent(tme, pilus.idx, sw, this->clone(), sw);
    return e;
  }
  else {
    return NULL;
  }
}

double ACell::pl_avg() 
{
  double plsum = 0;
  for (auto& pilus : pili) {
    plsum += pilus->leq;
  }
  return pili.empty() ? 0. : plsum/pili.size(); 
}

double ACell::l_total() 
{
  double lsum = 0;
  for (auto& pilus : pili) {
    lsum += pilus->leq;
  }
  return lsum;
}


double ACell::energy_pili()
{
  double en_total = 0.;
  for (auto pilus : this->pili)
  {
    double ens = pilus->energy();
    en_total += ens;
  }
  return en_total;
}

double ACell::energy() { 
  cout << "Cell Aell::energy" << endl;
  return energy_pili(); 
}

std::shared_ptr<Pili> ACell::get_pilus(int pidx)
{
  std::shared_ptr<Pili> pilus = nullptr;
  for (std::shared_ptr<Pili> p : this->pili)
  {
    if (p->idx == pidx)
    {
      pilus = p;
      break;
    }
  }
  if (pilus == nullptr) {
    throw std::invalid_argument( "No pilus with this idx" );
  }
  return pilus;
}

int ACell::get_pilus_vidx(std::shared_ptr<Pili> p)
{
  int vidx = -1;
  for (int i=0; i<this->pili.size(); i++)
  {
    if (pili[i]->idx == p->idx) {
      vidx = i;
      break;
    }
  }
  if (vidx == -1) {
    throw std::invalid_argument( "Pilus not in cell.pili" );
  }
  return vidx;
}

double ACell::pbrf() { 
  double sum = 0;
  for (auto pilus: pili) {
    sum += pilus->get_bound_retract_fraction();
  }
  return sum;
}

void swap(ACell& lhs, ACell& rhs) {
  using std::swap;
  swap(lhs.idx, rhs.idx);
  swap(lhs.npili, rhs.npili);
  swap(lhs.pili, rhs.pili);
}


std::ostream& operator<<(std::ostream& os, ACell& cell)
{
  os << cell.__str__();
  return os;
}



/////////////////////////////////////////////////////////////////////////////////


string CellEvent::__str__() {
  string form = "CellEvent: pilus %d of cell %d executed process %s at t=%12.6f\n";
  if (!trigger.empty())
    form += (boost::format("triggered %s\n") % trigger ).str();
  form += "Dump cell state:\n";
  ACell& cell = get_data();
  form += cell.__str__();
  return ( boost::format(form)
      % pidx % cell.idx % process % get_time() 
      ).str(); 
}


