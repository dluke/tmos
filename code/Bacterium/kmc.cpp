    
#include "kmc.hpp"

// process list
vector<string> plist = {"extension", "retraction", "ret_on", "ret_off", "ext_on", "ext_off",
  "resample", "detach", "dissolve", "spawn"};

// If process is detachment or extension/retraction of a bound pilus we need to goto MD step
bool md_condition(int procidx, int isbound) { 
  return (procidx == action::Detach) || (
    isbound == 1 && (procidx == action::Extension || procidx == action::Retraction));
}

void Kmc::postmd(ACell& cell, double tme) {
  attachstep(cell, tme);
  // updating lab frame pili anchors and axes is important for angle flexibility
  cell.update_anchors();
}

// Attchment is checked for at regular intervals, rather than as an MC process
// called from python main loop
void Kmc::attachstep(ACell& cell, double tme)
{
  
  // this could potentially generate multiple attachment events this step
  for (std::shared_ptr<Pili> pilus : cell.pili) {
    if (!pilus->isbound) {
      CellEvent* e = cell.attach(*pilus, tme);
      if (e != NULL) {
      // overwrite the process
        e->process = string("mdstep");
        this->set_this_event(e);
      }
    }
  }
}

// Kinetic MC step
kmctup Kmc::kmcstep(ACell& cell, double tme)
{
  vector<vector<double>> rates = cell.get_rates();
  vector<double> cell_rates = cell.get_cell_rates();
  int rsize = rates.empty() ? 0 : rates[0].size();
  int nrates = rates.size() * rsize + cell_rates.size();
  vector<double> cumsum(nrates);
  double csum = 0.;
  int i = 0;
  for (vector<double> v : rates)
  {
    for (int k = 0; k < v.size(); k++)
    {
      double rate = v[k];
      csum += rate;
      cumsum[i] = csum;
      i++;
    }
  }

  double pili_sum = csum;

  for (double rate : cell_rates) 
  {
    csum += rate;
    cumsum[i] = csum;
    i++;
  }

  double total_rate = csum; 
  double choice = r2() * csum;
  int procidx, pidx;
  string process_name;
  std::shared_ptr<Pili> pilus = nullptr;
  int unique_idx;
  if (choice < pili_sum)
  {
    // binary search
    vector<double>::iterator pit = std::lower_bound(
        cumsum.begin(), cumsum.end(), choice);
    int pindex = pit - cumsum.begin();
    pidx = pindex / rsize;
    procidx = (pindex % rsize); 
    pilus = cell.pili[pidx];
    unique_idx = pilus->idx;
    process_name = plist[procidx];
  }
  else 
  {
    procidx = action::Spawn;
    process_name = plist[action::Spawn];
    pidx = -1; // spawn event not yet associated with pilus
  }

  this->clear_events();

  CellEvent* e = NULL;
  switch (procidx) {
    case action::Extension:
      e = cell.elongate(*pilus, tme);
      break;
    case action::Retraction:
      e = cell.shrink(*pilus, tme);
      break;
    case action::Ret_on:
      pilus->ret_motor_on();
      break;
    case action::Ret_off:
      pilus->ret_motor_off();
      break;
    case action::Ext_on:
      // tracking extension cycles
      pilus->ext_motor_on();
      this->set_this_event(
          new CellEvent(tme, pilus->idx, process_name,
            std::move(cell.clone()), string("ext_on"))
        );
      break;
    case action::Ext_off:
      // tracking free pili maximum extension
      pilus->ext_motor_off();
      this->set_this_event(
          new CellEvent(tme, pilus->idx, process_name,
            std::move(cell.clone()), string("ext_off"))
        );
      break;
    case action::Resample:
      pilus->sample();
      this->set_this_event(
          new CellEvent(tme, pilus->idx, process_name,
            std::move(cell.clone()), string("resample"))
        );
      break;
    case action::Detach:
      e = cell.detach(*pilus, tme);
      break;
    case action::Dissolve:
      cell.prep_dissolve_pilus(pilus);
      break;
    case action::Spawn:
      pilus = cell.spawn_pilus();
      unique_idx = pilus->idx;
      cell.add_pilus(pilus);
      this->set_this_event(
          new CellEvent(tme, pilus->idx, process_name,
            std::move(cell.clone()), string("spawn"))
        );
  }
  
  // Save the cell state if appropriate
  if (e != NULL) {
    this->set_this_event(e);
  }

  // what if a pilus attaches and dissolves due to extension/retraction on the same step

  // if shrink/elongate triggered dissolve then call it
  // the shrink/elongate event carries the trigger="dissolve" event
  if (e != NULL && e->trigger == "dissolve") {
    cell.dissolve_pilus(pilus);
    pidx = -1; 
    if (cell.pilus_replacement_on) {
      auto pilus = cell.spawn_pilus();
      cell.add_pilus(pilus);
      this->set_this_event(
          new CellEvent(tme, pilus->idx, process_name, 
            std::move(cell.clone()), string("spawn")
            )
          );
    }
  }
  else if (!pilus->isbound && (procidx == action::Extension || procidx == action::Resample) )
  {
    // check attachment due to extension or resampling
    CellEvent* e = cell.attach(*pilus, tme);
    if (e != NULL) {
      // overwrite the process
      e->process = plist[procidx];
      this->set_this_event(e);
    }
  }
  
  // 
  int need_md = false;
  if (pilus != nullptr) {
    need_md = md_condition(procidx, pilus->isbound);
    if (need_md) {
      cell.last_pilus_ref = pilus;
    }
  }
  double trand = r2();
  double deltat = (1/total_rate) * log(1/trand);
  tme += deltat;

  kmctup kmc = std::make_tuple(tme, process_name, need_md, unique_idx);
  return kmc;
}

