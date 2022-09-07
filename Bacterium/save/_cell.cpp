
//////////////////////////////////////////////////////////////////////////////////////
// start cell

void Cell::add_pilus(int i, double pl)
{
  int isshrink = (r2() > 0.5) ? 1 : 0;
  double teta = pili_distrib();
  Vector3d p_axis = theta_axis(teta);

  // need to deallocate this memory, todo
  std::shared_ptr<Pili> pilus (new Pili(i, pl, pl, p_axis, isshrink));
  pili.push_back(std::move(pilus));
  // always unbound
  boundidx.erase(pilus->idx);
  freeidx.insert(pilus->idx);
}


Vector3d Cell::real_pos()
{
  return hexpt;
}

Vector3d Cell::get_trail()
{
  // rodl from tail to head
  return hexpt - this->rodl * this->axis;
}

Vector3d Cell::real_centre()
{
  return hexpt - 0.5 * rodl * this->axis;
}

Vector3d Cell::get_vector()
{
  return this->rodl * this->axis;
}
Vector3d Cell::get_axis()
{
  return this->axis;
}



double Cell::pili_distrib()

  double teta = theta(axis);
  // split function into two step process to avoid recreating this object
  std::normal_distribution<double> ndistribution(teta, pilivar);
  double outt = ndistribution(mtgen);
  return outt;
}

fltq Cell::pili_stretch()
{
  Vector3d fl = Vector3d();
  for (auto& pilus : pili)
  {
    
    fl += pilus->force_l() * pilus->axis;
  }
  // vector from the centre to the head
  // this axis is tail -> head
  Vector3d r = rodl/2. * axis;
  
  //calculate torque
  Vector3d paxis = axis.perp2d();
  Vector3d orthforce = dot(paxis, fl) * paxis;
  double torque = r.xycross(orthforce);

  fltq ft;
  ft.force = fl;
  ft.torque = torque;
  // turn off
  //ft.torque = 0.;

  return ft;
}



void Cell::update_pos(Vector3d dxy)
{
  // planar
  hexpt = hexpt + dxy;
}

void Cell::update_axis(double dangle)
{
  // dangle is the small angle to rotate the cell body based on current torque
  // this axis is tail -> head

  Vector3d head_axis = 0.5 * rodl * this->axis;

  Vector3d centre = hexpt - head_axis;
  head_axis = head_axis.rotate(dangle, e_z);
  // update head position
  this->hexpt = centre + head_axis;
  this->axis = (2./rodl) * head_axis;
  // apply periodic boundaries through Surface code
  for (auto& pilus : pili)
  {
    
    // we don't mess with the axis/axisEq of bound pili
    // they should get updated, pilus.update()
    //pilus->axisEq = pilus->axisEq.rotate( dangle, e_z);
    if (!pilus->isbound) {
      // always set free pili to have their equilibrium angles?
      pilus->axis = pilus->axisEq;
    }
  }
}

void Cell::update(Vector3d dxy, double dangle) {
  // order of pos and angle updates shouldn't matter
  this->update_pos(dxy);
  this->update_axis(dangle);
  this->update_pili();
}

void Cell::respawn(Pili& pilus, double tme)
{
  pilus.pl = Pili::d_free;
  pilus.leq = Pili::d_free;
  double teta = pili_distrib();
  pilus.axis = theta_axis(teta);
  pilus.axisEq = pilus.axis;
  this->detach(pilus, tme);

  pilus.isshrink = (r2() > 0.5) ? 1 : 0;
  pilus.free_time = tme;
  pilus.lifetime = tme; // track pilus time between respawn events
}


//

/*Cell Cell::copy() */
//{
  //Cell ncell = Cell(this->idx, this->hexpt, this->axis);
  //for (int i; i < npili; i++) {
    //ncell.pili[i] = *this->pili[i]->copy(); // presumable not correct ?
  //}
  //return ncell;
//}



void Cell::attach(Pili& pilus, double tme) {
  pilus.attach(tme, this->hexpt);
  boundidx.insert(pilus.idx);
  freeidx.erase(pilus.idx);
  // when the pilus attaches it should contract?
}


/// output utilities

tuple<double, double> Cell::energies()
{
  double total_en_s = 0.;
  double total_en_b = 0.;
  for (int idx : boundidx)
  {
    tuple<double, double> en_xx = pili[idx]->energies();
    total_en_s += get<0>(en_xx);
    total_en_b += get<1>(en_xx);
  }
  return std::make_tuple(total_en_s, total_en_b);
}

double Cell::pili_deviation(int pidx) 
{
  return angle(axis, pili[pidx]->axis);
}

std::string Cell::__str__() {
  return "Object Cell -- __str__ not implemented";
}


