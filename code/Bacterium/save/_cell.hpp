
class Cell: public ACell
{
public:

  // member variables
  Vector3d hexpt, axis;
  
  Cell(int idx, Vector3d hexpt, Vector3d axis) 
    : ACell(idx, def_npili), hexpt(hexpt), axis(axis) 
  {
    has_surface = false;
    this->common_init();
  }

  void add_pilus(int i, double pl) override;

  Vector3d real_pos() override;
  Vector3d get_trail() override;
  Vector3d real_centre() override;
  Vector3d get_vector() override;
  Vector3d get_axis() override;

  fltq pili_stretch();

  double pili_distrib() override;

  void update_pili() override;

  void update_pos(Vector3d dxy) override;

  void update_axis(double dangle) override;
  
  void update(Vector3d dxy, double dangle);

  void respawn(Pili& pilus, double tme);

  void attach(Pili& pilus, double tme);

  //Cell copy();
  
  // utilities
  unordered_set<int> get_bound_pili() { return boundidx; }

  tuple<double, double> energies();
  
  double teta() { return theta(axis); }

  // angle(axis, pili[pidx].axis)
  double pili_deviation(int pidx);
    
  virtual std::string __str__(); 
  

};


