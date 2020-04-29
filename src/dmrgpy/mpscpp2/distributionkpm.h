// compute the moments of a distribution using the KPM

static auto get_moments_distribution=[]()
{
  auto sites = get_sites(); // Get the different sites
  auto H = get_hamiltonian(sites) ; // get the Hamiltonian
  auto delta = get_float_value("kpm_delta");
  auto n = get_int_value("kpm_num_polynomials") ; // number of polynomials
  ofstream myfile;
  myfile.open("KPM_NUM_POLYNOMIALS.OUT");
  myfile << std::setprecision(8) << n << endl;
  myfile.close(); // close file
  auto psi = get_gs() ; // get the ground state
  // now read the operator to use
  auto m = get_mpo_operator("kpm_distribution_multioperator.in");
  ///////////////////////////////////
  ////// once the operators have been read, continue
  ///////////////////////////////////
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto kpmcutoff = get_float_value("kpm_cutoff") ; // bond dimension for KPM
  moments_vi_vj(m,psi,psi,n) ; //compute the KPM moments
  return 0 ;
} ;

