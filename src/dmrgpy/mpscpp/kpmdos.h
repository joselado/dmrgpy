
static auto get_moments_dos=[](auto sites, auto H)
{
  auto n = get_int_value("nkpm");
  auto delta = get_float_value("kpm_delta");
  // define on the fly the number of polynomials
  n = int(round(bandwidth(sites,H)/delta))*get_int_value("kpm_n_scale");
  ofstream myfile;
  myfile.open("KPM_NUM_POLYNOMIALS.OUT");
  myfile << std::setprecision(8) << n << endl;
  myfile.close(); // close file
  auto m = scale_hamiltonian(sites,H) ; // scale this Hamiltonian
  auto psi = MPS(sites); // get a random initial wavefunction
  psi = psi*(1.0/sqrt(overlapC(psi,psi))) ; // renormalize
  moments_vi_vj(m,psi,psi,n) ; //compute the KPM moments
  return 0 ;
} ;



