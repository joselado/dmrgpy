
static auto vev=[]() {
  auto sites = get_sites(); // 
  auto H = get_hamiltonian(sites) ; // get Hamiltonian
  auto psi = get_gs(sites,H) ; // get the ground state
  psi /= overlap(psi,psi); // normalize
  // now read the operator from tasks.in, this routine assumes
  // that there will be a multioperator
  auto A = get_mpo_operator("vev_multioperator.in"); // get the operator
  auto c = overlapC(psi,A,psi);
  ofstream ofile; // declare
  ofile.open("VEV.OUT");  // open file
  ofile << std::setprecision(20) << real(c) << "  " << imag(c) << endl ;
  ofile.close() ; // close file
  return 0; // dummy return
} ;
