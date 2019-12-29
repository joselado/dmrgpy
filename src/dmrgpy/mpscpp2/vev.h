
static auto vev=[]() {
  auto sites = get_sites(); // 
  auto psi = read_wf(get_str("starting_file_gs")) ; // read GS
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





static auto excited_vev=[]() {
  auto A = get_mpo_operator("vev_multioperator.in"); // get the operator
  auto sites = get_sites(); // Get the different sites
  auto H = get_hamiltonian(sites) ; // get the Hamiltonian
  auto sweeps = get_sweeps(); // get sweeps
  auto nexcited = get_int_value("nexcited") ; // get excited states
  auto wfs = get_excited(H,sites,sweeps,nexcited); // get excited states
  ofstream filevev;
  filevev.open("VEV.OUT"); // open file
  for(int i=0;i<nexcited;i++) {
              auto c = overlapC(wfs.at(i),A,wfs.at(i)); // compute overlap
	      filevev << std::setprecision(16) << real(c) << "  ";
	      filevev << std::setprecision(16) << imag(c) << "  ";
	      filevev << endl;
  };
  filevev.close(); // close file
};
