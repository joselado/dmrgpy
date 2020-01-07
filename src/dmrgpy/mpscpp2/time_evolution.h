// perform a time evolution

static auto quench=[]() {
  auto sites = get_sites();
  // now get the operators
  auto A1 = get_mpo_operator("dc_multioperator_i.in");
  auto A2 = get_mpo_operator("dc_multioperator_j.in");
  auto H = get_hamiltonian(sites) ; // get the ampo for the Hamiltonian
  auto psi = get_gs() ; // get the ground state
  int maxm = get_int_value("maxm") ; // bond dimension
  auto cutoff = get_float_value("cutoff") ; // cutoff
  // apply the first operator
  // now do the time evolution
  auto nt = get_int_value("tevol_nt"); // number of time steps
  auto dt = get_float_value("tevol_dt"); // delta of time
  int it;
  auto args = Args("Cutoff=",cutoff,"Maxm=",maxm);
  // compute ground state energy
  auto EGS = overlap(psi,H,psi)/overlap(psi,psi); // ground state enrgy
  // get the AutoMPO
  auto ampo = get_ampo(sites) ; // get the ampo for the Hamiltonian
  // shift by the ground state energy
  ampo += -EGS,"Id", 1; // minus ground state energy
  auto expH = toExpH<ITensor>(ampo,dt*Cplx_i); // get the exponential of the H
  ofstream fileevol; // file for the evolution
  fileevol.open("TIME_EVOLUTION.OUT"); // time evolution
  auto psi1 = exactApplyMPO(psi,A1,args) ;
  auto psi2 = exactApplyMPO(psi,A2,args) ;
//  normalize(psi1); // normalize
//  normalize(psi2); // normalize
  auto norm0 = sqrt(overlapC(psi1,psi1)) ;
  auto fittd = get_bool("tevol_fit_td") ; // use fitting method
  for (it=0;it<nt;it++) { // loop
	      if (fittd) fitApplyMPO(psi1,expH,psi1,args) ; // evolve
	      if (not fittd) psi1 = exactApplyMPO(expH,psi1,args); // evolve
              normalize(psi1); // normalize
	      psi1 *= norm0 ; // restore initial norm
	      auto z = overlapC(psi2,psi1) ; // overlap
//	      auto z = overlapC(psi,psi1) ; // overlap
	      // write in a file
	      fileevol << std::setprecision(8) << real(z) << "  "
                       << std::setprecision(8)<< imag(z) << endl;
  } ;
  fileevol.close(); // close file
};

