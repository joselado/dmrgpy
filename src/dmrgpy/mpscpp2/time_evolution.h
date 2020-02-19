// perform a time evolution

static auto evoloperator=[](auto H, auto dt) {
  auto sites = get_sites();
  auto ampo = AutoMPO(sites); // generate ampo
  ampo += 1.0,"Id", 1;
  std::complex<double> z(0.0,-dt);
  auto Iden = MPO(ampo); // identity
  auto out = sum_mpo(Iden,z*H) ;
  auto H2 = mult_mpo(H,H);
  auto H3 = mult_mpo(H,H2);
  out = sum_mpo(out,(0.5*z*z)*H2);
  out = sum_mpo(out,(z*z*z/6.0)*H2);
  return out;
}
;


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
  auto expH = MPO(ampo);
  if (get_bool("tevol_custom_exp")) 
    expH = evoloperator(MPO(ampo),dt) ; // create Hamiltonian
  else expH = toExpH<ITensor>(ampo,dt*Cplx_i); 
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
	      fileevol << std::setprecision(20) << real(z) << "  "
                       << std::setprecision(20)<< imag(z) << endl;
  } ;
  fileevol.close(); // close file
  writeToFile("psi_time_evolution.mps",psi1);
};



static auto evolution_AeiHtB=[]() {
  auto sites = get_sites();
  // now get the operators
  auto H = get_hamiltonian(sites) ; // get the ampo for the Hamiltonian
  int maxm = get_int_value("maxm") ; // bond dimension
  auto cutoff = get_float_value("cutoff") ; // cutoff
  auto nt = get_int_value("tevol_nt"); // number of time steps
  auto dt = get_float_value("tevol_dt"); // delta of time
  int it;
  auto args = Args("Cutoff=",cutoff,"Maxm=",maxm);
  // get the AutoMPO
  auto ampo = get_ampo(sites) ; // get the ampo for the Hamiltonian
  auto expH = MPO(ampo);
  if (get_bool("tevol_custom_exp")) 
    expH = evoloperator(MPO(ampo),dt) ; // create Hamiltonian
  else expH = toExpH<ITensor>(ampo,dt*Cplx_i); 
  ofstream fileevol; // file for the evolution
  fileevol.open("TIME_EVOLUTION.OUT"); // time evolution
  auto psi1 = read_wf(get_str("wfa_time_evolution.mps")) ;
  auto psi2 = read_wf(get_str("wfb_time_evolution.mps")) ;
  auto fittd = get_bool("tevol_fit_td") ; // use fitting method
  for (it=0;it<nt;it++) { // loop
	      if (fittd) fitApplyMPO(psi1,expH,psi1,args) ; // evolve
	      if (not fittd) psi1 = exactApplyMPO(expH,psi1,args); // evolve
	      auto z = overlapC(psi2,psi1) ; // overlap
//	      auto z = overlapC(psi,psi1) ; // overlap
	      // write in a file
	      fileevol << std::setprecision(20) << real(z) << "  "
                       << std::setprecision(20)<< imag(z) << endl;
  } ;
  fileevol.close(); // close file
  writeToFile("psi_time_evolution.mps",psi1);
};



static auto evolution_measure=[]() {
  auto sites = get_sites();
  // now get the operators
  auto A = get_mpo_operator("time_evolution_multioperator.in");
  auto H = get_hamiltonian(sites) ; // get the ampo for the Hamiltonian
  auto psi = read_wf("psi_evolve_and_measure.mps") ; 
  int maxm = get_int_value("maxm") ; // bond dimension
  auto cutoff = get_float_value("cutoff") ; // cutoff
  // apply the first operator
  // now do the time evolution
  auto nt = get_int_value("tevol_nt"); // number of time steps
  auto dt = get_float_value("tevol_dt"); // delta of time
  int it;
  auto args = Args("Cutoff=",cutoff,"Maxm=",maxm);
  // get the AutoMPO
  auto ampo = get_ampo(sites) ; // get the ampo for the Hamiltonian
  // shift by the ground state energy
  auto expH = MPO(ampo);
  if (get_bool("tevol_custom_exp")) 
    expH = evoloperator(MPO(ampo),dt) ; // create Hamiltonian
  else expH = toExpH<ITensor>(ampo,dt*Cplx_i); 
  ofstream fileevol; // file for the evolution
  fileevol.open("TIME_EVOLUTION.OUT"); // time evolution
  auto fittd = get_bool("tevol_fit_td") ; // use fitting method
  for (it=0;it<nt;it++) { // loop
	      if (fittd) fitApplyMPO(psi,expH,psi,args) ; // evolve
	      if (not fittd) psi = exactApplyMPO(expH,psi,args); // evolve
              auto z = overlapC(psi,A,psi) ;
	      // write in a file
	      fileevol << std::setprecision(20) << real(z) << "  "
                       << std::setprecision(20)<< imag(z) << endl;
  } ;
  fileevol.close(); // close file
  writeToFile("psi_evolve_and_measure.mps",psi);
};



