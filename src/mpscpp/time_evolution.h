// perform a time evolution

int quench(auto sites) {
  auto name1 = get_str("tevol_operator_i");  // first operator
  auto name2 = get_str("tevol_operator_j");  // second operator
  auto i1 = get_int_value("tevol_site_i");  // second operator
  auto i2 = get_int_value("tevol_site_j");  // second operator
  // now get the operators
  auto A1 = get_spin_operator(sites,i1,name1);
  auto A2 = get_spin_operator(sites,i2,name2);
  auto H = get_hamiltonian(sites) ; // get the ampo for the Hamiltonian
  auto psi = get_gs(sites,H) ; // get the ground state
  int maxm = get_int_value("maxm") ; // bond dimension for KPM
  auto cutoff = get_float_value("cutoff") ; // bond dimension for KPM
  // apply the first operator
  // now do the time evolution
  auto nt = get_int_value("tevol_nt"); // number of time steps
  auto dt = get_float_value("tevol_dt"); // delta of time
  int it;
  auto args = Args("Cutoff=",cutoff,"Maxm=",maxm);
  // get the AutoMPO
  auto ampo = get_ampo(sites) ; // get the ampo for the Hamiltonian
  auto expH = toExpH<ITensor>(ampo,dt*Cplx_i); // get the exponential of the H
  ofstream fileevol; // file for the evolution
  fileevol.open("TIME_EVOLUTION.OUT"); // time evolution
  auto psi1 = exactApplyMPO(psi,A1,args) ;
  auto psi2 = exactApplyMPO(psi,A2,args) ;
//  normalize(psi1); // normalize
//  normalize(psi2); // normalize
  auto norm0 = sqrt(overlapC(psi1,psi1)) ;
  for (it=0;it<nt;it++) { // loop
	      psi1 = exactApplyMPO(expH,psi1,args); // evolve
              normalize(psi1); // normalize
	      psi1 *= norm0 ; // restore initial norm
	      auto z = overlapC(psi2,psi1) ; // overlap
//	      auto z = overlapC(psi,psi1) ; // overlap
	      // write in a file
	      fileevol << std::setprecision(8) << real(z) << "  "
                       << std::setprecision(8)<< imag(z) << endl;
  } ;
  fileevol.close(); // close file
}

