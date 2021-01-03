static auto applyoperator=[]() {
  auto sites = get_sites();
  // now get the operators
  auto A = get_mpo_operator(get_str("applyoperator_multioperator"));
  auto psi0 = read_wf(get_str("applyoperator_wf0")) ; // get the WF
  int maxm = get_int_value("maxm") ; // bond dimension
  auto cutoff = get_float_value("cutoff") ; // cutoff
  auto args = Args("Cutoff",cutoff,"Maxm",maxm);
  auto psi1 = exactApplyMPO(psi0,A,args) ;
  writeToFile(get_str("applyoperator_wf1"),psi1);
};



static auto get_summps=[]() {
  // now get the MPS
  auto psi1 = read_wf("summps_wf1.mps") ; // get the WF
  auto psi2 = read_wf("summps_wf2.mps") ; // get the WF
  auto psi3 = sum_mps(psi1,psi2) ;
  writeToFile("summps_wf3.mps",psi3);
};


static auto overlap_aMb=[]() {
  // now get the MPS
  auto psi1 = read_wf("overlap_aMb_wf1.mps") ; // get the WF
  auto psi2 = read_wf("overlap_aMb_wf2.mps") ; // get the WF
  auto A = get_mpo_operator("overlap_aMb_M.in"); // get the operator
  auto c = overlapC(psi1,A,psi2) ; // compute overlap
  ofstream ofile; // declare
  ofile.open("OVERLAP_aMb.OUT");  // open file
  ofile << std::setprecision(20) << real(c) << "  " << imag(c) << endl ;
  ofile.close() ; // close file
};


