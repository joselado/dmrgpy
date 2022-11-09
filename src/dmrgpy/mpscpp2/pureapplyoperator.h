auto get_mpo(std::string name="MPO.mpo") {
  auto sites = get_sites();
  auto Op = Iden(sites); // get the identity
  readFromFile(name,Op);
  return Op ;
}




static auto pureapplyoperator=[]() {
  // now get the operator
  auto A = get_mpo(get_str("pureapplyoperator_operator"));
  auto psi0 = read_wf(get_str("pureapplyoperator_wf0")) ; // get the WF
  int maxm = get_int_value("maxm") ; // bond dimension
  auto cutoff = get_float_value("cutoff") ; // cutoff
  auto args = Args("Cutoff",cutoff,"Maxm",maxm);
  auto psi1 = exactApplyMPO(psi0,A,args) ;
  writeToFile(get_str("pureapplyoperator_wf1"),psi1);
};



static auto gen_pureoperator=[]() {
    auto A = get_mpo_operator(get_str("gen_pureoperator_operator_in"));
    writeToFile(get_str("gen_pureoperator_operator_out"),A);
};





static auto overlap_aMb_static=[]() {
  // now get the MPS
  auto psi1 = read_wf("overlap_aMb_wf1.mps") ; // get the WF
  auto psi2 = read_wf("overlap_aMb_wf2.mps") ; // get the WF
  auto A = get_mpo("overlap_aMb_pureoperator.mpo"); // get the MPO
  auto c = overlapC(psi1,A,psi2) ; // compute overlap
  ofstream ofile; // declare
  ofile.open("OVERLAP_aMb.OUT");  // open file
  ofile << std::setprecision(20) << real(c) << "  " << imag(c) << endl ;
  ofile.close() ; // close file
};


