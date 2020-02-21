static auto applyoperator=[]() {
  auto sites = get_sites();
  // now get the operators
  auto A = get_mpo_operator(get_str("applyoperator_multioperator.in"));
  auto psi0 = read_wf(get_str("applyoperator_wf0")) ; // get the WF
  int maxm = get_int_value("maxm") ; // bond dimension
  auto cutoff = get_float_value("cutoff") ; // cutoff
  auto args = Args("Cutoff=",cutoff,"Maxm=",maxm);
  auto psi1 = exactApplyMPO(psi0,A,args) ;
  writeToFile(get_str("applyoperator_wf1"),psi1);
};

