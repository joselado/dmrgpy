static auto conjugate_mps=[]() {
  // now get the MPS
  auto psi1 = read_wf("wf1.mps") ; // get the WF
  auto psi2 = conjMPS(psi1) ; // get the WF
  writeToFile("wf2.mps",psi2);
};



