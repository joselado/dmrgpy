auto compute_overlap() {
  ofstream myfile;
  auto wf1 = read_wf("overlap_wf1.mps") ; // read first wavefunction 
  auto wf2 = read_wf("overlap_wf2.mps") ; // read second wavefunction 
  auto out = overlapC(wf1,wf2) ; // compute overlap
  myfile.open("OVERLAP.OUT");
  myfile << std::setprecision(8) << real(out) << endl;  // compute
  myfile << std::setprecision(8) << imag(out) << endl;  // compute
}

