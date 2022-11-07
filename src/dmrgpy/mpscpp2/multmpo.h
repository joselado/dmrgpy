// function to multiply MPO

static auto multmpo_operator=[]() {
  auto A = get_mpo(get_str("multmpo_pureoperator_A"));
  auto B = get_mpo(get_str("multmpo_pureoperator_B"));
  auto C = mult_mpo(A,B);
  writeToFile(get_str("multmpo_pureoperator_C"),C) ;
};


static auto trace_mpo_operator=[]() {
  auto A = get_mpo(get_str("trace_pureoperator"));
  auto c = trace_mpo(A);
  ofstream ofile; // declare
  ofile.open("TRACE.OUT");  // open file
  ofile << std::setprecision(20) << real(c) << "  " << imag(c) << endl ;
  ofile.close() ; // close file
};



