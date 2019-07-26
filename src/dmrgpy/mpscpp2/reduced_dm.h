static auto reduced_dm=[]() {
  auto sites = get_sites(); // 
  auto H = get_hamiltonian(sites) ; // get Hamiltonian
  auto psi = get_gs(sites,H) ; // get the ground state
  psi /= overlap(psi,psi); // normalize
  //Given an MPS called "psi",
  //and assuming j > i
  
  //'gauge' the MPS to site i
  //any 'position' between i and j, inclusive, would work here
  auto i = get_int_value("index_i_DM");
//  auto j = get_int_value("index_j_DM");
  psi.position(i); 
  
  //index linking i to i+1:
  auto ir = commonIndex(psi.A(i),psi.A(i+1));
  
  auto rho = psi.A(i)*dag(prime(psi.A(i),Site,ir));
  for(int k = i+1; k <= psi.N(); ++k)
      {
      rho *= psi.A(k);
      rho *= dag(prime(psi.A(k),Link));
      }
//  rho *= psi.A(j);
//  rho *= dag(prime(psi.A(j)));
//  for(int k = j+1; k <= psi.N(); ++k)
//      {
//      rho *= psi.A(k);
//      rho *= dag(prime(psi.A(k),Link));
//      }
//  ;
  // function to write the data
  ofstream filerdm; // file for the DM
  filerdm.open("DM.OUT"); // open file
  filerdm.close(); //close file
  auto doPrint = [](auto z) { 
          ofstream filerdm; // file for the DM
          filerdm.open("DM.OUT", ios::out | ios::app ); // open file
	  filerdm << std::setprecision(8) << real(z) << "  " ;
	  filerdm << std::setprecision(8) << imag(z) << endl ;
          filerdm.close(); //close file
  };
//  auto doPrint = [](auto x) { cout << x << endl; };
  // evaluate the tensor in the function
  cout << rho << endl;
  rho.visit(doPrint); // loop over tensor
};
