
auto get_hopping(auto ampo) {
    ifstream jfile; // file to read
    jfile.open("hoppings.in"); // file with the hoppings
    int nt;
    jfile >> nt; // read the number of couplings
    int j1; // declare index
    int j2; // declare index
    auto tr=0.0;
    auto ti=0.0;
    // this function only considers spin independent hoppings
    for (int i=0;i<nt;++i) {
//      cout << i << endl ;
      jfile >> j1 >> j2 >> tr >> ti; // get the data
      if ((site_type(j1)==1) and (site_type(j2)==1))  {
          ampo += tr,"Cdagup",j1+1,"Cup",j2+1;
          ampo += ti*1i,"Cdagup",j1+1,"Cup",j2+1;
          ampo += tr,"Cdagdn",j1+1,"Cdn",j2+1; 
          ampo += ti*1i,"Cdagdn",j1+1,"Cdn",j2+1; 
      }
    } ;
    jfile.close() ;
    return ampo ;  // return the Hamiltonian with exchange added
}



// function to test that the hopping is well implemented
auto test_hopping(auto H, auto sites) {
	cout << "BEGIN test hopping" << endl;
	for (int i=0;i<sites.N();i++) { // loop
	  for (int j=0;j<sites.N();j++) { // loop
            if ((site_type(i)==1) and (site_type(j)==1))  {
              auto state1 = InitState(sites,"Emp");
              state1.set(i+1,"Up");
              auto state2 = InitState(sites,"Emp");
              state2.set(j+1,"Up");
              auto psi1 = MPS(state1);
              auto psi2 = MPS(state2);
              auto c = overlap(psi2,H,psi1);
	      if (abs(c)>0.0001)   cout << i <<"  "<< j <<"  "<< c << endl ;
	    }
	  }
	}
	cout << "END test hopping" << endl;
}

