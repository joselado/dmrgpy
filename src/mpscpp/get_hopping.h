
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
          ampo += ti,"Cdagup",j1+1,"CupI",j2+1;
          ampo += tr,"Cdagdn",j1+1,"Cdn",j2+1; 
          ampo += ti,"Cdagdn",j1+1,"CdnI",j2+1; 
      }
          // complex conjugate
//          ampo += tr,"Cdagup",j2+1,"Cup",j1+1;
//          ampo += -ti,"Cdagup",j2+1,"CupI",j1+1;
//          ampo += tr,"Cdagdn",j2+1,"Cdn",j1+1; 
//          ampo += -ti,"Cdagdn",j2+1,"CdnI",j1+1; }
    } ;
    jfile.close() ;
    return ampo ;  // return the Hamiltonian with exchange added
}
