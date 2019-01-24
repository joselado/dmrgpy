
static auto get_pairing = [](auto ampo) {
    ifstream jfile; // file to read
    jfile.open("pairing.in"); // file with the hoppings
    int nt;
    jfile >> nt; // read the number of couplings
    int j1; // declare index
    int j2; // declare index
    auto tr=0.0;
    auto ti=0.0;
    // this function only considers swave pairing
    for (int i=0;i<nt;++i) {
      jfile >> j1 >> j2 >> tr >> ti; // get the data
      if ((site_type(j1)==1) and (site_type(j2)==1))  {
          ampo += tr,"Cdn",j1+1,"Cup",j2+1;
          ampo += ti*1i,"Cdn",j1+1,"Cup",j2+1; 
          ampo += tr,"Cdagup",j1+1,"Cdagdn",j2+1;
          ampo += -ti*1i,"Cdagup",j1+1,"Cdagdn",j2+1; 
      }
    } ;
    jfile.close() ;
    return ampo ;  // return the Hamiltonian with exchange added
}
;
