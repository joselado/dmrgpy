
static auto get_hubbard =[](auto ampo) {
    ifstream jfile; // file to read
    jfile.open("hubbard.in"); // file with the hoppings
    int nt;
    jfile >> nt; // read the number of couplings
    int j1; // declare index
    int j2; // declare index
    auto U=0.0;
    // this function only considers spin independent hoppings
    for (int i=0;i<nt;++i) {
//      cout << i << endl ;
      jfile >> j1 >> j2 >> U; // get the data
      if ((site_type(j1)==1) and (site_type(j2)==1))  {
          ampo += U,"Nup",j1+1,"Nup",j2+1;
          ampo += U,"Nup",j1+1,"Ndn",j2+1;
          ampo += U,"Ndn",j1+1,"Ndn",j2+1;
          ampo += U,"Ndn",j1+1,"Nup",j2+1; }
    // spinless fermions
      if ((site_type(j1)==0) and (site_type(j2)==0))  
          ampo += U,"N",j1+1,"N",j2+1;
    } ;
    jfile.close() ;
    return ampo ;  // return the Hamiltonian with exchange added
}
;
