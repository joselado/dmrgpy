
static auto get_vijkl =[](auto ampo) {
    ifstream jfile; // file to read
    jfile.open("vijkl.in"); // file with the interactions
    int nt;
    jfile >> nt; // read the number of couplings
    int i,j,k,l; // declare index
    auto U=0.0;
    // this function only considers spin independent hoppings
    for (int ii=0;ii<nt;++ii) {
//      cout << i << endl ;
      jfile >> i >> j >> k >> l >> U; // get the data
    // spinless fermions
      if ((site_type(i)==0) and (site_type(j)==0) and 
		      (site_type(k)==0) and (site_type(l)==0)
		      )  
          ampo += U,"Cdag",i+1,"C",j+1,"Cdag",k+1,"C",l+1;
    } ;
    jfile.close() ;
    return ampo ;  // return the Hamiltonian with exchange added
}
;
