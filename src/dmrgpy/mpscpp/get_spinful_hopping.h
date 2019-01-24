
static auto get_spinful_hopping=[](auto ampo) {
    ifstream jfile; // file to read
    jfile.open("spinful_hoppings.in"); // file with the hoppings
    int nt;
    jfile >> nt; // read the number of couplings
    int j1; // declare index
    int j2; // declare index
    int indi,indj;
    auto tr=0.0;
    auto ti=0.0;
    auto namei = "aaa";
    auto namej = "aaa";
    // this function only considers spin independent hoppings
    for (int i=0;i<nt;++i) {
//      cout << i << endl ;
      jfile >> j1 >> j2 >> tr >> ti; // get the data
      indi = j1/2+1; // index
      indj = j2/2+1; // index
      if ((site_type(j1/2)==1) and (site_type(j2/2)==1))  {
	  if (j1%2==0) namei = "Cdagup";
	  if (j1%2==1) namei = "Cdagdn";
	  if (j2%2==0) namej = "Cup";
	  if (j2%2==1) namej = "Cdn";
          ampo += tr,namei,indi,namej,indj;
          ampo += ti*1i,namei,indi,namej,indj;
          //ampo += ti*1i,"Cdagdn",j1+1,"Cdn",j2+1; 
      }
    } ;
    jfile.close() ;
    return ampo ;  // return the Hamiltonian with exchange added
}
;