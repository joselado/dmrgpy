
static auto get_pairing=[](auto ampo) {
    ifstream jfile; // file to read
    jfile.open("pairing.in"); // file with the hoppings
    int nt;
    jfile >> nt; // read the number of couplings
    int j1; // declare index
    int j2; // declare index
    auto tr=0.0; // real part
    auto ti=0.0; // imaginary part
    // this function only considers spin independent hoppings
    for (int i=0;i<nt;++i) {
//      cout << i << endl ;
      jfile >> j1 >> j2 >> tr >> ti; // get the data
      // spinless fermions
      if ((site_type(j1)==0) and (site_type(j2)==0))  {
          ampo += tr,"C",j1+1,"C",j2+1;
          ampo += ti*1i,"C",j1+1,"C",j2+1;
          ampo += tr,"Cdag",j2+1,"Cdag",j1+1;
          ampo += -ti*1i,"Cdag",j2+1,"Cdag",j1+1;
      }
    } ;
    jfile.close() ;
    return ampo ;  // return the Hamiltonian with exchange added
}
;


// this function is not being used, it was only written for debugging
static auto get_pairing_MPO=[](auto sites) {
    auto ampo = AutoMPO(sites) ; // get the ampo
    auto H = MPO(ampo); // generate Hamiltonian
    ifstream jfile; // file to read
    jfile.open("pairing.in"); // file with the hoppings
    int nt;
    jfile >> nt; // read the number of couplings
    int j1; // declare index
    int j2; // declare index
    auto tr=0.0; // real part
    auto ti=0.0; // imaginary part
    // this function only considers spin independent hoppings
    int maxm = get_int_value("maxm") ; // bond dimension for KPM
    float cutoff = get_float_value("cutoff") ; // cutoff for KPM
    for (int i=0;i<nt;++i) {
//      cout << i << endl ;
      jfile >> j1 >> j2 >> tr >> ti; // get the data
      // spinless fermions
      if ((site_type(j1)==0) and (site_type(j2)==0))  {
	      auto g = tr + ti*1i; // coupling
	      auto gd = tr - ti*1i; // coupling
	      cout << "Pairing " << j1 << " " << j2 << " " << g << endl;
	      auto Ci = get_operator(sites,j1,"C"); 
	      auto Cj = get_operator(sites,j2,"C"); 
	      auto Cdagi = get_operator(sites,j1,"Cdag"); 
	      auto Cdagj = get_operator(sites,j2,"Cdag"); 
	      auto Delta = MPO(ampo); // new MPO
	      auto Deltad = MPO(ampo); // new MPO
	      nmultMPO(g*Ci,Cj,Delta,{"Maxm",maxm,"Cutoff",cutoff});
	      nmultMPO(gd*Cdagj,Cdagi,Deltad,{"Maxm",maxm,"Cutoff",cutoff});
	      H = sum(H,Delta); // add contribution
	      H = sum(H,Deltad); // add contribution
      }
    } ;
    jfile.close() ;
    return H ;  // return the Hamiltonian with exchange added
}
;

