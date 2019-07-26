
static auto  get_exchange=[](auto ampo) {
    ifstream jfile; // file to read
    jfile.open("exchange.in"); // file with the coupling
    int nj;
    jfile >> nj; // read the number of couplings
    std::complex<double> jxx; // declare js
    std::complex<double> jxy; // declare js
    std::complex<double> jxz; // declare js
    std::complex<double> jyx; // declare js
    std::complex<double> jyy; // declare js
    std::complex<double> jyz; // declare js
    std::complex<double> jzx; // declare js
    std::complex<double> jzy; // declare js
    std::complex<double> jzz; // declare js
    int j1,j2; // indexes for the sites
    int c1,c2; // indexes for the components
    auto tr=0.0; // value
    auto ti=0.0; // value
    auto name1 = "Sx";
    auto name2 = "Sx";
    for (int i=0;i<nj;++i) {
//      cout << i << endl ;
      jfile >> j1 >> j2 >> c1 >> c2 >> tr >> ti ; // everything
//      jxx >> jxy >> jxz >> // read this coupling
//      jyx >> jyy >> jyz >> // read this coupling
//      jzx >> jzy >> jzz ; // read this coupling
      j1 +=1 ; // numbering in Itensor starts in 1
      j2 +=1 ; // same
      // both are spins
      if ((site_type(j1-1)!=1) and (site_type(j2-1)!=1))  {
	  if (c1==0) name1 = "Sx" ;
	  if (c1==1) name1 = "Sy" ;
	  if (c1==2) name1 = "Sz" ;
	  if (c2==0) name2 = "Sx" ;
	  if (c2==1) name2 = "Sy" ;
	  if (c2==2) name2 = "Sz" ;
          ampo += tr,name1,j1,name2,j2;
          ampo += ti*1i,name1,j1,name2,j2;
          } ;
	// one fermion and one spin
//      if ((site_type(j1-1)!=1) and (site_type(j2-1)==1))  {
//          ampo += jxx,"Sx",j1,"Cdagdn",j2,"Cup",j2;
//          ampo += jxx,"Sx",j1,"Cdagup",j2,"Cdn",j2;
//          ampo += -jyy*1i,"Sy",j1,"Cdagdn",j2,"Cup",j2;
//          ampo += jyy*1i,"Sy",j1,"Cdagup",j2,"Cdn",j2;
//          ampo += jzz,"Sz",j1,"Cdagup",j2,"Cup",j2;
//          ampo += -jzz,"Sz",j1,"Cdagdn",j2,"Cdn",j2;
//          } ;
//	// one spin and one fermion, the other case
//      if ((site_type(j1-1)==1) and (site_type(j2-1)!=1))  {
//          ampo += jxx,"Sx",j2,"Cdagdn",j1,"Cup",j1;
//          ampo += jxx,"Sx",j2,"Cdagup",j1,"Cdn",j1;
//          ampo += -jyy*1i,"Sy",j2,"Cdagdn",j1,"Cup",j1;
//          ampo += jyy*1i,"Sy",j2,"Cdagup",j1,"Cdn",j1;
//          ampo += jzz,"Sz",j2,"Cdagup",j1,"Cup",j1;
//          ampo += -jzz,"Sz",j2,"Cdagdn",j1,"Cdn",j1;
//          } ;
    } ;
    jfile.close() ;
    return ampo ;  // return the Hamiltonian with exchange added
}
;
