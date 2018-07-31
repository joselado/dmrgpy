
auto get_exchange(auto ampo) {
    ifstream jfile; // file to read
    jfile.open("couplings.in"); // file with the coupling
    int nj;
    jfile >> nj; // read the number of couplings
    float jxx[nj]; // declare js
    float jxy[nj]; // declare js
    float jxz[nj]; // declare js
    float jyx[nj]; // declare js
    float jyy[nj]; // declare js
    float jyz[nj]; // declare js
    float jzx[nj]; // declare js
    float jzy[nj]; // declare js
    float jzz[nj]; // declare js
    int indexj1[nj]; // declare index
    int indexj2[nj]; // declare index
    for (int i=0;i<nj;++i) {
//      cout << i << endl ;
      jfile >> indexj1[i] >> indexj2[i] >>
      jxx[i] >> jxy[i] >> jxz[i] >> // read this coupling
      jyx[i] >> jyy[i] >> jyz[i] >> // read this coupling
      jzx[i] >> jzy[i] >> jzz[i] ; // read this coupling
    } ;
    jfile.close() ;

    // this function only works for spins so far
    int j1,j2; // indexes
    float jc ; // value
    const Cplx Cplx_i = Cplx(0,1);
    for(int j=0; j < nj; ++j)
        {
        j1 = indexj1[j]+1;
        j2 = indexj2[j]+1;
	// both are spins
	if ((site_type(j1-1)!=1) and (site_type(j2-1)!=1))  {
          ampo += jxx[j],"Sx",j1,"Sx",j2;
          ampo += jxy[j],"Sx",j1,"Sy",j2;
          ampo += jxz[j],"Sx",j1,"Sz",j2;
          ampo += jyx[j],"Sy",j1,"Sx",j2;
          ampo += jyy[j],"Sy",j1,"Sy",j2;
          ampo += jyz[j],"Sy",j1,"Sz",j2;
          ampo += jzx[j],"Sz",j1,"Sx",j2;
          ampo += jzy[j],"Sz",j1,"Sy",j2;
          ampo += jzz[j],"Sz",j1,"Sz",j2;
          } ;
	// one fermion and one spin
	if ((site_type(j1-1)!=1) and (site_type(j2-1)==1))  {
          ampo += jxx[j],"Sx",j1,"Cdagdn",j2,"Cup",j2;
          ampo += jxx[j],"Sx",j1,"Cdagup",j2,"Cdn",j2;
          ampo += -jyy[j],"Sy",j1,"Cdagdn",j2,"CupI",j2;
          ampo += jyy[j],"Sy",j1,"Cdagup",j2,"CdnI",j2;
          ampo += jzz[j],"Sz",j1,"Cdagup",j2,"Cup",j2;
          ampo += -jzz[j],"Sz",j1,"Cdagdn",j2,"Cdn",j2;
          } ;
	// one spin and one fermion, the other case
	if ((site_type(j1-1)==1) and (site_type(j2-1)!=1))  {
          ampo += jxx[j],"Sx",j2,"Cdagdn",j1,"Cup",j1;
          ampo += jxx[j],"Sx",j2,"Cdagup",j1,"Cdn",j1;
          ampo += -jyy[j],"Sy",j2,"Cdagdn",j1,"CupI",j1;
          ampo += jyy[j],"Sy",j2,"Cdagup",j1,"CdnI",j1;
          ampo += jzz[j],"Sz",j2,"Cdagup",j1,"Cup",j1;
          ampo += -jzz[j],"Sz",j2,"Cdagdn",j1,"Cdn",j1;
          } ;
	}
    return ampo ;  // return the Hamiltonian with exchange added
}
