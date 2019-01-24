// routine to measure the spin in different sites


int measure(auto psi, auto sites) {
  ofstream ofile1; // declare
  ofstream ofile2; // declare
  ofstream ofile3; // declare
  ofile1.open("MEASURE_SX.OUT"); 
  ofile2.open("MEASURE_SY.OUT"); 
  ofile3.open("MEASURE_SZ.OUT"); 
  int N=sites.N();
  int j ; // counter
  float sxj,syj,szj;
  for(int j=1; j <= N; ++j) 
        {
        //re-gauge psi to get ready to measure at position j
	if (site_type(j-1)==1) { // fermionic site
	  sxj = overlap(psi,get_operator(sites,j-1,"Sx"),psi);
	  syj = overlap(psi,get_operator(sites,j-1,"Sy"),psi);
	  szj = overlap(psi,get_operator(sites,j-1,"Sz"),psi);
	};
	if (site_type(j-1)!=1) { // spin site
          psi.position(j);
          ITensor ket = psi.A(j);
          ITensor bra = dag(prime(ket,Site));
          ITensor Sxjop = sites.op("Sx",j);
          ITensor Syjop = sites.op("Sy",j);
          ITensor Szjop = sites.op("Sz",j);
          //take an inner product 
          sxj = (bra*Sxjop*ket).cplx().real();
          syj = (bra*Syjop*ket).cplx().real();
          szj = (bra*Szjop*ket).cplx().real();
	} ;
        ofile1 << j-1 << "   "<< sxj << endl ;
        ofile2 << j-1 << "   "<< syj << endl ;
        ofile3 << j-1 << "   "<< szj << endl ;
        }
  ofile1.close() ; 
  ofile2.close() ; 
  ofile3.close() ; 
  // measure the occupation and anomalous term
  ofile1.open("MEASURE_N.OUT"); 
  ofile2.open("MEASURE_DELTA.OUT"); 
  ofile3.open("MEASURE_N2.OUT"); 
  for(int j=1; j <= N; ++j) 
        {
        //re-gauge psi to get ready to measure at position j
	if (site_type(j-1)==1) { // fermionic site
	  auto d = overlap(psi,get_occupation_operator(sites,j-1),psi);
	  auto d2 = overlap(psi,get_occupation2_operator(sites,j-1),psi);
	  auto delta = overlapC(psi,get_delta_operator(sites,j-1),psi);
          ofile1 << j-1 << "   "<< d << endl ;
          ofile3 << j-1 << "   "<< d2 << endl ;
          ofile2 << j-1 << "   "<< real(delta) << "  " << imag(delta) << endl ;
	} ;
	} ;
  ofile1.close() ; // close file
  ofile2.close() ; // close file
  ofile3.close() ; // close file
  return 0 ;
}
