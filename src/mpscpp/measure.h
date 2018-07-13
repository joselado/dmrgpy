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
  for(int j=1; j <= N; ++j) 
        {
        //re-gauge psi to get ready to measure at position j
        psi.position(j);

        ITensor ket = psi.A(j);
        ITensor bra = dag(prime(ket,Site));

        ITensor Sxjop = sites.op("Sx",j);
        ITensor Syjop = sites.op("Sy",j);
        ITensor Szjop = sites.op("Sz",j);

        //take an inner product 
        auto sxj = (bra*Sxjop*ket).cplx().real();
        auto syj = (bra*Syjop*ket).cplx().real();
        auto szj = (bra*Szjop*ket).cplx().real();
        ofile1 << j << "   "<< sxj << endl ;
        ofile2 << j << "   "<< syj << endl ;
        ofile3 << j << "   "<< szj << endl ;
        }
  ofile1.close() ; 
  ofile2.close() ; 
  ofile3.close() ; 
  return 0 ;
}
