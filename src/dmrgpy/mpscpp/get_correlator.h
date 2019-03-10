






// calculate a single correlator
static auto single_correlator=[](auto psi,auto sites, auto i,auto namei, int j, auto namej) {
  
  //Given an MPS or IQMPS called "psi",
  //constructed from a SiteSet "sites"
  
  //Replace "Op1" and "Op2" with the actual names
  //of the operators you want to measure

  auto op_i = sites.op(namei,i);
  auto op_j = sites.op(namej,j);
  
  //below we will assume j > i
  
  //'gauge' the MPS to site i
  //any 'position' between i and j, inclusive, would work here
  psi.position(i); 
  
  
  //index linking i to i+1:
  auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);
  
  auto C = psi.A(i)*op_i*dag(prime(psi.A(i),Site,ir));
  for(int k = i+1; k < j; ++k)
      {
      C *= psi.A(k);
      C *= dag(prime(psi.A(k),Link));
      }
  C *= psi.A(j);
  C *= op_j;
  //index linking j to j-1:
  auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);
  C *= dag(prime(psi.A(j),jl,Site));
  
  auto result = C.cplx().real(); //or C.cplx() if expecting complex
  return result ;

}
;



// calculate all the correlators
int get_correlator_old()   {
  auto sites = get_sites();
  auto H = get_hamiltonian(sites) ;
  auto psi = get_gs(sites,H) ;
  ifstream cfile; // declare
  ofstream ofile; // declare
  cfile.open("correlators.in");  // open file
  ofile.open("CORRELATORS.OUT");  // open file
  int nc; 
  cfile >> nc; // number of correlators
  float c; // declare float
  int i,j ;
  int maxm = get_int_value("maxm") ; // bond dimension for KPM
  auto cutoff = get_float_value("cutoff") ; // cutoff for KPM
  for (int ic=0;ic<nc;++ic) { // loop over correlators
    cfile >> i >> j; // index of correlators
    c = 0.0 ; // initialize
    // get the two operators
    if (site_type(i)>1) { // Spin operators
      auto opi = get_operator(sites,i,get_str("correlator_operator_i")) ;
      auto opj = get_operator(sites,j,get_str("correlator_operator_j")) ;
      auto psi1 = exactApplyMPO(psi,opi,{"Maxm",maxm,"Cutoff",cutoff}) ;
      auto psi2 = exactApplyMPO(psi,opj,{"Maxm",maxm,"Cutoff",cutoff}) ;
      c = overlap(psi1,psi2); // compute the overlap
    }
//       c = overlap(psi,get_sidotsj_operator(sites,i,j),psi) ; // add 
//       c += single_correlator(psi,sites,i+1,"Sx",j+1,"Sx"); // get this one
//       c += single_correlator(psi,sites,i+1,"Sy",j+1,"Sy"); // get this one
//       c += single_correlator(psi,sites,i+1,"Sz",j+1,"Sz"); // get this one
//    } ;
//  BE CAREFUL
//  fermionic operators have strings and cannot be computed in the 
//  previous way, the MPO has to be explicitly created with strings  
    if (site_type(i)<=1) { // fermions (spinful or spinless)
// so far only hopping terms have been implemented
       c = overlap(psi,get_hopping_operator(sites,i,j),psi) ; // add 
    } ;
    ofile << ic << "   " << c << endl ;
  };
  ofile.close() ;
  cfile.close() ;
}





// calculate all the correlators
int get_correlator()   {
  auto sites = get_sites();
  auto H = get_hamiltonian(sites) ;
  auto psi = get_gs(sites,H) ;
  ifstream cfile; // declare
  ofstream ofile; // declare
  cfile.open("correlators.in");  // open file
  ofile.open("CORRELATORS.OUT");  // open file
  int nc; 
  cfile >> nc; // number of correlators
  float c; // declare float
  int i,j ;
  int maxm = get_int_value("maxm") ; // bond dimension for KPM
  float cutoff = get_float_value("cutoff") ; // cutoff for KPM
  for (int ic=0;ic<nc;++ic) { // loop over correlators
    cfile >> i >> j; // index of correlators
    // get the two operators
    auto opi = get_operator(sites,i,get_str("correlator_operator_i")) ;
    auto opj = get_operator(sites,j,get_str("correlator_operator_j")) ;
    auto M = opi*1.0; // new MPO
    nmultMPO(opi,opj, M,{"Maxm",maxm,"Cutoff",cutoff}) ; // multiply MPO
    auto c = overlapC(psi,M,psi);
    ofile << ic << "   " << real(c) << "  " << imag(c) << endl ;
  };
  ofile.close() ;
  cfile.close() ;
}
