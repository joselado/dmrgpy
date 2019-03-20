// compute the KPM moments for matrix m and vectors vi and vj
static auto moments_vi_vj=[](auto m, auto vi, auto vj, int n) {
  // technique use to apply the mpo
//  auto fitmpo = get_bool("fitmpo_kpm") ;
//  fitmpo = false ; // this does not work ok
  ofstream myfile; // file for the moments
  ofstream entropyfile; // file for the entropies
  myfile.open("KPM_MOMENTS.OUT"); // open file
  entropyfile.open("KPM_ENTROPY.OUT"); // open file
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto kpmcutoff = get_float_value("kpm_cutoff") ; // bond dimension for KPM
  auto v = vi*1.0 ; // initialize
  auto am = vi*1.0 ; // initialize
  auto a = vi*0.0 ; // initialize
//  if (fitmpo) fitApplyMPO(v,m,a,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff}) ;
//  else  
  a = exactApplyMPO(v,m,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff}) ;
  auto ap = a*1.0 ; // initialize
  auto bk = overlapC(vj,v) ; // overlap
  auto bk1 = overlapC(vj,a) ; // overlap
  int cindex = v.N()/2 ; // central site
  myfile << std::setprecision(8) << real(bk) << "  "
                       << std::setprecision(8)<< imag(bk) << endl;
  myfile << std::setprecision(8) << real(bk1) << "  "
                       << std::setprecision(8)<< imag(bk1) << endl;
  entropyfile << entropy(v,cindex) << endl ;
  entropyfile << entropy(a,cindex) << endl ;
  int i ;
  for(i=0;i<n;i++) {
//    if (fitmpo) fitApplyMPO(a,m,ap,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff}) ;
//    else 
    ap = exactApplyMPO(a,m,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff}) ;
    ap = sum(2.0*ap,-1.0*am,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff}) ; // recursion relation
    bk = overlapC(vj,ap) ; // compute term 
    myfile << std::setprecision(8) << real(bk) << "  "
                       << std::setprecision(8)<< imag(bk) << endl;
    entropyfile << entropy(ap,cindex) << endl ;
    am = a*1.0; // next iteration
    a = ap*1.0; // next iteration
  } ;
  entropyfile.close();
  myfile.close();
  return 0 ;
} ;


