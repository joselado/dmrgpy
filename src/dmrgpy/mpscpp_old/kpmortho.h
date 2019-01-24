// compute the KPM moments for matrix m and vectors vi and vj
// orthogonalize them and compute the overlap over the matrix
int moments_kpm_ortho(auto m, auto wf0, auto A, auto B, int n) {
  ofstream myfile,myfile1;
  auto kpmwf = std::vector<MPS>(n+2) ; // storage for all the KPM vectors
//  myfile.open("KPM_MOMENTS.OUT"); // open file
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto vi = exactApplyMPO(wf0,B,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // initialize
  auto vj = exactApplyMPO(wf0,A,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // initialize
  auto v = vi*1.0 ; // initialize
  auto am = vi*1.0 ; // initialize
  auto a = exactApplyMPO(v,m,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // initialize
  auto ap = a*1.0 ; // initialize
  auto bk = overlapC(vj,v) ; // overlap
  auto bk1 = overlapC(vj,a) ; // overlap
  kpmwf.at(0) = 1.0*v ; // store
  kpmwf.at(1) = 1.0*a ; // store
//  myfile << std::setprecision(8) << real(bk) << "  "
//                       << std::setprecision(8)<< imag(bk) << endl;
//  myfile << std::setprecision(8) << real(bk1) << "  "
//                       << std::setprecision(8)<< imag(bk1) << endl;
  int i,j ;
  for(i=0;i<n;i++) {
    ap = sum(2.0*exactApplyMPO(a,m,{"Maxm",kpmmaxm,"Cutoff",1E-7}),-1.0*am,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // recursion relation
    kpmwf.at(2+i) = 1.0*ap ; // store
    bk = overlapC(vj,ap) ; // compute term 
//    myfile << std::setprecision(8) << real(bk) << "  "
//                       << std::setprecision(8)<< imag(bk) << endl;
    am = a*1.0; // next iteration
    a = ap*1.0; // next iteration
  } ;
//  myfile1.open("KPM_OVERLAPS_IMAG.OUT"); // open file
//  for(i=0;i<n+1;i++) {
//    for(j=0;j<n+1;j++) {
//    bk = overlapC(kpmwf.at(i),kpmwf.at(j)) ; // compute term 
//    myfile << std::setprecision(8) << real(bk) << "  " ;
//    myfile1 << std::setprecision(8) << imag(bk) << "  " ;
//    } ;
  return 0 ;
  // this does not work yet
} ;

