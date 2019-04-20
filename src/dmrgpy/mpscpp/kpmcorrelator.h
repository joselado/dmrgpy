// compute the KPM moments for matrix m and vectors vi and vj
// and a shift in the energy
static auto moments_vi_vj_shift=[](auto m, auto vi, auto vj, int n, auto shift) {
  ofstream myfile;
  myfile.open("KPM_MOMENTS.OUT"); // open file
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto v = vi*1.0 ; // initialize
  auto am = vi*1.0 ; // initialize
  auto a = exactApplyMPO(v,m,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // initialize
  a = sum(a,shift*v,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // shift
  auto ap = a*1.0 ; // initialize
  auto bk = overlapC(vj,v) ; // overlap
  auto bk1 = overlapC(vj,a) ; // overlap
  myfile << std::setprecision(20) << real(bk) << endl;
  myfile << std::setprecision(20) << real(bk1) << endl;
  int i ;
  for(i=0;i<n;i++) {
    ap = exactApplyMPO(a,m,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // apply
    ap = 2.0*sum(ap,shift*a,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // shift
    ap = sum(ap,-1.0*am,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // recursion relation
    bk = overlapC(vj,ap) ; // compute term 
    myfile << std::setprecision(20) << real(bk) << endl;
    am = a*1.0; // next iteration
    a = ap*1.0; // next iteration
  } ;
  return 0 ;
} ;






// compute the Chebyshev polynomials for a certain S+S- correlation
static auto get_moments_spismj_brute=[](auto sites, auto H, int n, int i, int j) {
  auto psi = get_gs(sites,H) ; // get the ground state
  auto m = scale_hamiltonian(sites,H) ; // scale this Hamiltonian
  auto ampo1 = AutoMPO(sites); 
  auto ampo2 = AutoMPO(sites); 
  ampo1 += 1.0,"S-",i ; // S- operator
  ampo2 += 1.0,"S-",j ; // S+ operator
  auto m1 = MPO(ampo1); // first operator
  auto m2 = MPO(ampo2); // second operator
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto psi1 = exactApplyMPO(psi,m1,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ;
  auto psi2 = exactApplyMPO(psi,m2,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ;
  moments_vi_vj(m,psi1,psi2,n) ; //compute the KPM moments
  return 0 ;
} ;





static auto get_moments_dynamical_correlator=[](auto sites, auto H)
{
  auto n = get_int_value("nkpm");
  auto delta = get_float_value("kpm_delta");
  // define on the fly the number of polynomials
  n = int(round(bandwidth(sites,H)/delta))*get_int_value("kpm_n_scale");
  ofstream myfile;
  myfile.open("KPM_NUM_POLYNOMIALS.OUT");
  myfile << std::setprecision(8) << n << endl;
  myfile.close(); // close file
  auto psi = get_gs(sites,H) ; // get the ground state
  auto m = scale_hamiltonian(sites,H) ; // scale this Hamiltonian
  auto m1 = Iden(sites) ; // identity
  auto m2 = Iden(sites) ; // identity
  // now read the operators to use
  // First the operator i
  if (get_bool("kpm_multioperator_i")) {
    m1 = get_multioperator("kpm_multioperator_i");
  } ;
  if (not get_bool("kpm_multioperator_i")) {
    m1 = get_operator(sites,get_int_value("site_i_kpm"),
		    get_str("kpm_operator_i")); // first operator
  } ;
  // afterwards the operator j
  if (get_bool("kpm_multioperator_j")) {
    m2 = get_multioperator("kpm_multioperator_j");
  } ;
  if (not get_bool("kpm_multioperator_j")) {
    m2 = get_operator(sites,get_int_value("site_j_kpm"),
		    get_str("kpm_operator_j")); // first operator
  } ;
  ///////////////////////////////////
  ////// once the operators have been read, continue
  ///////////////////////////////////
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto kpmcutoff = get_float_value("kpm_cutoff") ; // bond dimension for KPM
  if (check_task("orthogonal_kpm")) 
      moments_kpm_ortho(m,psi,m1,m2,n);  //compute the KPM moments
  else  {
    auto psi1 = exactApplyMPO(psi,m1,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff}) ;
    auto psi2 = exactApplyMPO(psi,m2,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff}) ;
    moments_vi_vj(m,psi1,psi2,n) ; } ; //compute the KPM moments
  return 0 ;
} ;








// compute the Chebyshev polynomials for a certain S+S- correlation
// using a smart energy window
static auto get_moments_spismj=[](auto sites, auto H, int n, int i, int j) {
  auto psi = get_gs(sites,H) ; // get the ground state
  float scale = get_float_value("kpm_scale") ; // energy scale
  auto e0 = overlap(psi,H,psi) ; // ground state energy
  auto shift = - e0 - scale/4 ; // shift to the Hamiltonian
  shift = shift/scale ; // now normalize to the scale
  auto m = H*(1.0/scale) ; // scale the Hamiltonian

  // scaling of the energy
  ofstream myfile;
  myfile.open("KPM_SCALE.OUT");
  myfile << std::setprecision(8) << 1.0/scale << endl;
  myfile << std::setprecision(8) << shift << endl;
  myfile.close() ; // close file
  //
  auto ampo1 = AutoMPO(sites); 
  auto ampo2 = AutoMPO(sites); 
  ampo1 += 1.0,"S-",i ; // S- operator
  ampo2 += 1.0,"S-",j ; // S+ operator
  auto m1 = MPO(ampo1); // first operator
  auto m2 = MPO(ampo2); // second operator
  int kpmmaxm = get_int_value("kpmmaxm") ; // bond dimension for KPM
  auto psi1 = exactApplyMPO(psi,m1,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ;
  auto psi2 = exactApplyMPO(psi,m2,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ;
  moments_vi_vj_shift(m,psi1,psi2,n,shift) ; //compute the KPM moments
  return 0 ;
} ;














