#include"kpmortho.h"


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
  myfile << std::setprecision(8) << real(bk) << endl;
  myfile << std::setprecision(8) << real(bk1) << endl;
  int i ;
  for(i=0;i<n;i++) {
    ap = exactApplyMPO(a,m,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // apply
    ap = 2.0*sum(ap,shift*a,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // shift
    ap = sum(ap,-1.0*am,{"Maxm",kpmmaxm,"Cutoff",1E-7}) ; // recursion relation
    bk = overlapC(vj,ap) ; // compute term 
    myfile << std::setprecision(8) << real(bk) << endl;
    am = a*1.0; // next iteration
    a = ap*1.0; // next iteration
  } ;
  return 0 ;
} ;






// compute the Chebyshev moments for the DOS
static auto get_moments_dos=[](auto sites, auto H, int n) {
  auto m = scale_hamiltonian(sites,H) ; // scale this Hamiltonian
  auto psi1 = MPS(sites) ; // random wavefunction
  psi1 = psi1*(1.0/sqrt(overlapC(psi1,psi1))) ;
  moments_vi_vj(m,psi1,psi1,n) ; //compute the KPM moments
  return 0 ;
} ;

// scale the Hamiltonian so it lies between -1 and 1
static auto scale_hamiltonian=[](auto sites, auto H) {
    auto psi = MPS(sites); // initialize
    auto sweeps = get_sweeps(); // get the sweeps
    auto emin = dmrg(psi,H,sweeps,{"Quiet=",true}); // get minimum energy
    auto emax = dmrg(psi,-1*H,sweeps,{"Quiet=",true}); // get maximum energy
    auto dosscale = 0.9/max(abs(emin),abs(emax)) ; // energy scale for DOS
    ofstream myfile;
    myfile.open("KPM_SCALE.OUT");
    myfile << std::setprecision(8) << dosscale << endl;
    myfile.close() ; // close file
    auto m = H*dosscale ; // scale Hamiltonian
    return m ; // return scaled Hamiltonian
}
;


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





static auto get_moments_dynamical_correlator=[](auto sites, auto H, int n,
      int i, int j, auto namei, auto namej) {
  auto psi = get_gs(sites,H) ; // get the ground state
  auto m = scale_hamiltonian(sites,H) ; // scale this Hamiltonian
//  auto ampo1 = AutoMPO(sites); 
//  auto ampo2 = AutoMPO(sites); 
//  ampo1 += 1.0,namei,i ; // operator
//  ampo2 += 1.0,namej,j ; // operator
  cout << "Site " << i << " name " << namei << endl;
  cout << "Site " << j << " name " << namej << endl;
  auto m1 = get_operator(sites,i,namei); // first operator
  auto m2 = get_operator(sites,j,namej); // first operator
//  auto m2 = MPO(ampo2); // second operator
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














