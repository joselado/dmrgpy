
float get_energy_fluctuation(auto psi1, auto H) {
    normalize(psi1); // normalize the wavefunction
    auto psi2 = exactApplyMPO(psi1,H,{"Maxm",get_int_value("maxm"),
                 "Cutoff",get_float_value("cutoff")}) ;
    float de;
    de = overlap(psi1,psi2);
    de = overlap(psi2,psi2) - de*de; // energy fluctuation
    return de; // return energy fluctuation
};

// orthogonalize vectors
static auto gram_schmidt=[](auto wfs) {
  for (int i=1;i<wfs.size();i++) {
     for (int j=0;j<i;j++) { // loop over previous
	     auto wf = overlapC(wfs.at(j),wfs.at(i))*wfs.at(j) ;
	     wf = sum_mps(wfs.at(i),-1.0*wf); // remove correction
             normalize(wf) ; // normalize
	     wfs.at(i) = wf ; // redefine wavefunction
     }
  }
  return wfs;
}
;

// function to get several excited states

static auto get_excited=[]() {
  ofstream myfile;
  int nexcited = get_int_value("nexcited") ; // number of excited states
  auto sites = get_sites(); // Get the different sites
  auto H = get_hamiltonian(sites) ; // get the Hamiltonian
  auto sweeps = get_sweeps(); // get the DMRG sweeps
  auto psi0 = MPS(sites); // first wavefunction
  auto en0 = dmrg(psi0,H,sweeps,{"Quiet=",true});
  auto psimax = MPS(sites); // second wavefunction
  int i; 
  auto wfs = std::vector<MPS>(1); // one vector
  auto psi1 = get_gs(); // initialize with the GS
  normalize(psi1);
  wfs.at(0) = psi1; // initialize with the GS
  psi1 = MPS(sites) ; // create new wave
  // lagrange multiplier
  float weight = bandwidth(sites,H)*get_float_value("scale_lagrange_excited"); 
  for (i=1;i<nexcited;i++)  { 
    // now compute a new excited state
    // new energy
    en0 = dmrg(psi1,H,wfs,sweeps,{"Quiet=",true,"Weight=",weight}); 
    normalize(psi1); // normalize the wavefunction
    // compute H|psi>
    // fluctuation in the energy
//    if (de<1e-2) { // if the fluctuation is small enough
      wfs.insert(wfs.end(),psi1); // store this wavefunction
      psi1 = MPS(sites) ; // new random wavefunction
//    };
  } ;
  // perform gram schmidt algorithm
  if (get_bool("excited_gram_schmidt")) wfs = gram_schmidt(wfs); 
  myfile.open("EXCITED.OUT"); // open file
  for (i=0;i<nexcited;i++)  { 
      psi1 = wfs.at(i); // get wavefunction
      auto de = get_energy_fluctuation(psi1,H);
      auto e = overlap(psi1,H,psi1);
      myfile << std::setprecision(20) << e << "  "; // write 
      myfile << std::setprecision(20) << de << "  " << endl; // write 
  };
  myfile.close();
  return wfs ;
}
;


