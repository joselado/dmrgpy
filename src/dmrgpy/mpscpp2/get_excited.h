
float get_energy_fluctuation(auto psi1, auto H) {
    normalize(psi1); // normalize the wavefunction
    auto psi2 = exactApplyMPO(psi1,H,{"Maxm",get_int_value("maxm"),
                 "Cutoff",get_float_value("cutoff")}) ;
    float de;
    de = overlap(psi1,psi2);
    de = overlap(psi2,psi2) - de*de; // energy fluctuation
    return de; // return energy fluctuation
};


// function to get several excited states

static auto get_excited=[]() {
  ofstream myfile;
  int nexcited = get_int_value("nexcited") ; // number of excited states
  auto sites = get_sites(); // Get the different sites
  auto H = get_hamiltonian(sites) ; // get the Hamiltonian
  auto sweeps = get_sweeps(); // get the DMRG sweeps
  myfile.open("EXCITED.OUT"); // open file
  auto psi0 = MPS(sites); // first wavefunction
  auto en0 = dmrg(psi0,H,sweeps,{"Quiet=",true});
  float de=0.0; // fluctuation in the energy
  de = get_energy_fluctuation(psi0,H); // fluctuation of the energy
  auto psimax = MPS(sites); // second wavefunction
  myfile << std::setprecision(20) << en0 << "  " << de << endl;
  int i; 
  auto wfs = std::vector<MPS>(1); // one vector
  wfs.at(0) = psi0; // initialize with the GS
  auto psi1 = MPS(sites) ; // create new wave
  // lagrange multiplier
  float weight = bandwidth(sites,H)*get_float_value("scale_lagrange_excited"); 
  for (i=1;i<nexcited;i++)  { 
    // now compute a new excited state
    // new energy
    en0 = dmrg(psi1,H,wfs,sweeps,{"Quiet=",true,"Weight=",weight}); 
    normalize(psi1); // normalize the wavefunction
    // compute H|psi>
    // fluctuation in the energy
    de = get_energy_fluctuation(psi1,H);
//    if (de<1e-2) { // if the fluctuation is small enough
      wfs.insert(wfs.end(),psi1); // store this wavefunction
      psi1 = MPS(sites) ; // new random wavefunction
      myfile << std::setprecision(20) << en0 << "  " << de << endl; // write 
//    };
  } ;
  return wfs ;
}
;


