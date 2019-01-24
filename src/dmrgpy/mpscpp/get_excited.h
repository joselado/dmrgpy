// function to get several excited states

static auto get_excited=[](auto H, auto sites, auto sweeps, int nexcited) {
  ofstream myfile;
  myfile.open("EXCITED.OUT"); // open file
  auto psi0 = MPS(sites);
  auto en0 = dmrg(psi0,H,sweeps,{"Quiet=",true});
  myfile << std::setprecision(8) << en0 << endl;
  int i; 
  auto wfs = std::vector<MPS>(nexcited);
  for (i=0;i<nexcited;i++) wfs.at(i) = psi0; // initialize with the GS
  auto psi1 = MPS(sites) ; // create new wave
  for (i=1;i<nexcited;i++)  { 
    // now compute a new excited state
    en0 = dmrg(psi1,H,wfs,sweeps,{"Quiet=",true,"Weight=",20.0}); // new energy
    wfs.at(i) = psi1 ; // store this wavefunction
    psi1 = MPS(sites) ; // new random wavefunction
    myfile << std::setprecision(8) << en0 << endl; // write energy 
  } ;
  return wfs ;
}
;

static auto get_new_excited=[](auto H, auto sites, auto sweeps, auto wfs) {
  auto psi1 = MPS(sites);
  auto en1 = dmrg(psi1,H,wfs,sweeps,{"Quiet=",true,"Weight=",20.0});
  return psi1 ;
}
;