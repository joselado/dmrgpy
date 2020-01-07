
// Get the ground state of this Hamiltonian, and write
// the wavefunction into a file

static auto get_gs=[]() {
    // read the GS from a file
    auto sites = get_sites(); // Get the different sites
    auto psi = MPS(sites);
    if (get_bool("gs_from_file"))  {
	    psi = read_wf(get_str("starting_file_gs")) ;
//  check if return this wavefunction 
           if (get_bool("skip_dmrg_gs")) return psi ;
    };
    auto sweeps = get_sweeps(); // get the DMRG sweeps
    auto H = get_hamiltonian(sites) ; // get the Hamiltonian
    double energy = dmrg(psi,H,sweeps); // ground state energy
    ofstream myfile; // create object
    writeToFile("psi_GS.mps",psi); // write the GS wavefunction
    writeToFile("sites.sites",sites); // write the sites
    myfile.open("GS_ENERGY.OUT"); // open file
    myfile << std::setprecision(20) << energy << endl; // write file
//  } ;
  return psi ; // return the ground state
}
;
