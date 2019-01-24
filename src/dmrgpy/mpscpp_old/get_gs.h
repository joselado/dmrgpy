
// Get the ground state of this Hamiltonian, and write
// the wavefunction into a file

auto get_gs(auto sites, auto H) {
    // read the GS from a file
    auto psi = MPS(sites);
    if (get_bool("gs_from_file")) 
	    psi = read_wf(get_str("starting_file_gs")) ;
//  if (check_task("restart")) { // restart the calculation
//    psi = read_wf() ; // read the wavefunction from a file
//  } 
//  else { // perform the calculation of the GS
    auto sweeps = get_sweeps(); // get the DMRG sweeps
    double energy = dmrg(psi,H,sweeps); // ground state energy
    ofstream myfile; // create object
    writeToFile("psi_GS.mps",psi); // write the GS wavefunction
    writeToFile("sites.sites",sites); // write the sites
    myfile.open("GS_ENERGY.OUT"); // open file
    myfile << std::setprecision(8) << energy << endl; // write file
//  } ;
  return psi ; // return the ground state
}
