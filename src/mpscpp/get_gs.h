
// Get the ground state of this Hamiltonian, and write
// the wavefunction into a file

auto get_gs(auto sites, auto H) {
  auto psi = MPS(sites);
  auto sweeps = get_sweeps(); // get the DMRG sweeps
  double energy = dmrg(psi,H,sweeps); // ground state energy
  ofstream myfile; // create object
//  writeToFile("sites_file",sites);
//  writeToFile("psi_GS.dmrg",psi);
//  writeToFile("sites.dmrg",sites);
  myfile.open("GS_ENERGY.OUT"); // open file
  myfile << std::setprecision(8) << energy << endl; // write file
  return psi ; // return the ground state
}
