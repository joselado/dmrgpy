// function to get a gap

static auto get_gap = [](auto H, auto sites, auto sweeps) {
  
  
  auto psi0 = MPS(sites);
  
  auto en0 = dmrg(psi0,H,sweeps,{"Quiet=",true});
  
  //
  // Make a vector of previous wavefunctions;
  // code will penalize future wavefunctions
  // for having any overlap with these
  //
  auto wfs = std::vector<MPS>(1);
  wfs.at(0) = psi0;
  
  auto psi1 = MPS(sites);
  
  //
  // Here the Weight option sets the energy penalty for
  // psi1 having any overlap with psi0
  //
  auto en1 = dmrg(psi1,H,wfs,sweeps,{"Quiet=",true,"Weight=",20.0});
  
  //
  // Print the final energies reported by DMRG
  //
//  printfln("\nGround State Energy = %.10f",en0);
 // printfln("\nExcited State Energy = %.10f",en1);
  
//  printfln("\nDMRG energy gap = %.10f",en1-en0);
  //
  // The overlap <psi0|psi1> should be very close to zero
  //
//  printfln("\nOverlap <psi0|psi1> = %.2E",overlap(psi0,psi1));
  ofstream myfile; // create object
  myfile.open("GAP.OUT"); // open file
  myfile << en1-en0 << endl; // write file
  return en1-en0 ;
};
