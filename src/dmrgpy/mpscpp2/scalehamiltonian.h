// scale the Hamiltonian so it lies between -1 and 1
static auto scale_hamiltonian=[](auto sites, auto H) {
//    auto psi = MPS(sites); // initialize
//    auto sweeps = get_sweeps(); // get the sweeps
//    auto emin = dmrg(psi,H,sweeps,{"Quiet=",true}); // get minimum energy
//    auto emax = dmrg(psi,-1*H,sweeps,{"Quiet=",true}); // get maximum energy
    auto emin = minimum_energy(sites,H) ;
    auto emax = maximum_energy(sites,H) ;
//    int maxm = get_int_value("maxm"); // bond dimension
//    float cutoff = get_float_value("cutoff"); // cutoff
//    auto args = Args({"Maxm", maxm, "Cutoff",cutoff}); // arguments
//    auto m = sum(-(emin+emax)/2*Iden(sites),H,args); // shift Hamiltonian 
    auto shift = -(emin+emax)/2; // shift in energy
    auto ampo = get_ampo(sites) ; // get the ampo for the Hamiltonian
    ampo += shift,"Id",1; // shift the Hamiltonian
    auto m = MPO(ampo); // generate the shifted Hamiltonian
    ofstream myfile;
    auto scale = (emax-emin)*get_float_value("kpm_scale"); // scale
    scale = 1.0/scale ;  // scale of the Hamiltonian
    myfile.open("KPM_SCALE.OUT");
    myfile << std::setprecision(16) << emin << "   "; // minimum energy
    myfile << std::setprecision(16) << emax << "   "; // maximum energy
    myfile << std::setprecision(16) << scale << endl; // maximum energy
    myfile.close() ; // close file
    m = m*scale ; // scale Hamiltonian
    return m ; // return scaled Hamiltonian
}
;

