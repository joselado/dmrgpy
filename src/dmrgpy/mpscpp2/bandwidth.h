// scale the Hamiltonian so it lies between -1 and 1
static auto bandwidth=[](auto sites, auto H) {
    auto psi = MPS(sites); // initialize
    auto sweeps = get_sweeps(); // get the sweeps
    auto emin = dmrg(psi,H,sweeps,{"Quiet=",true}); // get minimum energy
    auto emax = -dmrg(psi,-1*H,sweeps,{"Quiet=",true}); // get maximum energy
    return emax-emin ; // return the bandwidth
}
;


