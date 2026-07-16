
// Explicit cache replacing the previous hidden per-function statics.
// Kept as a single process-wide instance today (identical lifetime/semantics
// to the statics it replaces: populated on first use, reused after that
// within one mpscpp.x invocation), but named and resettable so a future
// persistent-process backend can give each independent chain/session its
// own instance instead of sharing this one across chains.
struct BandwidthCache
    {
    bool have_min = false; // whether the minimum energy has been computed
    bool have_max = false; // whether the maximum energy has been computed
    double emin = 0.0; // cached minimum energy
    double emax = 0.0; // cached maximum energy
    };

static BandwidthCache bandwidth_cache;

void inline
reset_bandwidth_cache()
    {
    bandwidth_cache = BandwidthCache(); // drop any cached energies
    }

auto minimum_energy=[](auto sites, auto H) {
    if (!bandwidth_cache.have_min) {
//      auto psi = MPS(sites); // initialize
//      auto sweeps = get_sweeps(); // get the sweeps
//      auto emin = dmrg(psi,H,sweeps,{"Quiet=",true}); // get minimum energy
      auto emin = get_gs_energy(H); // get the Hamiltonian
      bandwidth_cache.emin = emin; // saved energy
      bandwidth_cache.have_min = true; // called
    };
    return bandwidth_cache.emin ; // return energy
}
;

auto maximum_energy=[](auto sites, auto H) {
    if (!bandwidth_cache.have_max) {
      auto psi = MPS(sites); // initialize
      auto sweeps = get_sweeps(); // get the sweeps
      auto emax = -dmrg(psi,-1*H,sweeps,{"Quiet=",true}); // get maximum energy
      bandwidth_cache.emax = emax; // saved energy
      bandwidth_cache.have_max = true; // called
    };
    return bandwidth_cache.emax ; // return energy
}
;
// scale the Hamiltonian so it lies between -1 and 1
static auto bandwidth=[](auto sites, auto H) {
//    auto psi = MPS(sites); // initialize
//    auto sweeps = get_sweeps(); // get the sweeps
//    auto emin = dmrg(psi,H,sweeps,{"Quiet=",true}); // get minimum energy
//    auto emax = -dmrg(psi,-1*H,sweeps,{"Quiet=",true}); // get maximum energy
    auto emin = minimum_energy(sites,H) ;
    auto emax = maximum_energy(sites,H) ;
    return emax-emin ; // return the bandwidth
}
;
