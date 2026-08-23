// In-process Python/C++ interface for the ITensor v3 backend -- a straight
// port of mpscpp2/bindings.cc (same Chain session class, same pybind11
// surface exposed to Python) built against ITensor v3 instead of v2. See
// mpscpp3/chain_session.h for the actual v2->v3 API porting notes.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <limits> // std::numeric_limits<double>::quiet_NaN() (gs_energy_generalized's lam0 default)
#include "itensor/all.h"
#include "extra/all.h" // dmrgpy's own extra site types (spin-3/2, Z4, Boson-4, ...)

using namespace itensor;
using namespace std;

#include "check_task.h" // get_bool/get_str/get_int_value/... (get_sites.h needs these)
#include "get_sites.h" // SpinX, incl. the in-memory std::vector<int> constructor
#include "mo_terms.h" // OpFactor/MOTerm + build_ampo/build_mpo
#include "tebd.h" // bond_hamiltonians/build_tebd_gates/tebd_step (Chain::quench_tebd/evolve_and_measure_tebd)
#include "chain_session.h" // Chain: the session/handle model

namespace py = pybind11;

static std::string
dmrgcpp_version()
    {
    return "0.1.0-itensor-v3";
    }

// Exercises real (compiled, not header-only-inline) ITensor code: builds a
// tiny spin-1/2 chain, its identity MPO (the same AutoMPO/MPO construction
// path that get_hamiltonian.h and friends use), and contracts it against
// itself via overlapC -- the same pattern operators.h::trace_mpo uses. If
// pybind11's ABI ever clashed with the vendored ITensor's C++14 build, this
// is where it would surface (a crash or bad value), not just at import time.
static double
itensor_smoke_test(int nsites)
    {
    auto sites = SpinHalf(nsites);
    auto ampo = AutoMPO(sites);
    ampo += 1.0,"Id",1;
    auto Id = toMPO(ampo);
    auto tr = traceC(Id).real();
    return tr;
    }

// Python hands term lists over as plain nested lists/tuples (matching
// MultiOperator.op's existing shape: [coeff, [name,site], ...] per term) --
// pybind11/stl.h + pybind11/complex.h auto-convert those to this pair/vector
// shape with no bindings needed for OpFactor/MOTerm themselves. This little
// adapter converts to the C++-side MOTerm that build_ampo/build_mpo expect.
using PyOpFactor = std::pair<std::string,int>;
using PyTerm = std::pair<std::complex<double>,std::vector<PyOpFactor>>;

static std::vector<MOTerm>
terms_from_python(std::vector<PyTerm> const& pyterms)
    {
    std::vector<MOTerm> out;
    out.reserve(pyterms.size());
    for (auto const& pt : pyterms)
        {
        MOTerm mt;
        mt.coef = pt.first;
        mt.factors.reserve(pt.second.size());
        for (auto const& f : pt.second) mt.factors.push_back(OpFactor{f.first,f.second});
        out.push_back(std::move(mt));
        }
    return out;
    }

PYBIND11_MODULE(_dmrgcpp, m)
    {
    m.doc() = "dmrgpy in-process ITensor extension";
    m.def("version", &dmrgcpp_version,
          "Return the extension module version string");
    m.def("itensor_smoke_test", &itensor_smoke_test, py::arg("nsites") = 4,
          "Build a tiny spin-1/2 chain via the vendored ITensor library and "
          "return a trace, confirming the extension actually links against "
          "and executes ITensor code (not just imports).");
    m.def("set_realify_spin_terms",
          [](bool enabled) { realify_spin_terms_enabled() = enabled; },
          py::arg("enabled"),
          "Enable/disable (default: enabled) the Sy -> (S+ - S-)/(2i) "
          "rewrite every AutoMPO build goes through -- see mo_terms.h. It is "
          "an exact representation change whose only effect is that a real "
          "Hamiltonian written with Sy yields a real-valued MPO (and hence "
          "real DMRG/KPM/TDVP arithmetic) instead of a complex one; process-"
          "global, exposed so a test can check both paths agree numerically.");
    m.def("get_realify_spin_terms",
          []() { return realify_spin_terms_enabled(); },
          "Whether the Sy-elimination rewrite is currently enabled.");

    // module_local(): without it, pybind11 registers "MPS"/"MPO"/"Chain" in
    // a process-wide type registry keyed by typeid() -- and since
    // mpscpp2/bindings.cc's ITensor v2 build defines its own, ABI-incompatible
    // itensor::MPS/MPO/Chain with the *same* mangled names (same
    // unnamespaced C++ types, just a different vendored ITensor copy),
    // libstdc++'s cross-DSO RTTI equality (which falls back to strcmp'ing
    // type names, not just comparing addresses) makes pybind11 think it's
    // the exact same type already registered by mpscpp2 the moment both
    // extensions are imported into the same process -- confirmed directly:
    // without this, `from dmrgpy.mpscpp2 import _dmrgcpp` followed by
    // `from dmrgpy.mpscpp3 import _dmrgcpp` aborts with
    // `generic_type: type "Chain" is already registered!`. module_local()
    // keeps this module's registration in a separate, module-scoped table
    // instead of the shared global one, so it never collides with mpscpp2's
    // (or any other module's) same-named types -- see pybind11's own
    // internals.h (get_local_type_info vs get_global_type_info).
    py::class_<MPS>(m,"MPS",py::module_local()); // opaque handle, no Python-visible methods yet
    py::class_<MPO>(m,"MPO",py::module_local()); // opaque handle, for StaticOperator

    py::class_<Chain>(m,"Chain",py::module_local())
        .def(py::init<std::vector<int> const&>(), py::arg("site_types"))
        .def("set_sweep_params",&Chain::set_sweep_params,
             py::arg("maxm"),py::arg("nsweeps"),py::arg("cutoff"),py::arg("noise"))
        .def("set_mpomaxm",&Chain::set_mpomaxm,py::arg("mpomaxm"))
        .def("set_verbose",&Chain::set_verbose,py::arg("verbose"),
             "Enable/disable ITensor's per-sweep DMRG progress output "
             "(disabled by default)")
        .def("random_mps",&Chain::random_mps)
        .def("set_hamiltonian",[](Chain& self, std::vector<PyTerm> const& terms) {
                self.set_hamiltonian(terms_from_python(terms));
            }, py::arg("terms"))
        .def("gs_energy",&Chain::gs_energy,py::arg("skip_dmrg")=false)
        .def("gs_wavefunction",&Chain::gs_wavefunction,
             py::return_value_policy::copy)
        .def("set_wavefunction",&Chain::set_wavefunction,py::arg("wf"))
        .def("gs_energy_generalized",[](Chain& self, std::vector<PyTerm> const& terms_a,
                                         double lam0) {
                return self.gs_energy_generalized(terms_from_python(terms_a),lam0);
            }, py::arg("terms_a"),
               py::arg("lam0")=std::numeric_limits<double>::quiet_NaN(),
               "Generalized-eigenvalue DMRG: solves H|psi>=lambda*A|psi> for "
               "a Hermitian positive-definite metric operator A (terms_a), "
               "self-consistently updating the Lagrange multiplier lambda "
               "between DMRG sweeps -- see Chain::gs_energy_generalized's "
               "own comment for the algorithm. Returns lambda; A=identity "
               "reduces this exactly to plain gs_energy().")
        .def("excited_states",[](Chain& self, int n, double scale_lagrange,
                                  bool gram_schmidt) {
                auto out = self.excited_states(n,scale_lagrange,gram_schmidt);
                return py::make_tuple(out.energies,out.fluctuations,
                                      out.wavefunctions);
            }, py::arg("n"),py::arg("scale_lagrange")=1.0,
               py::arg("gram_schmidt")=false,
               "Returns (energies, fluctuations, wavefunctions)")
        .def("nhdmrg",[](Chain& self, std::vector<PyTerm> const& terms_h,
                          std::vector<PyTerm> const& terms_hadj,
                          int krylovdim, int restarts) {
                auto out = self.nhdmrg(terms_from_python(terms_h),
                    terms_from_python(terms_hadj),krylovdim,restarts);
                return py::make_tuple(out.energy,out.psil,out.psir);
            }, py::arg("terms_h"),py::arg("terms_hadj"),
               py::arg("krylovdim")=20,py::arg("restarts")=2,
               "Non-Hermitian DMRG (this file is the annotated original; "
               "mpscpp2 carries a v2-API back-port): "
               "optimizes a biorthogonal left/right eigenpair of the "
               "non-Hermitian operator given by terms_h, targeting the "
               "eigenvalue with smallest real part; terms_hadj must be the "
               "adjoint operator's terms (MultiOperator.get_dagger() on the "
               "Python side). Returns (energy, psil, psir)")
        .def("nhdmrg_generalized",[](Chain& self, std::vector<PyTerm> const& terms_h,
                          std::vector<PyTerm> const& terms_hadj,
                          std::vector<PyTerm> const& terms_a,
                          int krylovdim, int restarts, Cplx lam0) {
                auto out = self.nhdmrg_generalized(terms_from_python(terms_h),
                    terms_from_python(terms_hadj),terms_from_python(terms_a),
                    krylovdim,restarts,lam0);
                return py::make_tuple(out.energy,out.psil,out.psir);
            }, py::arg("terms_h"),py::arg("terms_hadj"),py::arg("terms_a"),
               py::arg("krylovdim")=20,py::arg("restarts")=2,
               py::arg("lam0")=Cplx(std::numeric_limits<double>::quiet_NaN(),0.0),
               "Non-Hermitian generalized-eigenvalue NH-DMRG: solves "
               "H|psi_R>=lambda*A|psi_R> for a possibly non-Hermitian "
               "operator (terms_h, with terms_hadj its adjoint -- same "
               "MultiOperator.get_dagger() convention as nhdmrg() above) "
               "and a Hermitian positive-definite metric operator A "
               "(terms_a) -- see Chain::nhdmrg_generalized's own comment "
               "for the algorithm. Returns (lambda, psil, psir) with "
               "<psil|psir>=1, same convention as nhdmrg().")
        .def("idmrg_ground_state",[](Chain& self, std::vector<PyTerm> const& terms_intra,
                                      std::vector<PyTerm> const& terms_inter,
                                      int maxm, double cutoff, int maxiter, double etol,
                                      int krylovdim, int restarts, double noise,
                                      int noise_iters) {
                auto out = self.idmrg_ground_state(terms_from_python(terms_intra),
                    terms_from_python(terms_inter),maxm,cutoff,maxiter,etol,
                    krylovdim,restarts,noise,noise_iters);
                return py::make_tuple(out.density,out.converged,out.niter_done);
            }, py::arg("terms_intra"),py::arg("terms_inter"),
               py::arg("maxm"),py::arg("cutoff"),py::arg("maxiter"),py::arg("etol"),
               py::arg("krylovdim")=30,py::arg("restarts")=2,
               py::arg("noise")=1e-4,py::arg("noise_iters")=40,
               "Infinite DMRG (iDMRG) ground-state energy density -- this "
               "Chain must have been constructed with site_types = the "
               "n_uc-site unit cell (not a full chain), see "
               "Chain::idmrg_ground_state's own comment for the algorithm "
               "and scope (Hermitian, n_uc<=2; fermionic terms ARE "
               "supported). Returns (density, converged, niter_done). "
               "Leaves a converged unit-cell snapshot on this Chain, so "
               "idmrg_onsite_expectation/idmrg_two_point_correlator/"
               "idmrg_local_excitation_gap can be called afterwards.")
        .def("idmrg_onsite_expectation",[](Chain& self, std::string const& opname, int p) {
                return self.idmrg_onsite_expectation(opname,p);
            }, py::arg("opname"),py::arg("p"),
               "<opname> at sublattice p (0..n_uc-1) of a converged "
               "idmrg_ground_state()'s infinite chain -- requires "
               "idmrg_ground_state to have been called first on this same "
               "Chain. See Chain::idmrg_onsite_expectation's own comment.")
        .def("idmrg_two_point_correlator",[](Chain& self, std::string const& opname_i,
                                              int p_i, std::string const& opname_j, int r) {
                return self.idmrg_two_point_correlator(opname_i,p_i,opname_j,r);
            }, py::arg("opname_i"),py::arg("p_i"),py::arg("opname_j"),py::arg("r"),
               "<opname_i(site p_i) opname_j(site p_i+r)> of a converged "
               "idmrg_ground_state()'s infinite chain, r in physical sites "
               "(r>=0) -- requires idmrg_ground_state to have been called "
               "first on this same Chain. Jordan-Wigner strings are threaded "
               "for fermionic operators; see "
               "Chain::idmrg_two_point_correlator's own comment.")
        .def("idmrg_local_excitation_gap",[](Chain& self, int niter) {
                return self.idmrg_local_excitation_gap(niter);
            }, py::arg("niter")=200,
               "The \"local superblock gap\": the second-lowest eigenvalue of "
               "the growing algorithm's own final 2-site effective "
               "Hamiltonian, minus its ground-state eigenvalue -- requires "
               "idmrg_ground_state to have been called first on this same "
               "Chain. A cheap, momentum-less cross-check, not a "
               "variationally optimal excited state; see "
               "Chain::idmrg_local_excitation_gap's own comment.")
        .def("idmrg_local_excitation_gap_detail",[](Chain& self, int niter) {
                auto out = self.idmrg_local_excitation_gap_detail(niter);
                return py::make_tuple(out.gap,out.e0_fresh,out.e0_stored);
            }, py::arg("niter")=200,
               "(gap, e0_fresh, e0_stored) -- idmrg_local_excitation_gap "
               "plus the ground eigenvalue it re-solved for and the one the "
               "growing algorithm's own final micro-step reported. A large "
               "disagreement between the two means the growth loop did not "
               "find its own local effective Hamiltonian's ground state; see "
               "Chain::idmrg_local_excitation_gap_detail's own comment.")
        .def("vms_onsite_expectation",[](Chain& self, std::string const& opname, int p) {
                return self.vms_onsite_expectation(opname,p);
            }, py::arg("opname"),py::arg("p"),
               "<opname> at site p of a converged SEQUENTIAL MULTI-SITE "
               "(n_uc>2) vumps_ground_state -- see Chain::vms_onsite_"
               "expectation. The grouped vumps_onsite_expectation covers "
               "n_uc<=2; the two representations are deliberately separate.")
        .def("vms_two_point_correlator",[](Chain& self, std::string const& opname_i,
                                            int p_i, std::string const& opname_j, int r) {
                return self.vms_two_point_correlator(opname_i,p_i,opname_j,r);
            }, py::arg("opname_i"),py::arg("p_i"),py::arg("opname_j"),py::arg("r"),
               "<opname_i(p_i) opname_j(p_i+r)> of a converged sequential "
               "multi-site (n_uc>2) vumps_ground_state, r in physical sites.")
        .def("vumps_ground_state",[](Chain& self, std::vector<PyTerm> const& terms_intra,
                                      std::vector<PyTerm> const& terms_inter,
                                      int D, double tol, int maxiter, int nrestarts,
                                      int niter_lanczos) {
                auto out = self.vumps_ground_state(terms_from_python(terms_intra),
                    terms_from_python(terms_inter),D,tol,maxiter,nrestarts,niter_lanczos);
                return py::make_tuple(out.e0,out.converged,out.niter_done,out.gauge_mismatch);
            }, py::arg("terms_intra"),py::arg("terms_inter"),
               py::arg("D"),py::arg("tol"),py::arg("maxiter"),py::arg("nrestarts")=4,
               py::arg("niter_lanczos")=60,
               "VUMPS ground state (fixed bond dimension D, thermodynamic "
               "limit) -- this Chain must have been constructed with "
               "site_types = the n_uc-site unit cell, see "
               "Chain::vumps_ground_state's own comment for the algorithm "
               "and scope (Hermitian, n_uc<=2, reach<=1 bonds only). "
               "Returns (e0, converged, niter_done, gauge_mismatch).")
        .def("vumps_excitation_energies",[](Chain& self, double k, int n) {
                return self.vumps_excitation_energies(k,n);
            }, py::arg("k"),py::arg("n")=1,
               "Lowest n tangent-space/quasiparticle excitation energies "
               "(above the ground state) at momentum k (radians per unit "
               "cell) -- requires vumps_ground_state to have been called "
               "first on this same Chain, see "
               "Chain::vumps_excitation_energies's own comment.")
        .def("vumps_onsite_expectation",[](Chain& self, std::string const& opname, int p) {
                return self.vumps_onsite_expectation(opname,p);
            }, py::arg("opname"),py::arg("p"),
               "<opname> at sub-site p (0..n_uc-1) of a converged "
               "vumps_ground_state() -- requires vumps_ground_state to have "
               "been called first on this same Chain, see "
               "Chain::vumps_onsite_expectation's own comment.")
        .def("vumps_two_point_correlator",[](Chain& self, std::string const& opname_i, int p_i,
                                              std::string const& opname_j, int r) {
                return self.vumps_two_point_correlator(opname_i,p_i,opname_j,r);
            }, py::arg("opname_i"),py::arg("p_i"),py::arg("opname_j"),py::arg("r"),
               "<opname_i(site p_i) opname_j(site p_i + r)> of a converged "
               "vumps_ground_state()'s infinite chain, r measured in "
               "physical sites (r>=0) -- requires vumps_ground_state to "
               "have been called first on this same Chain, see "
               "Chain::vumps_two_point_correlator's own comment.")
        .def("vumps_apply_mpo",[](Chain& self, std::vector<std::vector<Cplx>> const& W_bulk_flat,
                                   std::vector<int> const& Dw_left, std::vector<int> const& Dw_right,
                                   double cutoff, int maxdim) {
                auto out = self.vumps_apply_mpo(W_bulk_flat,Dw_left,Dw_right,cutoff,maxdim);
                py::array_t<std::complex<double>> AL({out.D,out.d_g,out.D});
                py::array_t<std::complex<double>> AR({out.D,out.d_g,out.D});
                py::array_t<std::complex<double>> C({out.D,out.D});
                py::array_t<std::complex<double>> AC({out.D,out.d_g,out.D});
                std::copy(out.AL.begin(),out.AL.end(),AL.mutable_data());
                std::copy(out.AR.begin(),out.AR.end(),AR.mutable_data());
                std::copy(out.C.begin(),out.C.end(),C.mutable_data());
                std::copy(out.AC.begin(),out.AC.end(),AC.mutable_data());
                return py::make_tuple(out.D,out.d_g,AL,AR,C,AC,out.eta);
            }, py::arg("W_bulk_flat"),py::arg("Dw_left"),py::arg("Dw_right"),
               py::arg("cutoff"),py::arg("maxdim")=0,
               "Apply a periodic (bounded) MPO to this Chain's own converged "
               "VUMPS snapshot -- requires vumps_ground_state/"
               "vumps_load_uniform_state to have been called first on this "
               "same Chain, see Chain::vumps_apply_mpo's own comment for the "
               "W_bulk_flat convention and the bounded-operator scope "
               "restriction. maxdim<=0 means no cap. Returns "
               "(D,d_g,AL,AR,C,AC,eta), the new mixed-gauge uniform iMPS -- "
               "feed to vumps_load_uniform_state to make onsite_expectation/"
               "two_point_correlator see it.")
        .def("vumps_imps_sum",[](Chain& self, int D_b, std::vector<Cplx> const& AL_b,
                                  double cutoff, int maxdim) {
                auto out = self.vumps_imps_sum(D_b,AL_b,cutoff,maxdim);
                py::array_t<std::complex<double>> AL({out.D,out.d_g,out.D});
                py::array_t<std::complex<double>> AR({out.D,out.d_g,out.D});
                py::array_t<std::complex<double>> C({out.D,out.D});
                py::array_t<std::complex<double>> AC({out.D,out.d_g,out.D});
                std::copy(out.AL.begin(),out.AL.end(),AL.mutable_data());
                std::copy(out.AR.begin(),out.AR.end(),AR.mutable_data());
                std::copy(out.C.begin(),out.C.end(),C.mutable_data());
                std::copy(out.AC.begin(),out.AC.end(),AC.mutable_data());
                return py::make_tuple(out.D,out.d_g,AL,AR,C,AC,out.eta);
            }, py::arg("D_b"),py::arg("AL_b"),py::arg("cutoff"),py::arg("maxdim")=0,
               "Direct sum of this Chain's own converged VUMPS snapshot with "
               "a second, externally-supplied converged state (D_b,AL_b, "
               "flat row-major (D_b,d_g,D_b)) -- requires vumps_ground_state/"
               "vumps_load_uniform_state to have been called first on this "
               "same Chain, see Chain::vumps_imps_sum's own comment for the "
               "physical scope restriction (raises for two ordinary, "
               "equally-normalized VUMPS ground states). maxdim<=0 means no "
               "cap. Returns (D,d_g,AL,AR,C,AC,eta), same convention as "
               "vumps_apply_mpo.")
        .def("vumps_load_uniform_state",[](Chain& self, int D, int d_g,
                                            std::vector<Cplx> const& AL,
                                            std::vector<Cplx> const& AR,
                                            std::vector<Cplx> const& C) {
                self.vumps_load_uniform_state(D,d_g,AL,AR,C);
            }, py::arg("D"),py::arg("d_g"),py::arg("AL"),py::arg("AR"),py::arg("C"),
               "Loads an externally-computed mixed-gauge uniform iMPS (e.g. "
               "vumps_apply_mpo/vumps_imps_sum's own AL/AR/C output) into "
               "this Chain's own VUMPS snapshot, so vumps_onsite_expectation/"
               "vumps_two_point_correlator/a further apply_mpo/imps_sum call "
               "see it -- see Chain::vumps_load_uniform_state's own comment.")
        .def("vumps_get_snapshot",[](Chain& self) {
                int D=0, d_g=0; std::vector<Cplx> AL, AR, C;
                self.vumps_get_snapshot(D,d_g,AL,AR,C);
                py::array_t<std::complex<double>> AL_arr({D,d_g,D});
                py::array_t<std::complex<double>> AR_arr({D,d_g,D});
                py::array_t<std::complex<double>> C_arr({D,D});
                std::copy(AL.begin(),AL.end(),AL_arr.mutable_data());
                std::copy(AR.begin(),AR.end(),AR_arr.mutable_data());
                std::copy(C.begin(),C.end(),C_arr.mutable_data());
                return py::make_tuple(D,d_g,AL_arr,AR_arr,C_arr);
            }, "Reads back this Chain's own current VUMPS snapshot AL/AR/C "
               "-- requires vumps_ground_state/vumps_load_uniform_state to "
               "have been called first, see Chain::vumps_get_snapshot's own "
               "comment. Returns (D,d_g,AL,AR,C).")
        .def("td_dynamical_correlator_window",
             [](Chain& self, int n_window, std::string const& opname_A,
                std::string const& opname_B, double dt, int nt,
                std::vector<int> const& x_values, int maxdim, double cutoff,
                int niter, bool connected, int p_i) {
                auto out = self.td_dynamical_correlator_window(
                    n_window,opname_A,opname_B,dt,nt,x_values,maxdim,cutoff,
                    niter,connected,p_i);
                py::array_t<double> ts_arr(out.ts.size());
                std::copy(out.ts.begin(),out.ts.end(),ts_arr.mutable_data());
                py::array_t<int> xs_arr(out.xs.size());
                std::copy(out.xs.begin(),out.xs.end(),xs_arr.mutable_data());
                py::array_t<std::complex<double>> S_arr({(int)out.ts.size(),(int)out.xs.size()});
                std::copy(out.S.begin(),out.S.end(),S_arr.mutable_data());
                return py::make_tuple(ts_arr,xs_arr,S_arr);
            }, py::arg("n_window"),py::arg("opname_A"),py::arg("opname_B"),
               py::arg("dt"),py::arg("nt"),py::arg("x_values"),
               py::arg("maxdim"),py::arg("cutoff"),py::arg("niter"),
               py::arg("connected")=true,py::arg("p_i")=0,
               "Real-time IBC-window dynamical correlator S(x,t) (Milsted/"
               "Vanderstraeten, arXiv:1804.09163, Sec. V.1) -- requires "
               "idmrg_ground_state to have been called first on this same "
               "Chain (uses its own converged environment snapshot, see "
               "Chain::td_dynamical_correlator_window's own comment). "
               "Returns (ts, xs, S), S shaped (len(ts),len(xs)) complex.")
        .def("vev",[](Chain& self, std::vector<PyTerm> const& terms,
                       MPS const& wf, int npow) {
                return self.vev(terms_from_python(terms),wf,npow);
            }, py::arg("terms"),py::arg("wf"),py::arg("npow")=1)
        .def("apply_operator",[](Chain& self, std::vector<PyTerm> const& terms,
                                  MPS const& wf) {
                return self.apply_operator(terms_from_python(terms),wf);
            }, py::arg("terms"),py::arg("wf"))
        .def("correlation_matrix",[](Chain& self, MPS const& wf) {
                int n = self.num_sites();
                auto flat = self.correlation_matrix(wf);
                py::array_t<std::complex<double>> arr({n,n});
                std::copy(flat.begin(),flat.end(),arr.mutable_data());
                return arr;
            }, py::arg("wf"), "Returns <Cdag_i C_j> as an (N,N) array")
        .def("four_correlation_tensor",
            [](Chain& self, MPS const& wf, bool accelerate) {
                int n = self.num_sites();
                auto flat = self.four_correlation_tensor(wf,accelerate);
                py::array_t<std::complex<double>> arr({n,n,n,n});
                std::copy(flat.begin(),flat.end(),arr.mutable_data());
                return arr;
            }, py::arg("wf"), py::arg("accelerate")=true,
               "Returns <Cdag_i C_j Cdag_k C_l> as an (N,N,N,N) array")
        .def("four_correlation_tensor_spinful",
            [](Chain& self, MPS const& wf, bool accelerate) {
                int n = 2*self.num_sites(); // flat modes: 2 per native site
                auto flat = self.four_correlation_tensor_spinful(wf,accelerate);
                py::array_t<std::complex<double>> arr({n,n,n,n});
                std::copy(flat.begin(),flat.end(),arr.mutable_data());
                return arr;
            }, py::arg("wf"), py::arg("accelerate")=true,
               "Native-spinful-site (Electron/Hubbard) version of "
               "four_correlation_tensor: returns <Cdag_i C_j Cdag_k C_l> "
               "as a (2N,2N,2N,2N) array over flat fermionic modes "
               "(mode 2*s=up, 2*s+1=down at physical site s)")
        .def("four_correlation_tensor_sweep",
            [](Chain& self, MPS const& wf, bool accelerate) {
                int n = self.num_sites();
                auto flat = self.four_correlation_tensor_sweep(wf,accelerate);
                py::array_t<std::complex<double>> arr({n,n,n,n});
                std::copy(flat.begin(),flat.end(),arr.mutable_data());
                return arr;
            }, py::arg("wf"), py::arg("accelerate")=true,
               "Single-sweep, environment-reuse version of "
               "four_correlation_tensor (see chain_session.h's "
               "four_correlation_tensor_sweep for the algorithm, following "
               "ITensorCorrelators.jl): returns the same <Cdag_i C_j Cdag_k "
               "C_l> (N,N,N,N) array, computed without rebuilding a fresh "
               "AutoMPO per (i,j,k,l) tuple")
        .def("four_correlation_tensor_fold",
            [](Chain& self, MPS const& wf,
               std::vector<std::string> const& cdag_names,
               std::vector<int> const& cdag_sites,
               std::vector<std::string> const& c_names,
               std::vector<int> const& c_sites, bool accelerate) {
                int nm = int(c_names.size());
                auto flat = self.four_correlation_tensor_fold(
                    wf,cdag_names,cdag_sites,c_names,c_sites,accelerate);
                py::array_t<std::complex<double>> arr({nm,nm,nm,nm});
                std::copy(flat.begin(),flat.end(),arr.mutable_data());
                return arr;
            }, py::arg("wf"), py::arg("cdag_names"), py::arg("cdag_sites"),
               py::arg("c_names"), py::arg("c_sites"),
               py::arg("accelerate")=true,
               "<Cdag_i C_j Cdag_k C_l> over flat fermionic modes by local "
               "operator folds -- no MPO per tuple. Covers spinless and "
               "native spinful (Electron) sites alike; see chain_session.h")
        .def("kpm_dynamical_correlator",
            [](Chain& self, std::vector<PyTerm> const& terms_i,
               std::vector<PyTerm> const& terms_j, int kpmmaxm,
               double kpm_scale, bool kpm_accelerate, int kpm_n_scale,
               double delta, double kpm_cutoff) {
                auto out = self.kpm_dynamical_correlator(
                    terms_from_python(terms_i),terms_from_python(terms_j),
                    kpmmaxm,kpm_scale,kpm_accelerate,kpm_n_scale,delta,kpm_cutoff);
                return py::make_tuple(out.moments,out.emin,out.emax,
                                      out.scale,out.num_polynomials);
            }, py::arg("terms_i"),py::arg("terms_j"),py::arg("kpmmaxm"),
               py::arg("kpm_scale"),py::arg("kpm_accelerate"),
               py::arg("kpm_n_scale"),py::arg("delta"),py::arg("kpm_cutoff"),
               "Returns (moments, emin, emax, scale, num_polynomials)")
        .def("general_kpm",[](Chain& self, std::vector<PyTerm> const& terms_x,
                               MPS const& wfa, MPS const& wfb, int kpmmaxm,
                               bool kpm_accelerate, int num_polynomials,
                               double kpm_cutoff) {
                return self.general_kpm(terms_from_python(terms_x),wfa,wfb,
                    kpmmaxm,kpm_accelerate,num_polynomials,kpm_cutoff);
            }, py::arg("terms_x"),py::arg("wfa"),py::arg("wfb"),
               py::arg("kpmmaxm"),py::arg("kpm_accelerate"),
               py::arg("num_polynomials"),py::arg("kpm_cutoff"))
        .def("nhkpm_moments",[](Chain& self, std::vector<PyTerm> const& terms_hs,
                               std::vector<PyTerm> const& terms_hs_dag,
                               MPS const& wfa, MPS const& wfb, int n,
                               int kpmmaxm, double kpmcutoff) {
                return self.nhkpm_moments(terms_from_python(terms_hs),
                    terms_from_python(terms_hs_dag),wfa,wfb,n,kpmmaxm,kpmcutoff);
            }, py::arg("terms_hs"),py::arg("terms_hs_dag"),py::arg("wfa"),
               py::arg("wfb"),py::arg("n"),py::arg("kpmmaxm"),py::arg("kpmcutoff"),
               "Non-Hermitian KPM biorthogonal Chebyshev moments; "
               "terms_hs/terms_hs_dag are the z-shifted, rescaled operator "
               "(z*Id-H)/E_max and its dagger, see nonhermitian/kpm.py")
        .def("kpm_dynamical_correlator_truncated",
            [](Chain& self, std::vector<PyTerm> const& terms_i,
               std::vector<PyTerm> const& terms_j, int kpmmaxm,
               double kpm_scale, bool kpm_accelerate, int kpm_n_scale,
               double delta, double kpm_cutoff, int trunc_dK,
               int trunc_nsweeps, double trunc_threshold) {
                auto out = self.kpm_dynamical_correlator_truncated(
                    terms_from_python(terms_i),terms_from_python(terms_j),
                    kpmmaxm,kpm_scale,kpm_accelerate,kpm_n_scale,delta,kpm_cutoff,
                    trunc_dK,trunc_nsweeps,trunc_threshold);
                return py::make_tuple(out.moments,out.emin,out.emax,
                                      out.scale,out.num_polynomials);
            }, py::arg("terms_i"),py::arg("terms_j"),py::arg("kpmmaxm"),
               py::arg("kpm_scale"),py::arg("kpm_accelerate"),
               py::arg("kpm_n_scale"),py::arg("delta"),py::arg("kpm_cutoff"),
               py::arg("trunc_dK"),py::arg("trunc_nsweeps"),py::arg("trunc_threshold"),
               "Energy-truncated KPM dynamical correlator (Holzner et al. "
               "PRB 83, 195115 (2011), Sec. III-B), a wholly independent "
               "method from kpm_dynamical_correlator -- see chain_session.h. "
               "Returns (moments, emin, emax, scale, num_polynomials)")
        .def("scaled_hamiltonian_gs_anchored",
            [](Chain& self, double kpm_scale, double safety) {
                auto out = self.scaled_hamiltonian_gs_anchored(kpm_scale,safety);
                return py::make_tuple(out.scaled_H,out.emin,out.emax,out.scale);
            }, py::arg("kpm_scale"),py::arg("safety")=0.025,
               "Ground-state-anchored KPM rescaling (Eq. 21b) -- exposed "
               "mainly for direct testing of kpm_energy_truncate(); see "
               "chain_session.h for the full derivation. Returns "
               "(scaled_H, emin, emax, scale)")
        .def("kpm_energy_truncate",
            [](Chain& self, MPS const& psi, MPO const& H, int dK,
               int n_sweeps, double threshold) {
                auto out = self.kpm_energy_truncate(psi,H,dK,n_sweeps,threshold);
                return py::make_tuple(out.first,out.second.avg_truncated_weight,
                                      out.second.state_change_norm);
            }, py::arg("psi"),py::arg("H"),py::arg("dK"),py::arg("n_sweeps"),
               py::arg("threshold"),
               "Project high-rescaled-energy components out of a Chebyshev "
               "vector (Holzner et al. PRB 83, 195115 (2011), Sec. III-B); "
               "see chain_session.h's own comment for the algorithm. `H` "
               "must be the same rescaled Hamiltonian used to build `psi`'s "
               "own Chebyshev recursion. Returns (new_psi, "
               "avg_truncated_weight, state_change_norm)")
        .def("reduced_dm",[](Chain& self, MPS const& wf, int site) {
                auto flat = self.reduced_dm(wf,site);
                int dim = self.site_dim(site);
                py::array_t<std::complex<double>> arr({dim,dim});
                std::copy(flat.begin(),flat.end(),arr.mutable_data());
                return arr;
            }, py::arg("wf"),py::arg("site"),
               "Returns the reduced density matrix at a (1-based) site as a "
               "(dim,dim) array")
        .def("exponential_apply",[](Chain& self, std::vector<PyTerm> const& terms,
                                     MPS const& wf, std::complex<double> tau,
                                     int nsteps) {
                return self.exponential_apply(terms_from_python(terms),wf,tau,nsteps);
            }, py::arg("terms"),py::arg("wf"),py::arg("tau"),py::arg("nsteps"),
               "Applies exp(tau*H) to wf via nsteps repeated 2nd-order Taylor steps")
        .def("set_hamiltonian_mpo",&Chain::set_hamiltonian_mpo,py::arg("H"),
               "Set the Hamiltonian from an ALREADY-BUILT MPO (build_operator's "
               "own return value, or anything multiply_operators/sum_operators/"
               "scale_operator produced from one) instead of from a symbolic "
               "term list -- see Chain::set_hamiltonian_mpo's own comment for "
               "when that is worth doing.")
        .def("build_operator",[](Chain& self, std::vector<PyTerm> const& terms) {
                return self.build_operator(terms_from_python(terms));
            }, py::arg("terms"))
        .def("apply_pure_operator",&Chain::apply_pure_operator,
             py::arg("A"),py::arg("wf"))
        .def("multiply_operators",&Chain::multiply_operators,
             py::arg("A"),py::arg("B"))
        .def("sum_operators",&Chain::sum_operators,
             py::arg("A"),py::arg("B"))
        .def("scale_operator",&Chain::scale_operator,
             py::arg("A"),py::arg("z"))
        .def("trace_operator",&Chain::trace_operator,py::arg("A"))
        .def("hermitian_operator",&Chain::hermitian_operator,py::arg("A"))
        .def("overlap_aMb_operator",&Chain::overlap_aMb_operator,
             py::arg("wf1"),py::arg("A"),py::arg("wf2"))
        .def("bond_entropy",&Chain::bond_entropy,py::arg("wf"),py::arg("b"))
        .def("quench",[](Chain& self, std::vector<PyTerm> const& terms_h,
                          std::vector<PyTerm> const& terms_i,
                          std::vector<PyTerm> const& terms_j,
                          int nt, double dt, bool fit_td) {
                auto out = self.quench(terms_from_python(terms_h),
                    terms_from_python(terms_i),terms_from_python(terms_j),
                    nt,dt,fit_td);
                return py::make_tuple(out.correlator,out.final_wf);
            }, py::arg("terms_h"),py::arg("terms_i"),py::arg("terms_j"),
               py::arg("nt"),py::arg("dt"),py::arg("fit_td")=true,
               "Returns (correlator, final_wf)")
        .def("quench_tdvp",[](Chain& self, std::vector<PyTerm> const& terms_h,
                          std::vector<PyTerm> const& terms_i,
                          std::vector<PyTerm> const& terms_j,
                          int nt, double dt) {
                auto out = self.quench_tdvp(terms_from_python(terms_h),
                    terms_from_python(terms_i),terms_from_python(terms_j),
                    nt,dt);
                return py::make_tuple(out.correlator,out.final_wf);
            }, py::arg("terms_h"),py::arg("terms_i"),py::arg("terms_j"),
               py::arg("nt"),py::arg("dt"),
               "TDVP counterpart of quench(). Returns (correlator, final_wf)")
        .def("evolve_and_measure_tdvp",
            [](Chain& self, std::vector<PyTerm> const& terms_h,
               std::vector<PyTerm> const& terms_op, MPS const& wf,
               int nt, double dt) {
                auto out = self.evolve_and_measure_tdvp(terms_from_python(terms_h),
                    terms_from_python(terms_op),wf,nt,dt);
                return py::make_tuple(out.correlator,out.final_wf);
            }, py::arg("terms_h"),py::arg("terms_op"),py::arg("wf"),
               py::arg("nt"),py::arg("dt"),
               "TDVP counterpart of evolve_and_measure(). Returns (correlator, final_wf)")
        .def("quench_tdvp_gse",[](Chain& self, std::vector<PyTerm> const& terms_h,
                          std::vector<PyTerm> const& terms_i,
                          std::vector<PyTerm> const& terms_j,
                          int nt, double dt, int gse_sweeps, int krylov_order,
                          double gse_cutoff) {
                auto out = self.quench_tdvp_gse(terms_from_python(terms_h),
                    terms_from_python(terms_i),terms_from_python(terms_j),
                    nt,dt,gse_sweeps,krylov_order,gse_cutoff);
                return py::make_tuple(out.correlator,out.final_wf);
            }, py::arg("terms_h"),py::arg("terms_i"),py::arg("terms_j"),
               py::arg("nt"),py::arg("dt"),py::arg("gse_sweeps"),
               py::arg("krylov_order"),py::arg("gse_cutoff"),
               "One-site-TDVP-with-global-subspace-expansion counterpart of "
               "quench_tdvp(). Returns (correlator, final_wf)")
        .def("evolve_and_measure_tdvp_gse",
            [](Chain& self, std::vector<PyTerm> const& terms_h,
               std::vector<PyTerm> const& terms_op, MPS const& wf,
               int nt, double dt, int gse_sweeps, int krylov_order,
               double gse_cutoff) {
                auto out = self.evolve_and_measure_tdvp_gse(terms_from_python(terms_h),
                    terms_from_python(terms_op),wf,nt,dt,gse_sweeps,krylov_order,
                    gse_cutoff);
                return py::make_tuple(out.correlator,out.final_wf);
            }, py::arg("terms_h"),py::arg("terms_op"),py::arg("wf"),
               py::arg("nt"),py::arg("dt"),py::arg("gse_sweeps"),
               py::arg("krylov_order"),py::arg("gse_cutoff"),
               "One-site-TDVP-with-global-subspace-expansion counterpart of "
               "evolve_and_measure_tdvp(). Returns (correlator, final_wf)")
        .def("quench_tebd",[](Chain& self, std::vector<PyTerm> const& terms_h,
                          std::vector<PyTerm> const& terms_i,
                          std::vector<PyTerm> const& terms_j,
                          int nt, double dt) {
                auto out = self.quench_tebd(terms_from_python(terms_h),
                    terms_from_python(terms_i),terms_from_python(terms_j),
                    nt,dt);
                return py::make_tuple(out.correlator,out.final_wf);
            }, py::arg("terms_h"),py::arg("terms_i"),py::arg("terms_j"),
               py::arg("nt"),py::arg("dt"),
               "2nd-order-Trotter TEBD counterpart of quench()/quench_tdvp(). "
               "Only valid for a strictly nearest-neighbor terms_h (raises "
               "for any term spanning 3+ distinct sites). "
               "Returns (correlator, final_wf)")
        .def("evolve_and_measure_tebd",
            [](Chain& self, std::vector<PyTerm> const& terms_h,
               std::vector<PyTerm> const& terms_op, MPS const& wf,
               int nt, double dt) {
                auto out = self.evolve_and_measure_tebd(terms_from_python(terms_h),
                    terms_from_python(terms_op),wf,nt,dt);
                return py::make_tuple(out.correlator,out.final_wf);
            }, py::arg("terms_h"),py::arg("terms_op"),py::arg("wf"),
               py::arg("nt"),py::arg("dt"),
               "TEBD counterpart of evolve_and_measure()/evolve_and_measure_tdvp(). "
               "Only valid for a strictly nearest-neighbor terms_h. "
               "Returns (correlator, final_wf)")
        .def("tdvp_step",&Chain::tdvp_step,
             py::arg("H"),py::arg("wf"),py::arg("dt"),py::arg("num_center")=2,
             py::arg("niter")=50,
             "One TDVP step of size dt (may be complex -- see "
             "TDVP/README.md's own \"t\" convention) applied to an "
             "already-built MPO H and MPS wf, one-site (num_center=1) or "
             "two-site (num_center=2, the default). Lets a caller drive "
             "the evolution one variable-sized step at a time, unlike "
             "quench_tdvp/evolve_and_measure_tdvp above, which loop "
             "internally over a fixed number of equal, real dt steps "
             "(used by the \"TDZ\" complex-time-evolution dynamical "
             "correlator, see tdz.py, whose per-step contour increment "
             "varies with t). One-site TDVP doesn't grow bond dimension "
             "on its own -- pair with global_subspace_expand() below.")
        .def("metts_vev",[](Chain& self, std::vector<std::vector<PyTerm>> const& terms_ops,
                             double T, int nsamples, int nwarmup,
                             double dbeta_half_step,
                             std::vector<std::string> const& basis_ops,
                             unsigned long seed, int niter) {
                std::vector<std::vector<MOTerm>> terms_ops_cpp;
                terms_ops_cpp.reserve(terms_ops.size());
                for (auto const& terms_op : terms_ops)
                    terms_ops_cpp.push_back(terms_from_python(terms_op));
                return self.metts_vev(terms_ops_cpp,T,nsamples,
                    nwarmup,dbeta_half_step,basis_ops,seed,niter);
            }, py::arg("terms_ops"),py::arg("T"),py::arg("nsamples")=200,
               py::arg("nwarmup")=20,py::arg("dbeta_half_step")=0.05,
               py::arg("basis_ops")=std::vector<std::string>{"Sz","Sx"},
               py::arg("seed")=0,py::arg("niter")=30,
               "Finite-temperature <A> for one or more operators via METTS "
               "sampling (White & Stoudenmire, arXiv:1002.1305 -- see "
               "Chain::metts_vev's own comment and pyitensor/metts.py for "
               "the full algorithm). terms_ops is a list of operators, all "
               "measured on the same sampled Markov chain. Returns "
               "(means, stderrs), lists matching terms_ops.")
        .def("metts_dynamical_correlator",
            [](Chain& self, std::vector<PyTerm> const& terms_a,
               std::vector<PyTerm> const& terms_b,
               double T, int nt, double dt, int nsamples, int nwarmup,
               double dbeta_half_step,
               std::vector<std::string> const& basis_ops,
               unsigned long seed, int niter, int tdvp_niter) {
                return self.metts_dynamical_correlator(
                    terms_from_python(terms_a),terms_from_python(terms_b),
                    T,nt,dt,nsamples,nwarmup,dbeta_half_step,basis_ops,
                    seed,niter,tdvp_niter);
            }, py::arg("terms_a"),py::arg("terms_b"),py::arg("T"),
               py::arg("nt")=200,py::arg("dt")=0.1,py::arg("nsamples")=100,
               py::arg("nwarmup")=20,py::arg("dbeta_half_step")=0.05,
               py::arg("basis_ops")=std::vector<std::string>{"Sz","Sx"},
               py::arg("seed")=0,py::arg("niter")=30,py::arg("tdvp_niter")=50,
               "Real-time finite-temperature correlator C_AB(t)=<A(t)B>_T "
               "via dynamical METTS sampling (Wang, McClarty, Dankova, "
               "Honecker & Wietek, arXiv:2405.18484, Sec. II -- see "
               "Chain::metts_dynamical_correlator's own comment and "
               "pyitensor/metts.py's metts_dynamical_correlator for the "
               "full algorithm). Returns (means, stderrs), length-nt "
               "arrays over t=0,dt,...,(nt-1)*dt.")
        .def("global_subspace_expand",&Chain::global_subspace_expand,
             py::arg("H"),py::arg("phi"),py::arg("krylov_order"),
             py::arg("cutoff"),py::arg("maxdim")=0,
             "Krylov-subspace global subspace expansion (arXiv:2005.06104), "
             "growing phi's bond dimension using H so one-site TDVP can "
             "keep up with two-site TDVP's own SVD-driven growth.")
        .def("evolve_taylor_step",&Chain::evolve_taylor_step,
             py::arg("H"),py::arg("wf"),py::arg("z"),
             "Applies one Taylor-expanded exp(z*H) step (z may be "
             "complex) to an already-built MPO H and MPS wf -- the "
             "MPO-Taylor (non-TDVP) analogue of tdvp_step() above, used "
             "by \"TDZ\" (tdz.py) as a cross-check / non-TDVP "
             "alternative here, and as mpscpp2's only route to TDZ "
             "(mpscpp2 has no TDVP).")
        .def("evolve_and_measure",
            [](Chain& self, std::vector<PyTerm> const& terms_h,
               std::vector<PyTerm> const& terms_op, MPS const& wf,
               int nt, double dt, bool fit_td) {
                auto out = self.evolve_and_measure(terms_from_python(terms_h),
                    terms_from_python(terms_op),wf,nt,dt,fit_td);
                return py::make_tuple(out.correlator,out.final_wf);
            }, py::arg("terms_h"),py::arg("terms_op"),py::arg("wf"),
               py::arg("nt"),py::arg("dt"),py::arg("fit_td")=true,
               "Returns (correlator, final_wf)")
        .def("overlap",&Chain::overlap_mps,py::arg("wf1"),py::arg("wf2"))
        .def("overlap_aMb",[](Chain& self, MPS const& wf1,
                               std::vector<PyTerm> const& terms, MPS const& wf2) {
                return self.overlap_aMb(wf1,terms_from_python(terms),wf2);
            }, py::arg("wf1"),py::arg("terms"),py::arg("wf2"))
        .def("sum_mps",&Chain::sum_mps,py::arg("wf1"),py::arg("wf2"))
        .def("conjugate",&Chain::conjugate,py::arg("wf"))
        .def("cvm_dynamical_correlator",
            [](Chain& self, std::vector<PyTerm> const& terms_i,
               std::vector<PyTerm> const& terms_j, double omega, double eta,
               double energy, double tol, int max_it) {
                return self.cvm_dynamical_correlator(
                    terms_from_python(terms_i),terms_from_python(terms_j),
                    omega,eta,energy,tol,max_it);
            }, py::arg("terms_i"),py::arg("terms_j"),py::arg("omega"),
               py::arg("eta"),py::arg("energy"),py::arg("tol"),py::arg("max_it"))
        .def("apply_inverse",[](Chain& self, std::vector<PyTerm> const& terms,
                                 MPS const& wf, double tol, int max_it) {
                return self.apply_inverse(terms_from_python(terms),wf,tol,max_it);
            }, py::arg("terms"),py::arg("wf"),py::arg("tol"),py::arg("max_it"))
        ;
    }
