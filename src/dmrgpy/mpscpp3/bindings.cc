// In-process Python/C++ interface for the ITensor v3 backend -- a straight
// port of mpscpp2/bindings.cc (same Chain session class, same pybind11
// surface exposed to Python) built against ITensor v3 instead of v2. See
// mpscpp3/chain_session.h for the actual v2->v3 API porting notes.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include "itensor/all.h"
#include "extra/all.h" // dmrgpy's own extra site types (spin-3/2, Z4, Boson-4, ...)

using namespace itensor;
using namespace std;

#include "check_task.h" // get_bool/get_str/get_int_value/... (get_sites.h needs these)
#include "get_sites.h" // SpinX, incl. the in-memory std::vector<int> constructor
#include "mo_terms.h" // OpFactor/MOTerm + build_ampo/build_mpo
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
        .def("excited_states",[](Chain& self, int n, double scale_lagrange,
                                  bool gram_schmidt) {
                auto out = self.excited_states(n,scale_lagrange,gram_schmidt);
                return py::make_tuple(out.energies,out.fluctuations,
                                      out.wavefunctions);
            }, py::arg("n"),py::arg("scale_lagrange")=1.0,
               py::arg("gram_schmidt")=false,
               "Returns (energies, fluctuations, wavefunctions)")
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
        .def("build_operator",[](Chain& self, std::vector<PyTerm> const& terms) {
                return self.build_operator(terms_from_python(terms));
            }, py::arg("terms"))
        .def("apply_pure_operator",&Chain::apply_pure_operator,
             py::arg("A"),py::arg("wf"))
        .def("multiply_operators",&Chain::multiply_operators,
             py::arg("A"),py::arg("B"))
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
