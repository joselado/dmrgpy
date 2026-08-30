"""The catalogue of calculations swept by run_benchmarks.py.

Kept separate from run_benchmarks.py (CLI + orchestration) and
reporttk.py (LaTeX/plots) so that "which capabilities does the benchmark
cover" is one readable list rather than something to reconstruct from the
sweep loop.

Two families of calculation live here:

* the **finite** family, swept over chain length N on a uniform S=1/2
  Heisenberg chain -- ground state, static correlator, KPM and TD
  dynamical correlators, excited states, entanglement entropy, real-time
  evolution after a quench, and the same ground state confined to a
  conserved Sz sector;
* the **infinite** family, swept over bond dimension chi on the same
  model in the thermodynamic limit (infinitechain.py) -- energy density
  and <Sz> under both gs_method="idmrg" and gs_method="vumps".

Every calculation is described by one CALCS entry, which says which
backends can run it, how to prepare it (untimed) and what to time, and
what to compare the answer against. Adding a capability to the benchmark
means adding one entry here plus a label in reporttk.CALC_LABELS.

Accuracy references come from exact diagonalization for the finite
family and from closed-form results for the infinite one (the
Bethe-ansatz energy density, and <Sz>=0 by symmetry), so a backend that
is fast because it is computing the wrong thing shows up as such.
"""
import numpy as np

from dmrgpy import spinchain
from dmrgpy import infinitechain
from dmrgpy import timedependent


# Bethe-ansatz ground-state energy per site of the uniform S=1/2
# Heisenberg chain in the thermodynamic limit (same constant used by
# examples/idmrg/heisenberg_infinite_python_VS_v3).
EXACT_DENSITY = 0.25 - np.log(2.0)


# ---------------------------------------------------------------------
# Problem construction
# ---------------------------------------------------------------------

def make_chain(n, itensor_version, params):
    """A uniform S=1/2 Heisenberg chain of length n on the requested
    backend -- the standard test Hamiltonian used across this codebase's
    own examples/tests (see e.g. test_benchmarks.py's Heisenberg test)."""
    spins = ["1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins, itensor_version=itensor_version)

    # spinchain.set_exchange() was removed; this is exactly the
    # Hamiltonian it built for this fj -- it summed the coupling over
    # BOTH orderings of every pair, so this is 2*sum_<ij> S_i.S_j. Kept
    # in that form so benchmark timings and energies stay comparable with
    # every run recorded before the removal.
    h = 0
    for i in range(n):
        for j in range(n):
            if abs(i - j) == 1: h = h + sc.SS(i, j)
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = params["maxm"], params["nsweeps"], params["cutoff"]
    return sc


def make_infinite_chain(chi, itensor_version, gs_method, params):
    """The same Heisenberg model directly in the thermodynamic limit, on
    a one-site unit cell (infinitechain.py). `chi` is the bond dimension,
    which is this family's sweep axis -- there is no chain length to
    sweep."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = gs_method
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0]
                       + ic.SzC[0] * ic.SzR[0])
    ic.maxm = chi
    ic.maxiter, ic.etol = params["inf_maxiter"], params["inf_etol"]
    return ic


def quench_state(sc, n, mode="DMRG"):
    """Ground state of a staggered Sz field (a Neel product state, unique
    and exactly representable, so DMRG and ED agree on it), with the
    chain's Heisenberg Hamiltonian restored afterwards -- the initial
    condition for the real-time evolution benchmark, following
    examples/time_evolution/tdvp_VS_ED_time_evolution."""
    h1 = sc.hamiltonian.copy()
    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs(mode=mode)
    sc.set_hamiltonian(h1)
    return wf


# ---------------------------------------------------------------------
# Finite-chain calculations
#
# Every prep_* takes the shared, already-converged chain for this
# (backend, N) plus the repeat index, and returns the *thunk to time*.
# It is called once per repeat, so anything a previous repeat cached
# cannot leak into the next one's timing: the DMRG session caches its
# ground-state energy (groundstate.py's send-cache), so a second
# gs_energy() on the same chain returns instantly and would otherwise be
# recorded as this backend's best time.
# ---------------------------------------------------------------------

def prep_gs(sc, iv, n, params, rep):
    # rep 0 uses the shared chain, whose ground state has deliberately
    # not been computed yet, so this timing doubles as the preparation
    # every other calculation below depends on. Later repeats need a
    # fresh chain (see the caching note above).
    target = sc if rep == 0 else make_chain(n, iv, params)
    return lambda: target.gs_energy(mode="DMRG")


def prep_static(sc, iv, n, params, rep):
    # <Sz_0 Sz_i> for every site, on top of the already-converged ground
    # state, so this times the correlator alone and not the DMRG sweep.
    return lambda: np.array([sc.vev(sc.Sz[0] * sc.Sz[i]).real for i in range(n)])


def prep_dynamic(sc, iv, n, params, rep):
    return lambda: sc.get_dynamical_correlator(
            mode="DMRG", submode="KPM", name=(sc.Sz[0], sc.Sz[0]),
            n=params["kpm_n"], es=params["es"])[1].real


def prep_dynamic_td(sc, iv, n, params, rep):
    # The same spectral function through a completely different engine:
    # real-time evolution plus a Fourier transform (timedependent.py),
    # i.e. TDVP on itensor_version 3/"python" and the MPO-Taylor
    # propagator on 2 -- worth timing next to KPM precisely because the
    # two submodes have different cost structures.
    return lambda: sc.get_dynamical_correlator(
            mode="DMRG", submode="TD", name=(sc.Sz[0], sc.Sz[0]),
            es=params["es"], delta=params["td_delta"])[1].real


def prep_excited(sc, iv, n, params, rep):
    # Energy of the first excited state via the overlap-penalty method
    # (excited.py); n=2 states are computed, the ground state among them.
    return lambda: float(np.real(sc.get_excited(n=2, mode="DMRG")[1]))


def prep_entropy(sc, iv, n, params, rep):
    wf = sc.get_gs()
    return lambda: np.array([sc.get_bond_entropy(wf, b - 1, b)
                             for b in range(1, n)])


def prep_evolution(sc, iv, n, params, rep):
    # A fresh chain: this benchmark quenches the Hamiltonian, and doing
    # that on the shared chain would reset the ground state every other
    # calculation here is built on.
    own = make_chain(n, iv, params)
    wf = quench_state(own, n)
    nt, dt = params["nt"], params["dt"]
    return lambda: timedependent.evolve_and_measure(
            own, operator=own.Sz[0], nt=nt, dt=dt, wf=wf)[1].real


def prep_sector(sc, iv, n, params, rep):
    # Conserved-quantum-number mode (set_conserved_sector), available on
    # itensor_version=3 and "python" only. Sz=0 holds the Heisenberg
    # ground state, so this must reproduce the *same* energy as the dense
    # run above -- the interesting number is the time.
    own = make_chain(n, iv, params)
    own.set_conserved_sector(Sz=0)
    return lambda: own.gs_energy(mode="DMRG")


# ---------------------------------------------------------------------
# Infinite-chain calculations (swept over bond dimension, not length)
# ---------------------------------------------------------------------

def _prep_inf_gs(ic, chi, iv, gs_method, params, rep):
    target = ic if rep == 0 else make_infinite_chain(chi, iv, gs_method, params)
    return lambda: target.gs_energy()


def _prep_inf_vev(ic, chi, iv, gs_method, params, rep):
    # <Sz> on the converged infinite state. Cheap next to the growth
    # loop, but not free: it builds the transfer-matrix fixed points, and
    # it is the observable that catches a wrongly gauged unit cell (an
    # energy-only check cannot -- see CLAUDE.md's iDMRG notes).
    return lambda: float(np.real(ic.vev("Sz", 0)))


def prep_idmrg_gs(ic, chi, iv, params, rep):
    return _prep_inf_gs(ic, chi, iv, "idmrg", params, rep)


def prep_idmrg_vev(ic, chi, iv, params, rep):
    return _prep_inf_vev(ic, chi, iv, "idmrg", params, rep)


def prep_vumps_gs(ic, chi, iv, params, rep):
    return _prep_inf_gs(ic, chi, iv, "vumps", params, rep)


def prep_vumps_vev(ic, chi, iv, params, rep):
    return _prep_inf_vev(ic, chi, iv, "vumps", params, rep)


# ---------------------------------------------------------------------
# The catalogue
#
#   key        row/plot identifier, also the --calcs name
#   group      "finite" (swept over N) or "infinite" (swept over chi)
#   method     infinite family only: which gs_method the entry uses
#   backends   None = every available backend; otherwise a whitelist
#   ref        key into the reference dict (None = timing only)
#   ref_kind   "scalar" or "array", how the accuracy column compares
#
# The prose caveats that go with some of these (the evolution column
# timing a different algorithm on itensor_version=2, the sector column
# existing on two backends only, ...) live in reporttk.CALC_NOTES, so
# that no LaTeX ends up in this file.
# ---------------------------------------------------------------------

CALCS = [
    dict(key="gs", group="finite", prep=prep_gs, backends=None,
         ref="gs", ref_kind="scalar"),
    dict(key="static", group="finite", prep=prep_static, backends=None,
         ref="static", ref_kind="array"),
    dict(key="dynamic", group="finite", prep=prep_dynamic, backends=None,
         ref=None, ref_kind=None),
    dict(key="dynamic_td", group="finite", prep=prep_dynamic_td, backends=None,
         ref=None, ref_kind=None),
    dict(key="excited", group="finite", prep=prep_excited, backends=None,
         ref="excited", ref_kind="scalar"),
    dict(key="entropy", group="finite", prep=prep_entropy, backends=None,
         ref="entropy", ref_kind="array"),
    dict(key="evolution", group="finite", prep=prep_evolution, backends=None,
         ref="evolution", ref_kind="array"),
    dict(key="sector", group="finite", prep=prep_sector,
         backends=("v3", "python"), ref="gs", ref_kind="scalar"),
    dict(key="idmrg_gs", group="infinite", method="idmrg", prep=prep_idmrg_gs,
         backends=("v3", "python"), ref="density", ref_kind="scalar"),
    dict(key="idmrg_vev", group="infinite", method="idmrg", prep=prep_idmrg_vev,
         backends=("v3", "python"), ref="zero", ref_kind="scalar"),
    dict(key="vumps_gs", group="infinite", method="vumps", prep=prep_vumps_gs,
         backends=("v3", "python"), ref="density", ref_kind="scalar"),
    dict(key="vumps_vev", group="infinite", method="vumps", prep=prep_vumps_vev,
         backends=("v3", "python"), ref="zero", ref_kind="scalar"),
]

CALC_KEYS = [c["key"] for c in CALCS]
FINITE_KEYS = [c["key"] for c in CALCS if c["group"] == "finite"]
INFINITE_KEYS = [c["key"] for c in CALCS if c["group"] == "infinite"]


def spec(key):
    for c in CALCS:
        if c["key"] == key:
            return c
    raise KeyError(key)


def allows(key, backend_key):
    """Whether this calculation exists on that backend at all -- the
    conserved sector and the infinite chain are implemented only on
    itensor_version=3 and "python", and a blank column there is by
    construction, not a failure."""
    allowed = spec(key)["backends"]
    return allowed is None or backend_key in allowed


# ---------------------------------------------------------------------
# Reference values
# ---------------------------------------------------------------------

def _entropy_from_vector(v, n, d=2):
    """Bond entanglement entropies of a full ED state vector, one per
    bond, by SVD of the reshaped amplitude tensor -- the ED counterpart
    of entropy.py's compute_entropy_single(), which is
    self._session-only and so has no ED path of its own."""
    out = []
    for b in range(1, n):
        s = np.linalg.svd(np.asarray(v).reshape(d ** b, -1), compute_uv=False)
        p = s ** 2
        p = p[p > 1e-16]
        out.append(float(-np.sum(p * np.log(p))))
    return np.array(out)


def finite_references(n, params, calc_keys):
    """Exact-diagonalization reference values for the finite-family
    calculations in `calc_keys`, at chain length n.

    Computed on the pure-Python backend, since ED needs no compiled
    extension to even construct the chain and the ED path itself is
    identical for every itensor_version. Returns {calc_ref_key: value}
    plus an "ok"/"err" status; a reference that raises leaves that
    calculation's accuracy column blank rather than aborting the sweep."""
    wanted = set()
    for key in calc_keys:
        s = spec(key)
        if s["group"] == "finite" and s["ref"] is not None:
            wanted.add(s["ref"])
    out = dict(n=n, ok=True, err=None)
    if not wanted:
        return out
    try:
        sc = make_chain(n, "python", params)
        if "gs" in wanted:
            out["gs"] = sc.gs_energy(mode="ED")
        if "static" in wanted:
            out["static"] = np.array([sc.vev(sc.Sz[0] * sc.Sz[i], mode="ED").real
                                      for i in range(n)])
        if "excited" in wanted:
            out["excited"] = float(np.real(sc.get_excited(n=2, mode="ED")[1]))
        if "entropy" in wanted:
            out["entropy"] = _entropy_from_vector(
                    np.asarray(sc.get_ED_obj().get_gs_array()).flatten(), n)
        if "evolution" in wanted:
            ed = make_chain(n, "python", params)
            wf = quench_state(ed, n, mode="ED")
            out["evolution"] = timedependent.evolve_and_measure(
                    ed, operator=ed.Sz[0], nt=params["nt"], dt=params["dt"],
                    wf=wf, mode="ED")[1].real
    except Exception as e:
        out["ok"], out["err"] = False, repr(e)
    return out


def infinite_reference(ref_kind_key):
    """Closed-form reference for an infinite-chain calculation: the
    Bethe-ansatz energy density, or zero for <Sz>."""
    if ref_kind_key == "density":
        return EXACT_DENSITY
    if ref_kind_key == "zero":
        return 0.0
    return None
