"""One-site profiles and two-point entries in one sweep on
itensor_version="python" (pyitensor/chain.py's sweep_measure,
expect_sweep and correlation_sweep, the algorithm of ITensorMPS.jl's
expect and correlation_matrix), against the per-term MPO path they
replace, and the callers routed through them (correlator.sweep_vevs)
against their old per-term loops."""
from contextlib import contextmanager
from functools import lru_cache

import numpy as np
import pytest

from dmrgpy import bosonchain, correlator, fermionchain, kondocloud, spinchain

TOL = 1e-12


@contextmanager
def per_term():
    """The callers on their old path, one vev per operator."""
    old = correlator.ONE_SWEEP
    correlator.ONE_SWEEP = False
    try:
        yield
    finally:
        correlator.ONE_SWEEP = old


def _solve(chain, h, sector=None):
    chain.set_hamiltonian(h)
    if sector:
        chain.set_conserved_sector(**sector)
    chain.gs_energy()
    return chain


@lru_cache(maxsize=None)
def spins(kind="half"):
    """A complex state: XXZ with a Dzyaloshinskii-Moriya term, and fields
    that break Sz conservation unless a sector is asked for."""
    labels = {"half": ["S=1/2"]*8, "mixed": ["S=1/2", "S=1"]*3,
              "sector": ["S=1/2"]*8}[kind]
    sc = spinchain.Spin_Chain(labels, itensor_version="python", maxm=40,
                              nsweeps=8)
    n = sc.ns
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
        h = h + 0.7*sc.Sz[i]*sc.Sz[i+1]
        h = h + 0.4*(sc.Sx[i]*sc.Sy[i+1] - sc.Sy[i]*sc.Sx[i+1])
    if kind == "sector":
        h = h + 0.2*sc.Sz[0]
        return _solve(sc, h, sector=dict(Sz=0))
    h = h + 0.3*sc.Sx[0] + 0.2*sc.Sz[n-1] + 0.1*sc.Sy[2]
    return _solve(sc, h)


@lru_cache(maxsize=None)
def spinless(kind="plain"):
    """Spinless fermions with a complex hopping. "parity" adds a pairing
    term, so the state is not one of definite particle number."""
    fc = fermionchain.Fermionic_Chain(8, itensor_version="python", maxm=40,
                                      nsweeps=8)
    n = fc.ns
    h = 0
    for i in range(n-1):
        h = h + np.exp(0.7j)*fc.Cdag[i]*fc.C[i+1]
        if kind == "parity":
            h = h + 0.3*fc.C[i]*fc.C[i+1]
    h = h + h.get_dagger()
    h = h + sum(1.5*fc.N[i]*fc.N[i+1] for i in range(n-1)) + 0.2*fc.N[0]
    return _solve(fc, h, sector=dict(Nf=3) if kind == "sector" else None)


@lru_cache(maxsize=None)
def interleaved():
    fc = fermionchain.Spinful_Fermionic_Chain(4, itensor_version="python",
                                              maxm=40, nsweeps=8)
    h = 0
    for i in range(3):
        h = h + np.exp(0.3j)*(fc.Cdagup[i]*fc.Cup[i+1]
                              + fc.Cdagdn[i]*fc.Cdn[i+1])
    h = h + 0.2*fc.Cdagup[1]*fc.Cdn[1]  # spin mixing, so <Sx> is not 0
    h = h + h.get_dagger()
    h = h + 2.0*sum(fc.Nup[i]*fc.Ndn[i] for i in range(4)) + 0.3*fc.Nup[0]
    return _solve(fc, h)


@lru_cache(maxsize=None)
def native():
    fc = fermionchain.Spinful_Fermionic_Chain_Native(4, itensor_version="python",
                                                     maxm=40, nsweeps=8)
    h = 0
    for i in range(3):
        h = h + np.exp(0.3j)*(fc.Cdagup[i]*fc.Cup[i+1]
                              + fc.Cdagdn[i]*fc.Cdn[i+1])
    h = h + 0.2*fc.Cdagup[1]*fc.Cdn[1]
    h = h + h.get_dagger()
    h = h + 2.0*sum(fc.Nup[i]*fc.Ndn[i] for i in range(4)) + 0.3*fc.Nup[0]
    return _solve(fc, h)


@lru_cache(maxsize=None)
def bosons():
    bc = bosonchain.Bosonic_Chain(6, maxnb=[3]*6, itensor_version="python",
                                  maxm=40, nsweeps=8)
    n = bc.ns
    h = 0
    for i in range(n-1):
        h = h - np.exp(0.5j)*bc.Adag[i]*bc.A[i+1]
    h = h + 0.2*bc.A[2]  # no particle-number conservation
    h = h + h.get_dagger()
    h = h + sum(0.8*bc.N[i]*bc.N[i] - 0.9*bc.N[i] for i in range(n))
    return _solve(bc, h)


CHAINS = {
    "spins": lambda: spins("half"),
    "mixed-spins": lambda: spins("mixed"),
    "spin-sector": lambda: spins("sector"),
    "spinless": lambda: spinless("plain"),
    "spinless-parity": lambda: spinless("parity"),
    "spinless-sector": lambda: spinless("sector"),
    "interleaved": interleaved,
    "native": native,
    "bosons": bosons,
}

ONE_SITE = [
    ("spins", "Sz"), ("spins", "Sx"), ("spins", "Sy"), ("spins", "S+"),
    ("spins", ("Sx", "Sz")),
    ("mixed-spins", "Sx"), ("mixed-spins", "S-"),
    ("spin-sector", "Sz"),
    ("spinless", "N"), ("spinless", ("N", "N")), ("spinless", ("Cdag", "C")),
    ("spinless-parity", "N"),
    ("spinless-sector", "N"),
    ("interleaved", "N"),
    ("native", "Nup"), ("native", ("Cdagup", "Cdn")),
    ("bosons", "N"), ("bosons", "A"), ("bosons", ("N", "N")),
]

TWO_SITE = [
    ("spins", "Sz", "Sz"), ("spins", "Sx", "Sx"), ("spins", "S+", "S-"),
    ("spins", "S+", "S+"), ("spins", "Sx", "Sz"), ("spins", ("Sx", "Sz"), "Sy"),
    ("mixed-spins", "S+", "S-"), ("mixed-spins", "Sx", "Sz"),
    ("spin-sector", "S+", "S-"), ("spin-sector", "Sz", "Sz"),
    ("spinless", "Cdag", "C"), ("spinless", "C", "Cdag"), ("spinless", "C", "C"),
    ("spinless", "Cdag", "Cdag"), ("spinless", "N", "N"),
    ("spinless-parity", "C", "C"), ("spinless-parity", "Cdag", "C"),
    ("spinless-sector", "Cdag", "C"), ("spinless-sector", "N", "N"),
    ("interleaved", "Cdag", "C"), ("interleaved", "C", "Cdag"),
    ("interleaved", "N", "N"),
    ("native", "Cdagup", "Cup"), ("native", "Cup", "Cdagdn"),
    ("native", ("Cdagup", "Cdn"), ("Cdagdn", "Cup")), ("native", "Nup", "Ndn"),
    ("bosons", "A", "Adag"), ("bosons", "Adag", "A"), ("bosons", "A", "A"),
    ("bosons", "N", "N"),
]


def _state(case):
    chain = CHAINS[case]()
    return chain, chain._session, chain.get_gs().cpp_handle


def _ids(cases):
    return ["-".join("*".join(x) if isinstance(x, tuple) else x for x in c)
            for c in cases]


@pytest.mark.parametrize("case,name", ONE_SITE, ids=_ids(ONE_SITE))
def test_expect_sweep_is_the_per_site_loop(case, name):
    _, s, wf = _state(case)
    out = s.expect_sweep(wf, name)
    ref = s.expect_per_site(wf, name)
    assert np.max(np.abs(out - ref)) < TOL


@pytest.mark.parametrize("case,a,b", TWO_SITE, ids=_ids(TWO_SITE))
def test_correlation_sweep_is_the_per_pair_loop(case, a, b):
    _, s, wf = _state(case)
    out = s.correlation_sweep(wf, a, b)
    ref = s.correlation_per_pair(wf, a, b)
    assert np.max(np.abs(out - ref)) < TOL


def test_the_values_are_not_trivial():
    """The complex spin state has a complex S+ profile and a non-Hermitian
    S+ S+ matrix, so a wrong conjugate or a dropped lower triangle shows."""
    _, s, wf = _state("spins")
    assert np.max(np.abs(s.expect_sweep(wf, "S+").imag)) > 1e-4
    pp = s.correlation_sweep(wf, "S+", "S+")
    assert np.max(np.abs(pp - pp.conj().T)) > 1e-4
    _, s, wf = _state("spinless-parity")
    cc = s.correlation_sweep(wf, "C", "C")
    assert np.max(np.abs(cc)) > 1e-4
    assert np.max(np.abs(cc + cc.T)) < TOL  # C_i C_j = -C_j C_i


def test_selected_pairs_and_several_requests_in_one_sweep():
    _, s, wf = _state("spinless")
    pairs = [(3, 1), (0, 5), (2, 2), (5, 0), (7, 6), (3, 1)]
    out = s.correlation_sweep(wf, "C", "Cdag", pairs=pairs)
    ref = s.correlation_per_pair(wf, "C", "Cdag", pairs=pairs)
    assert np.max(np.abs(out - ref)) < TOL
    v1, v2 = s.sweep_measure(wf, ones=[("N", [4, 1, 4]), (("N", "N"), [0])],
                             twos=[("Cdag", "C", pairs), ("N", "N", pairs)])
    assert np.max(np.abs(v1[0] - s.expect_per_site(wf, "N", [4, 1, 4]))) < TOL
    assert np.max(np.abs(v1[1] - s.expect_per_site(wf, ("N", "N"), [0]))) < TOL
    assert np.max(np.abs(v2[0] - s.correlation_per_pair(wf, "Cdag", "C",
                                                        pairs))) < TOL
    assert np.max(np.abs(v2[1] - s.correlation_per_pair(wf, "N", "N",
                                                        pairs))) < TOL


def test_cdag_c_matrix_is_correlation_matrix_sweep():
    _, s, wf = _state("interleaved")
    ref = s.correlation_matrix_per_pair(wf)
    out = s.correlation_matrix_sweep(wf)
    assert np.max(np.abs(out - ref)) < TOL
    assert np.max(np.abs(out - s.correlation_sweep(wf, "Cdag", "C"))) < 1e-13


def test_unnormalized_state_scales_as_its_norm_squared():
    chain, s, _ = _state("spins")
    wf = (chain.get_gs()*2.0).cpp_handle
    assert np.max(np.abs(s.expect_sweep(wf, "Sx")
                         - s.expect_per_site(wf, "Sx"))) < TOL
    assert np.max(np.abs(s.correlation_sweep(wf, "Sx", "Sz")
                         - s.correlation_per_pair(wf, "Sx", "Sz"))) < TOL
    wf1 = chain.get_gs().cpp_handle
    assert np.max(np.abs(s.expect_sweep(wf, "Sz")
                         - 4*s.expect_sweep(wf1, "Sz"))) < TOL


def test_a_state_with_no_known_centre():
    _, s, wf = _state("spins")
    psi = wf.copy()
    psi.center = None
    assert np.max(np.abs(s.expect_sweep(psi, "Sx")
                         - s.expect_per_site(wf, "Sx"))) < TOL
    assert np.max(np.abs(s.correlation_sweep(psi, "S+", "S-")
                         - s.correlation_per_pair(wf, "S+", "S-"))) < TOL


def test_fermion_odd_operators_are_refused():
    _, s, wf = _state("spinless-parity")
    with pytest.raises(ValueError, match="odd"):
        s.expect_sweep(wf, "C")
    with pytest.raises(ValueError, match="mixes"):
        s.correlation_sweep(wf, "Cdag", "N")
    s.expect_per_site(wf, "C")  # the per-term path takes it


def test_an_operator_that_breaks_the_sector_still_raises():
    chain, s, wf = _state("spin-sector")
    with pytest.raises(ValueError, match="conserved-sector"):
        s.expect_sweep(wf, "Sx")
    with pytest.raises(ValueError, match="conserved-sector"):
        s.correlation_sweep(wf, "Sx", "Sx")
    with pytest.raises(ValueError, match="conserved-sector"):
        s.correlation_sweep(wf, "S+", "Sz", pairs=[(0, 3)])
    with pytest.raises(ValueError):
        chain.get_magnetization()  # as before: vev(Sx) raises
    fc, s, wf = _state("spinless-sector")
    with pytest.raises(ValueError, match="conserved-sector"):
        s.correlation_sweep(wf, "C", "C")


def test_sweep_runs_on_jax():
    pytest.importorskip("jax")
    from dmrgpy.pyitensor import backend
    _, s, wf = _state("spinless")
    ref1 = s.expect_per_site(wf, "N")
    ref2 = s.correlation_per_pair(wf, "C", "Cdag")
    try:
        backend.set_backend("jax")
        out1 = s.expect_sweep(wf, "N")
        out2 = s.correlation_sweep(wf, "C", "Cdag")
    finally:
        backend.set_backend("numpy")
    assert np.max(np.abs(np.asarray(out1) - ref1)) < 1e-10
    assert np.max(np.abs(np.asarray(out2) - ref2)) < 1e-10


# -- the callers, against their old per-term loops ------------------------

def _same(new, old):
    new, old = np.asarray(new), np.asarray(old)
    assert new.shape == old.shape
    assert np.max(np.abs(new - old)) < TOL


@pytest.mark.parametrize("case", ["spins", "mixed-spins"])
def test_get_magnetization(case, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # it writes MAGNETIZATION.OUT
    chain = CHAINS[case]()
    new = chain.get_magnetization()
    with per_term():
        old = chain.get_magnetization()
    _same(new, old)
    wf = chain.get_gs()*2.0  # vev normalizes, and so does the route
    new = chain.get_magnetization(wf=wf)
    with per_term():
        old = chain.get_magnetization(wf=wf)
    _same(new, old)


PAIRS = [(0, 3), (3, 0), (2, 2), (1, 5), (5, 1), (6, 7)]


@pytest.mark.parametrize("case,name", [
    ("spins", "SS"), ("spins", "ZZ"), ("spins", "XY"), ("spins", "+-"),
    ("spins", "++"), ("spins", "XX"), ("mixed-spins", "SS"),
    ("spin-sector", "SS"), ("spin-sector", "ZZ"), ("spin-sector", "+-"),
    ("spinless", "cdc"), ("spinless", "ccd"), ("spinless", "cc"),
    ("spinless", "density"), ("spinless-parity", "cc"),
    ("spinless-sector", "cdc"), ("interleaved", "cdc"),
])
def test_get_correlator(case, name):
    chain = CHAINS[case]()
    pairs = [p for p in PAIRS if max(p) < chain.ns]
    new = correlator.get_correlator(chain, pairs=pairs, name=name)
    with per_term():
        old = correlator.get_correlator(chain, pairs=pairs, name=name)
    _same(new, old)


@pytest.mark.parametrize("case", ["spins", "spin-sector", "spinless",
                                  "spinless-sector"])
def test_chain_get_correlator(case):
    """Many_Body_Chain.get_correlator (Sz Sz on a spin chain, Cdag C on a
    spinless one) takes the sweep too."""
    chain = CHAINS[case]()
    if type(chain) not in (spinchain.Spin_Chain, fermionchain.Fermionic_Chain):
        pytest.skip("get_correlator covers these two classes only")
    pairs = [p for p in PAIRS if max(p) < chain.ns]
    new = chain.get_correlator(pairs=pairs)
    with per_term():
        old = chain.get_correlator(pairs=pairs)
    _same(new, old)


@pytest.mark.parametrize("case", ["spinless", "spinless-sector",
                                  "spinless-parity", "interleaved", "bosons"])
def test_density_loops(case):
    chain = CHAINS[case]()
    calls = [chain.get_density, chain.get_density_fluctuation]
    if case == "interleaved":
        calls += [chain.get_density_spinless,
                  chain.get_density_fluctuation_spinless,
                  chain.get_magnetization]
    for call in calls:
        new = call()
        with per_term():
            old = call()
        _same(new, old)


def test_ed_and_other_keywords_keep_the_old_path():
    chain = CHAINS["spinless"]()
    _same(chain.get_density(mode="ED"),
          [chain.vev(chain.N[i], mode="ED").real for i in range(chain.ns)])
    assert correlator.sweep_vevs(chain, chain.N, mode="ED") is None
    assert correlator.sweep_vevs(chain, chain.N, npow=1) is None
    assert correlator.sweep_vevs(chain, chain.N) is not None


@pytest.mark.parametrize("chain_kind,connected,half", [
    ("native", False, False), ("native", True, False),
    ("native", False, True), ("interleaved", True, False),
    ("interleaved", False, False),
])
def test_get_kondo_cloud(chain_kind, connected, half):
    fc = kondocloud.anderson_lead(6, ed=-0.5, U=1.0, Gamma=0.5,
                                  chain=chain_kind, half_filling=half,
                                  itensor_version="python")
    fc.maxm = 40
    fc.nsweeps = 8
    fc.cutoff = 1e-12
    new = fc.get_kondo_cloud(connected=connected)
    with per_term():
        old = fc.get_kondo_cloud(connected=connected)
    for key in ("SdSi", "SdSd", "Sigma", "mu", "nd", "sum_rule"):
        _same(new[key], old[key])
