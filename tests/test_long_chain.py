"""Functional coverage at a much larger scale (L~30) than the rest of
this test suite: exact diagonalization is infeasible here (2^30 states
for a spin-1/2/spinless-fermion chain, far more for the spinful/Hubbard
case), so these tests only exercise the compiled C++ DMRG backend(s)
directly, checking ground-state energy and a static correlator against
golden regression values -- and are skipped if no backend is compiled.

These golden values were cross-checked two ways before being hardcoded:
1. Against ITensor v2 and v3 run independently on the current tree (see
   below) -- they agree to ~1e-8 or tighter for all three systems.
2. Against a real compiled DMRG run of the snapshot from ~2 days ago
   (commit 2701a07, before the v2/v3 split and the in-process pybind11
   rewrite): same physics, same energies/correlators to ~1e-4-1e-6,
   confirming the intervening MultiOperator/backend rewrite didn't
   change any results at this scale either.

Each test computes the ground state once and reuses it (wf.dot(op*wf))
for every correlator point, instead of the "set_hamiltonian + vev per
point" pattern used elsewhere in this suite: at this system size,
re-running the full DMRG minimization per observable would make these
tests take minutes instead of seconds.
"""
import pytest

from dmrgpy import spinchain, fermionchain
from dmrgpy import cppext

pytestmark = pytest.mark.skipif(
    not (cppext.available(2) or cppext.available(3)),
    reason="long-chain tests need a compiled C++ DMRG backend (ED is infeasible at this size)",
)

TOL = 1e-3


def _available_versions():
    return [v for v in (2, 3) if cppext.available(v)]


def test_long_spin_chain_energy_and_correlator():
    """30-site Heisenberg chain with a small uniform field (to avoid
    the exact Sz degeneracies a plain Heisenberg chain would have)."""
    n = 30
    spins = ["S=1/2" for _ in range(n)]

    expected_e = -13.111355744816176
    expected_corr = {
        0: 0.24999999999805514,
        5: -0.0417185121466323,
        10: 0.014746630276378143,
        15: -0.01179528244287141,
        20: 0.006401055911225504,
        25: -0.006039454604966251,
    }

    for version in _available_versions():
        sc = spinchain.Spin_Chain(spins)
        h = 0
        for i in range(n - 1):
            h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
        for i in range(n):
            h = h + 0.05 * sc.Sz[i]
        sc.set_hamiltonian(h)
        sc.setup_cpp(version)
        sc.maxm = 60
        sc.nsweeps = 20

        wf = sc.get_gs(mode="DMRG")
        e = sc.gs_energy(mode="DMRG")
        assert e == pytest.approx(expected_e, abs=TOL)

        for j, val in expected_corr.items():
            c = wf.dot(sc.Sz[0] * sc.Sz[j] * wf).real
            assert c == pytest.approx(val, abs=TOL)


def test_long_spinless_fermion_chain_energy_and_correlator():
    """30-site spinless fermion chain: hopping, nearest-neighbor
    repulsion, and a chemical potential."""
    n = 30

    expected_e = -24.26702413637703
    expected_corr = {
        0: 0.5902465297213813,
        5: -0.1018613192405062,
        10: -0.0023644321458177898,
        15: 0.037545877217042875,
        20: -0.0024941370885131716,
        25: -0.025144759501689963,
    }

    for version in _available_versions():
        fc = fermionchain.Fermionic_Chain(n)
        h = 0
        for i in range(n - 1):
            h = h + fc.Cdag[i] * fc.C[i + 1]
        h = h + h.get_dagger()
        V = 0.5
        for i in range(n - 1):
            h = h + V * fc.N[i] * fc.N[i + 1]
        mu = -0.5
        for i in range(n):
            h = h + mu * fc.N[i]
        fc.set_hamiltonian(h)
        fc.setup_cpp(version)
        fc.maxm = 60
        fc.nsweeps = 20

        wf = fc.get_gs(mode="DMRG")
        e = fc.gs_energy(mode="DMRG")
        assert e == pytest.approx(expected_e, abs=TOL)

        for j, val in expected_corr.items():
            c = wf.dot(fc.Cdag[0] * fc.C[j] * wf).real
            assert c == pytest.approx(val, abs=TOL)


def test_long_spinful_fermion_chain_energy_and_correlator():
    """15-orbital (30 physical site) Hubbard chain: hopping plus a
    particle-hole symmetric on-site U. Like
    test_spinful_fermion_chain.py, this includes a small explicit
    Zeeman-like term to lift the exact up/down degeneracy the plain
    Hamiltonian would otherwise have -- without it, an unconstrained
    DMRG search can converge to either degenerate ground state, making
    <Nup_0>-type correlators ill-defined (see that file's docstring)."""
    n = 15

    expected_e = -19.680735817134664
    expected_corr = {
        0: 0.22180972667943885,
        3: 0.09886975583018276,
        6: 0.04248518281217159,
        9: -0.01588329449500078,
        12: -0.022249370620915483,
    }

    for version in _available_versions():
        fc = fermionchain.Spinful_Fermionic_Chain(n)
        h = 0
        for i in range(n - 1):
            h = h + fc.Cdagup[i] * fc.Cup[i + 1]
            h = h + fc.Cdagdn[i] * fc.Cdn[i + 1]
        U = 1.0
        for i in range(n):
            h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
        h = h + h.get_dagger()
        eps = 0.3
        h = h + eps * (fc.Nup[0] - fc.Ndn[0])
        fc.set_hamiltonian(h)
        fc.setup_cpp(version)
        fc.maxm = 100
        fc.nsweeps = 25
        fc.noise = 5e-2

        wf = fc.get_gs(mode="DMRG")
        e = fc.gs_energy(mode="DMRG")
        assert e == pytest.approx(expected_e, abs=TOL)

        for j, val in expected_corr.items():
            c = wf.dot(fc.Cdagup[0] * fc.Cup[j] * wf).real
            assert c == pytest.approx(val, abs=TOL)
