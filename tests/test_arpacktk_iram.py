"""Regression coverage for algebra/arpacktk.py's Implicitly Restarted
Arnoldi Method (IRAM, ported from ARPACK), the default MPS Arnoldi
solver used by get_excited_states (non-Hermitian, n>=2) and by
degeneracy.py's eigenvalue_degeneracy. Before this file existed, none of
that code path had any pytest coverage at all -- a real crash/silent-
wrong-answer bug in mpsiram_shift_invert would have shipped silently
through run_tests.py staying green (found and fixed via manual stress
testing, see the PR history).
"""
import numpy as np
import pytest

from dmrgpy import spinchain, fermionchain
from dmrgpy.algebra import arnolditk
from dmrgpy.algebra.arpacktk import mpsiram, mpsiram_shift_invert


def _heisenberg_chain(n=6, seed=1):
    rng = np.random.RandomState(seed)
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    h = 0
    for i in range(n - 1):
        h = h + rng.random() * sc.Sx[i] * sc.Sx[i + 1]
        h = h + rng.random() * sc.Sy[i] * sc.Sy[i + 1]
        h = h + rng.random() * sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    return sc, h


def _staggered_fermion_chain(n=6, seed=4):
    """Interacting fermion chain with a staggered imaginary potential --
    non-Hermitian, complex spectrum, exact complex-conjugate degenerate
    ground state (used elsewhere in this repo, e.g.
    examples/non_hermitian/nhdmrg_VS_ED_VS_arnoldi)."""
    rng = np.random.RandomState(seed)
    fc = fermionchain.Fermionic_Chain(n)
    mh = np.zeros((n, n), dtype=complex)
    for i in range(n - 1):
        mh[i, i + 1] = 1.0
        mh[i + 1, i] = 1.0
    for i in range(n):
        mh[i, i] = 1j * (-1) ** i
    h = 0
    for i in range(n):
        for j in range(n):
            h = h + mh[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1):
        h = h + (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    fc.set_hamiltonian(h)
    return fc, h


def test_mpsiram_hermitian_ground_state_ed():
    sc, h = _heisenberg_chain(n=6, seed=1)
    arnolditk.arnoldimode = "ED"
    es_ed = np.sort(np.array(sc.get_excited(mode="ED", n=1)).real)
    es, wfs = mpsiram(sc, h, which="SR", nev=1, tol=1e-8)
    assert abs(es[0].real - es_ed[0]) < 1e-6
    wf = wfs[0]
    assert abs(wf.dot(wf) - 1.0) < 1e-8 # normalized


def test_mpsiram_simultaneous_multi_eigenvalue_ed():
    """A genuine nev=3 simultaneous block search (not one-at-a-time)."""
    sc, h = _heisenberg_chain(n=6, seed=1)
    arnolditk.arnoldimode = "ED"
    es_ed = np.sort(np.array(sc.get_excited(mode="ED", n=3)).real)
    es, wfs = mpsiram(sc, h, which="SR", nev=3, tol=1e-6)
    assert len(es) == 3
    for e_iram, e_ed in zip(sorted(es.real), es_ed):
        assert abs(e_iram - e_ed) < 1e-4


def test_mpsiram_non_hermitian_conjugate_pair_ed():
    fc, h = _staggered_fermion_chain(n=6, seed=4)
    arnolditk.arnoldimode = "ED"
    es_ed = np.array(sorted(fc.get_excited(mode="ED", n=6), key=lambda x: x.real))
    es, wfs = mpsiram(fc, h, which="SR", nev=1, tol=1e-6)
    assert abs(es[0].real - es_ed[0].real) < 1e-4
    assert min(abs(es[0] - x) for x in es_ed) < 1e-4


def test_mpsiram_shift_invert_finds_nearest_eigenvalue_ed():
    sc, h = _heisenberg_chain(n=6, seed=1)
    arnolditk.arnoldimode = "ED"
    es_ed = np.sort(np.array(sc.get_excited(mode="ED", n=6)).real)
    target = es_ed[2] + 0.01
    es, wfs = mpsiram_shift_invert(sc, h, e=target, nev=1, delta=1e-6, maxn=40)
    closest = es_ed[np.argmin(np.abs(es_ed - target))]
    assert abs(es[0].real - closest) < 1e-3


def test_mpsiram_shift_invert_deflated_second_state_is_a_genuine_eigenvalue_ed():
    """A recursive (deflated) second shift-invert search must still land
    on *some* genuine eigenvalue of H, not a spurious combination --
    even though it is not guaranteed to land on the exact conjugate
    partner of the first (see module note on deflation and non-Hermitian
    operators: right eigenvectors of a non-Hermitian H need not be
    mutually orthogonal under the plain Hermitian inner product -- e.g.
    confirmed directly for this model's ground-state conjugate pair,
    |<v1|v2>| ~= 0.40 -- so .dot()-based deflation, used identically by
    arnolditk's own wfskip mechanism, is only a heuristic here, not a
    rigorous one)."""
    fc, h = _staggered_fermion_chain(n=6, seed=4)
    arnolditk.arnoldimode = "ED"
    es_ed = np.array(sorted(fc.get_excited(mode="ED", n=6), key=lambda x: x.real))
    e0 = es_ed[0].real
    es, wfs = mpsiram_shift_invert(fc, h, e=1.2 * e0, nev=2, delta=1e-3, maxn=40)
    for e in es:
        assert min(abs(e - x) for x in es_ed) < 1e-2


def test_excited_states_non_hermitian_default_matches_ed():
    """get_excited_states' non-Hermitian route (excited.py's
    excited_states_non_hermitian, IRAM by default) must agree with ED."""
    fc, h = _staggered_fermion_chain(n=6, seed=4)
    es_ed = np.array(sorted(fc.get_excited(mode="ED", n=3), key=lambda x: x.real))
    es, wfs = fc.get_excited_states(n=3, purify=False)
    es_sorted = np.array(sorted(es, key=lambda x: x.real))
    for e_iram, e_ed in zip(es_sorted, es_ed):
        assert abs(e_iram.real - e_ed.real) < 1e-3


def test_eigenvalue_degeneracy_non_degenerate_target_ed():
    """degeneracy.py's eigenvalue_degeneracy (shift-invert IRAM, called
    recursively with deflation) targeted at a genuinely non-degenerate
    eigenvalue must report a degeneracy close to 1.

    Note: this deliberately avoids this model's own exact
    complex-conjugate ground-state pair as a target -- right
    eigenvectors of a non-Hermitian H need not be mutually orthogonal
    under the plain Hermitian inner product (confirmed directly here,
    |<v1|v2>| ~= 0.40 for that pair), so the .dot()-based deflation this
    routine relies on (identical to arnolditk's own wfskip mechanism) is
    only a heuristic, not a rigorous one, for resolving non-orthogonal
    degenerate clusters -- a pre-existing characteristic of the
    deflation approach itself, not specific to this solver."""
    from dmrgpy.degeneracy import eigenvalue_degeneracy
    fc, h = _staggered_fermion_chain(n=4, seed=3)
    arnolditk.arnoldimode = "ED"
    es_ed = np.array(sorted(fc.get_excited(mode="ED", n=6), key=lambda x: x.real))
    # the "-0.25" level in this model is an isolated, non-degenerate
    # eigenvalue (confirmed via the ED spectrum)
    target = es_ed[np.argmin(np.abs(es_ed.real - (-0.25)))].real
    deg = eigenvalue_degeneracy(fc, h, target, n=1, delta=5e-2)
    assert 0.8 < deg < 1.5


@pytest.mark.parametrize("itensor_version", ["python"])
def test_mpsiram_dmrg_mode_matches_ed(itensor_version):
    """Real MPS/MPO (pyitensor backend, no compiled extension needed)
    ground state via IRAM must agree with ED."""
    sc, h = _heisenberg_chain(n=4, seed=1)
    sc.setup_python()
    arnolditk.arnoldimode = "DMRG"
    es_ed = np.sort(np.array(sc.get_excited(mode="ED", n=1)).real)
    es, wfs = mpsiram(sc, h, which="SR", nev=1, tol=1e-6)
    assert abs(es[0].real - es_ed[0]) < 1e-4
