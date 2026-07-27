"""Regression coverage for the generalized-eigenvalue DMRG solver
H|psi>=lambda*A|psi> (A Hermitian positive definite): pyitensor's
dmrg_generalized() (src/dmrgpy/pyitensor/dmrg.py) and its ITensor v3
C++ port (Chain::gs_energy_generalized, mpscpp3/chain_session.h) -- see
either docstring for the self-consistent Lagrange-multiplier algorithm.
Both are exposed identically as Many_Body_Chain.gs_energy_generalized(),
so every accuracy test below is parametrized over itensor_version in
("python", 3) and runs against both implementations; the v3 case is
skipped automatically when the compiled extension isn't available
(mirrors how the rest of this suite handles a missing compiler -- see
_helpers.py's own versions= kwarg for the same idea applied to plain
gs_energy()).

Same test problem/ground-truth convention as test_arpacktk_iram.py's
mode-2 (ARPACK generalized eigenproblem) tests, so all three solvers are
directly comparable -- see also
examples/groundstate/dmrg_generalized_benchmark for a head-to-head
wall-time/accuracy comparison across all three.
"""
import numpy as np
import pytest
import scipy.linalg as sla

from dmrgpy import cppext, fermionchain
from dmrgpy.multioperator import identity as mo_identity

ITENSOR_VERSIONS = [
    "python",
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3), reason="ITensor v3 extension not compiled")),
]


def _generalized_fermion_problem(n=4, seed=2):
    """A*x = lambda*M*x: A a Hermitian interacting fermion Hamiltonian,
    M = 1 + 0.5*sum(N_i) a Hermitian positive-definite weight operator
    (diagonal, eigenvalues in [1, 1+0.5n], never zero)."""
    rng = np.random.RandomState(seed)
    fc = fermionchain.Fermionic_Chain(n)
    a = 0
    for i in range(n - 1):
        a = a + rng.random() * (fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i])
    for i in range(n - 1):
        a = a + rng.random() * (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    m = 1
    for i in range(n):
        m = m + 0.5 * fc.N[i]
    fc.set_hamiltonian(a)
    return fc, a, m


def _generalized_ground_truth(fc, a, m):
    Amat = fc.get_ED_obj().MO2matrix(a)
    Mmat = fc.get_ED_obj().MO2matrix(m)
    Amat = Amat.toarray() if hasattr(Amat, "toarray") else np.array(Amat)
    Mmat = Mmat.toarray() if hasattr(Mmat, "toarray") else np.array(Mmat)
    return np.sort(sla.eigh(Amat, Mmat, eigvals_only=True))


def _setup(fc, itensor_version, maxm=40, nsweeps=12, cutoff=1e-12):
    if itensor_version == "python":
        fc.setup_python()
    else:
        fc.setup_cpp(version=itensor_version)
    fc.maxm = maxm
    fc.nsweeps = nsweeps
    fc.cutoff = cutoff
    return fc


@pytest.mark.parametrize("itensor_version", ITENSOR_VERSIONS)
def test_gs_energy_generalized_matches_ed(itensor_version):
    """The generalized-eigenvalue DMRG ground state must agree with
    scipy.linalg.eigh's generalized Hermitian-definite eigensolver (the
    same ground truth test_arpacktk_iram.py's mode-2 tests use)."""
    fc, a, m = _generalized_fermion_problem(n=4, seed=2)
    w = _generalized_ground_truth(fc, a, m)
    _setup(fc, itensor_version)
    lam = fc.gs_energy_generalized(m)
    assert abs(lam - w[0]) < 1e-6


@pytest.mark.parametrize("itensor_version", ITENSOR_VERSIONS)
def test_gs_energy_generalized_matches_ed_larger_chain(itensor_version):
    fc, a, m = _generalized_fermion_problem(n=6, seed=2)
    w = _generalized_ground_truth(fc, a, m)
    _setup(fc, itensor_version, maxm=60, nsweeps=14)
    lam = fc.gs_energy_generalized(m)
    assert abs(lam - w[0]) < 1e-4


@pytest.mark.parametrize("itensor_version", ITENSOR_VERSIONS)
def test_gs_energy_generalized_identity_matches_plain_dmrg(itensor_version):
    """A = identity must reduce the generalized problem exactly to the
    chain's own plain ground-state energy <H> (dmrg_generalized()'s own
    docstring: A=Id makes lambda the ordinary Rayleigh quotient, i.e.
    plain DMRG)."""
    fc_plain, a_plain, _ = _generalized_fermion_problem(n=4, seed=3)
    _setup(fc_plain, itensor_version)
    e_plain = fc_plain.gs_energy()

    fc_gen, a_gen, _ = _generalized_fermion_problem(n=4, seed=3)
    _setup(fc_gen, itensor_version)
    lam = fc_gen.gs_energy_generalized(mo_identity())
    assert abs(lam - e_plain) < 1e-6


@pytest.mark.parametrize("itensor_version", ITENSOR_VERSIONS)
def test_gs_energy_generalized_lam0_hint_does_not_change_result(itensor_version):
    """A caller-supplied starting lambda estimate (lam0) is only meant to
    speed up convergence, not change the converged answer -- a poor
    lam0 (0.0, far from the true eigenvalue) must still converge to the
    same lambda as the default (data-driven) starting guess."""
    fc, a, m = _generalized_fermion_problem(n=4, seed=2)
    w = _generalized_ground_truth(fc, a, m)
    _setup(fc, itensor_version, nsweeps=16)
    lam = fc.gs_energy_generalized(m, lam0=0.0)
    assert abs(lam - w[0]) < 1e-4


def test_gs_energy_generalized_requires_supported_backend():
    """The generalized solver only exists on itensor_version 'python' and
    3 so far (see CLAUDE.md's "Implement in pyitensor first" scoping,
    since extended to v3) -- mpscpp2 (itensor_version=2) and julia_live
    don't have this session method yet, so this must fail loudly rather
    than silently doing the wrong thing on those backends."""
    fc, a, m = _generalized_fermion_problem(n=4, seed=2)
    fc.setup_cpp(version=2)
    with pytest.raises(NotImplementedError):
        fc.gs_energy_generalized(m)


def test_gs_energy_generalized_requires_compiled_v3_extension():
    """itensor_version=3 with no compiled extension for it (self._session
    left as None by sites.py's initialize(), the same state mode.py's own
    get_mode() falls back to ED for elsewhere) must raise a clear
    RuntimeError, not an AttributeError from calling into a None
    session."""
    fc, a, m = _generalized_fermion_problem(n=4, seed=2)
    fc.setup_cpp(version=3)
    fc._session = None  # simulate "extension not compiled" regardless of this env
    with pytest.raises(RuntimeError):
        fc.gs_energy_generalized(m)


@pytest.mark.parametrize("itensor_version", ITENSOR_VERSIONS)
def test_gs_energy_generalized_requires_hamiltonian(itensor_version):
    fc = fermionchain.Fermionic_Chain(4)
    if itensor_version == "python":
        fc.setup_python()
    else:
        fc.setup_cpp(version=itensor_version)
    with pytest.raises(RuntimeError):
        fc.gs_energy_generalized(mo_identity())
