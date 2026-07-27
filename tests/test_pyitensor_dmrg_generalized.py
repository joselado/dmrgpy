"""Regression coverage for pyitensor/dmrg.py's dmrg_generalized(): a
two-site DMRG solver for the generalized eigenproblem
H|psi>=lambda*A|psi> (A Hermitian positive definite), exposed as
Chain.gs_energy_generalized() (pyitensor/chain.py) and
Many_Body_Chain.gs_energy_generalized() (groundstate.py/
manybodychain.py). See dmrg_generalized()'s own docstring for the
self-consistent Lagrange-multiplier algorithm: each outer sweep solves
the ordinary eigenproblem of H-lambda*A, then updates lambda to the
freshly-swept state's generalized Rayleigh quotient
<psi|H|psi>/<psi|A|psi>.

Same test problem/ground-truth convention as
test_arpacktk_iram.py's mode-2 (ARPACK generalized eigenproblem) tests,
so the two solvers are directly comparable -- see also
examples/groundstate/dmrg_generalized_benchmark for a head-to-head
wall-time/accuracy comparison.
"""
import numpy as np
import pytest
import scipy.linalg as sla

from dmrgpy import fermionchain
from dmrgpy.multioperator import identity as mo_identity


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


def _setup(fc, maxm=40, nsweeps=12, cutoff=1e-12):
    fc.setup_python()
    fc.maxm = maxm
    fc.nsweeps = nsweeps
    fc.cutoff = cutoff
    return fc


def test_gs_energy_generalized_matches_ed():
    """The pyitensor generalized-eigenvalue DMRG ground state must agree
    with scipy.linalg.eigh's generalized Hermitian-definite eigensolver
    (the same ground truth test_arpacktk_iram.py's mode-2 tests use)."""
    fc, a, m = _generalized_fermion_problem(n=4, seed=2)
    w = _generalized_ground_truth(fc, a, m)
    _setup(fc)
    lam = fc.gs_energy_generalized(m)
    assert abs(lam - w[0]) < 1e-6


def test_gs_energy_generalized_matches_ed_larger_chain():
    fc, a, m = _generalized_fermion_problem(n=6, seed=2)
    w = _generalized_ground_truth(fc, a, m)
    _setup(fc, maxm=60, nsweeps=14)
    lam = fc.gs_energy_generalized(m)
    assert abs(lam - w[0]) < 1e-4


def test_gs_energy_generalized_identity_matches_plain_dmrg():
    """A = identity must reduce the generalized problem exactly to the
    chain's own plain ground-state energy <H> (dmrg_generalized()'s own
    docstring: A=Id makes lambda the ordinary Rayleigh quotient, i.e.
    plain DMRG)."""
    fc_plain, a_plain, _ = _generalized_fermion_problem(n=4, seed=3)
    _setup(fc_plain)
    e_plain = fc_plain.gs_energy()

    fc_gen, a_gen, _ = _generalized_fermion_problem(n=4, seed=3)
    _setup(fc_gen)
    lam = fc_gen.gs_energy_generalized(mo_identity())
    assert abs(lam - e_plain) < 1e-6


def test_gs_energy_generalized_lam0_hint_does_not_change_result():
    """A caller-supplied starting lambda estimate (lam0) is only meant to
    speed up convergence, not change the converged answer -- a poor
    lam0 (0.0, far from the true eigenvalue) must still converge to the
    same lambda as the default (data-driven) starting guess."""
    fc, a, m = _generalized_fermion_problem(n=4, seed=2)
    w = _generalized_ground_truth(fc, a, m)
    fc.setup_python()
    fc.maxm, fc.nsweeps, fc.cutoff = 40, 16, 1e-12
    lam = fc.gs_energy_generalized(m, lam0=0.0)
    assert abs(lam - w[0]) < 1e-4


def test_gs_energy_generalized_requires_python_backend():
    """The generalized solver only exists in pyitensor so far -- see
    CLAUDE.md's "Implement in pyitensor first" scoping; mpscpp2/mpscpp3
    sessions don't have Chain.gs_energy_generalized (or its analogue)
    yet, so this must fail loudly rather than silently doing the wrong
    thing on those backends."""
    fc, a, m = _generalized_fermion_problem(n=4, seed=2)
    with pytest.raises(NotImplementedError):
        fc.gs_energy_generalized(m)


def test_gs_energy_generalized_requires_hamiltonian():
    fc = fermionchain.Fermionic_Chain(4)
    fc.setup_python()
    with pytest.raises(RuntimeError):
        fc.gs_energy_generalized(mo_identity())
