"""Spin-operator algebra checks.

These exercise MultiOperator's +, *, get_dagger() and is_hermitian()/
is_zero_operator() on genuine spin operators, independent of the
Jordan-Wigner path (which is covered separately in
test_jordan_wigner.py).
"""
from dmrgpy import spinchain, fermionchain

from reference_data import (
    SPIN_COMMUTATOR_IS_ZERO,
    FERMION_HERMITIAN_H_IS_HERMITIAN,
)


def test_spin_commutator_is_zero():
    """[Sx, Sy] - i*Sz must vanish identically for a spin-1/2 site."""
    sc = spinchain.Spin_Chain(["S=1/2", "S=1/2", "S=1/2"])
    comm = sc.Sx[0] * sc.Sy[0] - sc.Sy[0] * sc.Sx[0] - 1j * sc.Sz[0]
    assert bool(sc.is_zero_operator(comm)) == SPIN_COMMUTATOR_IS_ZERO


def test_fermion_hamiltonian_is_hermitian():
    """A hopping Hamiltonian built as H + H.get_dagger() must be Hermitian.

    This exercises get_dagger() (rewritten to avoid its internal O(n^2))
    together with is_hermitian()'s sympy-based simplify() path, on a
    genuinely fermionic (not spin) operator.
    """
    L = 5
    fc = fermionchain.Fermionic_Chain(L)
    H = 0
    for i in range(L - 1):
        H = H + fc.Cdag[i] * fc.C[i + 1]
    H = H + H.get_dagger()
    assert bool(H.is_hermitian()) == FERMION_HERMITIAN_H_IS_HERMITIAN
