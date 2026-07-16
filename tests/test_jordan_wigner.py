"""Jordan-Wigner regression tests.

multioperator.jordan_wigner() and get_dagger() were rewritten (to
accumulate terms into a flat list instead of the public "+" operator,
which made them accidentally O(n^2)) without intending to change their
output. These tests pin that output against golden values computed
with the pre-optimization implementation (see reference_data.py).
"""
from dmrgpy import fermionchain

from reference_data import JW_TERMS_4FERMION, JW_TERMS_CDAGC, DAGGER_TERMS_CDAG1_C4
from _helpers import normalize_terms


def test_jordan_wigner_two_fermion_terms():
    """Cdag_i * C_j at various distances and orderings (i==j, adjacent,
    far apart, i<j and i>j) must Jordan-Wigner-transform to the same
    Fermi-string of F operators as the original implementation."""
    fc = fermionchain.Fermionic_Chain(6)
    cases = [(0, 0), (0, 1), (1, 0), (0, 3), (3, 0), (2, 5)]
    for (i, j) in cases:
        op = fc.Cdag[i] * fc.C[j]
        got = normalize_terms(op.to_terms())
        expected = JW_TERMS_CDAGC["%d_%d" % (i, j)]
        assert got == expected, "Cdag_%d*C_%d: got %r, expected %r" % (i, j, got, expected)


def test_jordan_wigner_four_fermion_term():
    """A 4-fermion product (Cdag_0 C_2 Cdag_3 C_1) must Jordan-Wigner-
    transform identically to the original implementation."""
    fc = fermionchain.Fermionic_Chain(6)
    op = fc.Cdag[0] * fc.C[2] * fc.Cdag[3] * fc.C[1]
    got = normalize_terms(op.to_terms())
    assert got == JW_TERMS_4FERMION


def test_get_dagger_matches_reference():
    """(Cdag_1 * C_4)^dagger, after Jordan-Wigner, must match the
    original implementation's output (checks both the operator-name
    dagger map and the reversed term order)."""
    fc = fermionchain.Fermionic_Chain(6)
    op = fc.Cdag[1] * fc.C[4]
    opd = op.get_dagger()
    got = normalize_terms(opd.to_terms())
    assert got == DAGGER_TERMS_CDAG1_C4
