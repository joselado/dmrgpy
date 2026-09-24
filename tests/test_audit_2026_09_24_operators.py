"""Regression tests for the operators cluster of the 2026-09-24 audit.

Each test here locks in finding 1 or finding 2 of
`docs/audit_2026_09_24_hole_hunt.md`, which records the original symptom,
the reproduction that was executed and the reviewer's analysis. Both are
about what the symbolic layer believes an operator's adjoint or exchange
sign to be:

- #1: `canonical.py` graded the bare post-Jordan-Wigner names `A`/`Adag`
  (and the spinful ones) as even against the pre-transform `C`/`Cdag`,
  but with C_j = F_0...F_{j-1} A_j the exchange sign of A_i and C_j
  depends on which one sits on the lower site. A term mixing the two is
  now left as written, so the proofs refuse it instead of proving an
  anti-Hermitian operator Hermitian or an operator of norm 4 zero.
- #2: `get_dagger()` passed the anti-Hermitian `ISy` through unchanged,
  so every correlator daggering its first operator came out exactly
  minus itself. The dagger now carries a phase, ISy^dagger = -ISy.

Chains are 3 or 4 sites and run on `itensor_version="python"`, which is
always available.
"""

import numpy as np
import pytest

from dmrgpy import fermionchain, bosonchain, spinchain, multioperator
from dmrgpy.multioperatortk import canonical


# ---------------------------------------------------------------- helpers

def fermion_chain(n=4):
    """Spinless first-neighbor hopping chain on the pure-Python backend."""
    fc = fermionchain.Fermionic_Chain(n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 16, 10
    return fc


def spin_chain_in_a_field(n=4, hx=0.35, hz=0.15):
    """The record's finding-2 chain: S=1/2 Heisenberg plus a tilted field,
    so that <GS|ISy_0 ... Sz_1|GS> has nonzero Lehmann weights."""
    sc = spinchain.Spin_Chain([2] * n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
              + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n): h = h + hx * sc.Sx[i] + hz * sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 12
    return sc


def terms(MO):
    """The term list as comparable (coefficient, factors) pairs."""
    return [(complex(t[0]), tuple((o[0], o[1]) for o in t[1:])) for t in MO.op]


# ------------------------------------- #1 C-type and A-type names mixed

def test_mixed_anti_hermitian_operator_is_not_proven_hermitian():
    """P = Adag_0 C_1 + A_0 Cdag_1 is exactly anti-Hermitian, ||P-P^dag||
    = 5.66 on dense matrices, and the canonical form used to prove it
    Hermitian. The chain-level check, which probes numerically when the
    proof does not land, has to say the same."""
    fc = fermion_chain()
    P = fc.Adag[0] * fc.C[1] + fc.A[0] * fc.Cdag[1]
    assert not P.is_hermitian()
    assert not fc.is_hermitian(P)


def test_mixed_commutator_of_norm_four_is_not_proven_zero():
    """X = C_1 A_0 - A_0 C_1 has Frobenius norm 4 (A_0 anticommutes with
    C_1, whose string runs through site 0), and was proven zero, which
    also zeroed operator_norm() and is_zero_operator()."""
    fc = fermion_chain()
    X = fc.C[1] * fc.A[0] - fc.A[0] * fc.C[1]
    assert not X.is_zero()
    assert not fc.is_zero_operator(X, ntries=3)


def test_simplify_keeps_the_value_of_a_mixed_term():
    """simplify() used to flip the sign of Cdag_3 Adag_0, one sign per
    C-type factor written left of an A-type factor on a lower site."""
    fc = fermion_chain()
    np.random.seed(3)
    w1 = fc.random_mps()
    w2 = fc.random_mps()
    T = fc.Cdag[3] * fc.Adag[0]
    literal = fc.aMb(w1, T, w2)
    canonical_value = fc.aMb(w1, T.simplify(), w2)
    assert abs(literal) > 1e-2 # a sign flip would be visible
    assert canonical_value == pytest.approx(literal, abs=1e-10)


def test_td_pair_mixing_representations_is_not_a_dagger_pair():
    """The record's TD pair: A and B^dagger differ by the sign of one term
    (Adag_0 Cdag_1 = -Cdag_1 Adag_0), and proving them a dagger pair sent
    TD down the one-evolution Re F shortcut, 81 per cent off the peak.
    A genuine dagger of the same mixed B is still recognized, since
    identical spellings collect whatever their names."""
    fc = fermion_chain()
    B = fc.C[2] * fc.C[3] + fc.A[0] * fc.C[1]
    A = fc.Cdag[3] * fc.Cdag[2] + fc.Adag[0] * fc.Cdag[1]
    assert not canonical.is_dagger_pair(A, B)
    assert canonical.is_dagger_pair(B.get_dagger(), B)


def test_single_representation_proofs_are_kept():
    """The refusal is narrow: giving the A-type names no parity at all
    would also have refused the boson hopping and every
    Jordan-Wigner-transformed operator, which are right as they are."""
    fc = fermion_chain()
    C, Cd, N = fc.C, fc.Cdag, fc.N
    h = 0
    for i in range(3): h = h + Cd[i] * C[i + 1] + Cd[i + 1] * C[i]
    for i in range(3): h = h + 0.7 * N[i] * N[i + 1]
    h = h + 0.3 * (Cd[0] * Cd[2] + C[2] * C[0])
    assert h.is_hermitian()
    # A/Adag/F only, the form the MPS backends are handed
    assert multioperator.jordan_wigner(h).is_hermitian()
    bc = bosonchain.Bosonic_Chain(3, maxnb=[3, 3, 3], itensor_version="python")
    assert (bc.Adag[0] * bc.A[1] + bc.Adag[1] * bc.A[0]).is_hermitian()
    assert (bc.Adag[2] * bc.A[0] - bc.A[0] * bc.Adag[2]).is_zero()


# ----------------------------------------------- #2 the phase of ISy

def test_isy_dagger_is_minus_isy():
    """ISy is i*Sy, so its adjoint is -ISy, and daggering twice has to
    give it back."""
    sc = spinchain.Spin_Chain([2] * 3, itensor_version="python")
    ISy = sc.get_operator("ISy", 0)
    assert terms(ISy.get_dagger()) == terms(-ISy)
    assert terms(ISy.get_dagger().get_dagger()) == terms(ISy)
    assert (ISy + ISy.get_dagger()).is_zero()
    # ISy is still off the parity table, so no symbolic proof either way
    assert not ISy.is_hermitian()


def test_chain_hermiticity_probe_sees_the_phase_of_isy():
    """The chain-level check falls back to a numerical probe of
    op-op.get_dagger() for ISy, so it used to call ISy Hermitian and
    1j*ISy = -Sy not Hermitian, exactly the wrong way round."""
    sc = spin_chain_in_a_field(n=3)
    ISy = sc.get_operator("ISy", 0)
    assert not sc.is_hermitian(ISy)
    assert sc.is_hermitian(1j * ISy)


@pytest.mark.parametrize("submode,kwargs",
        [("KPM", dict(delta=0.2)), ("TD", dict(delta=0.4, dt=0.1))])
def test_isy_correlator_equals_its_1j_sy_spelling(submode, kwargs):
    """get_operator("ISy",0) and 1j*Sy[0] are the same operator, so they
    have to give the same correlator; with ISy in the first slot KPM
    (the default) and TD returned exactly minus the 1j*Sy one. Both
    spellings run on one chain, so they share the ground state."""
    sc = spin_chain_in_a_field(n=4)
    es = np.linspace(-1.0, 6.0, 40)
    ys = []
    for A in (sc.get_operator("ISy", 0), 1j * sc.Sy[0]):
        _x, y = sc.get_dynamical_correlator(mode="DMRG", submode=submode,
                name=[A, sc.Sz[1]], es=es, **kwargs)
        ys.append(np.asarray(y))
    assert np.max(np.abs(ys[1])) > 1e-2 # a sign flip would be visible
    assert np.max(np.abs(ys[0] - ys[1])) < 1e-6
