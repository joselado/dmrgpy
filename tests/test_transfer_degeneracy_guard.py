"""The dominant-fixed-point guard must reject ties, not periodicity.

`idmrg._check_dominant_eigenvalue_nondegenerate` exists to stop a single
"dominant eigenvector" being picked arbitrarily out of a genuinely tied
eigenspace -- the cat-state case `idmrg.imps_sum` can produce. It used to
test for that by MAGNITUDE, which also rejects a perfectly well-posed
period-2 state: a transfer matrix's peripheral spectrum is
`rho * exp(2 pi i k / p)`, so a period-2 state carries an eigenvalue at
`-rho`, distinct from `+rho` and with its own distinct eigenvector.

That is exactly what a half-filled free-fermion chain has (it is critical
with `2k_F = pi`, so the phase-pi eigenvalue's magnitude approaches 1 as
the correlation length diverges). Measured on that model at D=16, every
firing of the magnitude test had the signature |lambda| = (1, 0.999999999)
with arg = (0, +-pi), and it took out ~42% of individual VUMPS attempts and
~16% of whole ground-state solves.
"""
import numpy as np
import pytest

from dmrgpy.pyitensor import idmrg


def _check(w, perron):
    return idmrg._check_dominant_eigenvalue_nondegenerate(
        np.array(w, dtype=complex), "test", perron=perron)


def test_period_two_spectrum_is_accepted():
    """+1 and -1 are distinct eigenvalues: well-posed, must not raise."""
    order = _check([1.0, -0.999999999, 0.5], perron=True)
    assert np.isclose(np.array([1.0, -0.999999999, 0.5])[order[0]], 1.0)


def test_period_two_picks_the_perron_root_even_when_it_is_not_first():
    """The magnitude tie must be broken toward the positive eigenvalue --
    the only member of the peripheral spectrum whose eigenvector is the
    positive(-semidefinite) density matrix every caller divides by the
    trace of. Here -1 has the LARGER magnitude, so a plain argsort by
    magnitude would return it."""
    w = [-1.0, 0.999999999, 0.5]
    order = _check(w, perron=True)
    assert np.real(np.array(w)[order[0]]) > 0


def test_a_genuine_tie_still_raises():
    """The cat-state case the guard exists for: the SAME eigenvalue twice."""
    with pytest.raises(RuntimeError):
        _check([1.0, 1.0, 0.5], perron=True)


def test_a_genuine_near_tie_still_raises():
    with pytest.raises(RuntimeError):
        _check([1.0, 1.0 - 1e-13, 0.5], perron=True)


def test_complex_conjugate_peripheral_pair_is_accepted():
    """A period-3-like state: exp(+-2 pi i/3) at equal magnitude, with a
    real positive Perron root that must be selected."""
    w = [np.exp(2j * np.pi / 3), np.exp(-2j * np.pi / 3), 1.0, 0.2]
    order = _check(w, perron=True)
    assert np.isclose(np.array(w)[order[0]], 1.0)


def test_the_mixed_overlap_caller_keeps_the_magnitude_test():
    """_dominant_eigenvalue_mixed extracts a per-unit-cell factor |eta|^N,
    which oscillates rather than converging when two eigenvalues share a
    magnitude -- so there a magnitude tie IS ambiguous and must still
    raise, period-2 or not."""
    with pytest.raises(RuntimeError):
        _check([1.0, -0.9999999999, 0.5], perron=False)


def test_the_fixed_point_of_a_period_two_state_is_positive():
    """End to end through _dominant_fixed_point rather than the guard
    alone: a transfer matrix built to have +1 and -1 peripheral eigenvalues
    must yield the POSITIVE fixed point, since callers treat it as a
    density matrix and normalize by its trace."""
    chi = 2
    # E4 has indices (l,L,r,R) and acts rho[r,R] -> out[l,L], so an
    # elementwise map out[a,b] = M[a,b]*rho[a,b] is E4[a,b,c,d] =
    # M[a,b]*delta_ac*delta_bd. Its eigenvalues are M's entries and its
    # eigenvectors the elementary matrices, so M[0,0]=1 gives a Perron
    # eigenvector [[1,0],[0,0]] -- positive semidefinite, trace 1 -- and
    # M[1,1] puts a distinct peripheral eigenvalue at -1.
    M = np.array([[1.0, 0.4], [0.4, -0.999999999]])
    E = np.zeros((chi, chi, chi, chi), dtype=complex)
    for a in range(chi):
        for b in range(chi):
            E[a, b, a, b] = M[a, b]
    rho, eta = idmrg._dominant_fixed_point(
        [E], "right", "test", "unused message {}")
    assert np.real(eta) > 0
    ev = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    assert np.min(ev) > -1e-10, "fixed point is not positive semidefinite"
