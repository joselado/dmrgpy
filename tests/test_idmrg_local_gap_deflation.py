"""The deflated excited-state solve behind `local_excitation_gap`
(`pyitensor/idmrg.py`'s `_lowest_two_eigenvalues`, and its ITensor v3 C++
port inside `Chain::idmrg_local_excitation_gap`).

The bug these lock in: the excited state was found as the smallest
eigenvalue of `P H P` with the bare projector `P = I - |psi0><psi0|`. That
operator agrees with `H` on `psi0`'s orthogonal complement, but it also
carries `psi0` itself at eigenvalue *exactly zero* -- so as soon as `H`'s
own ground eigenvalue is POSITIVE, zero is the smallest eigenvalue of
`P H P` and is precisely what a smallest-eigenvalue solver returns.
Deflation excludes `psi0` from the Krylov space in exact arithmetic only;
rounding regrows the component and a restarted solver locks onto it. The
reported gap is then `0 - e0`, i.e. exactly minus the stored superblock
energy -- reported from an external Kondo-chain session at `maxm=16` as
-358.003 meV against a stored energy of +0.358002587681, digit for digit.

Whether that stored energy is positive is nobody's choice: the per-site
energy baseline (`_subtract_energy_baseline`) leaves it as a small residual
boundary term of either sign, which is why the same model was fine at
`maxm=8` (energy -1.05, negative) and wrong at `maxm=16`. So the sign of a
number no test looked at decided whether the answer was right -- hence a
test that sweeps it deliberately.

Testing `_lowest_two_eigenvalues` against a dense random Hermitian matrix
of known spectrum, rather than through a chain, is what makes that sweep
possible at all: the stored superblock's energy cannot be dialled from the
public API, while `H + shift*I` dials it exactly.
"""
import numpy as np
import pytest

from dmrgpy.pyitensor.idmrg import _lowest_two_eigenvalues


def _hermitian(n, seed):
    rng = np.random.default_rng(seed)
    m = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    return (m + m.conj().T) / 2


# The middle two shifts put the ground eigenvalue just above zero and well
# above it -- the regime the bare-projector deflation got wrong.
@pytest.mark.parametrize("ground", [None, 0.35, 5.0, 50.0])
def test_two_lowest_eigenvalues_are_exact_at_any_ground_energy_sign(ground):
    n = 200
    h0 = _hermitian(n, seed=0)
    w, v = np.linalg.eigh(h0)
    shift = 0.0 if ground is None else ground - w[0]
    h = h0 + shift * np.eye(n)
    e0, e1 = _lowest_two_eigenvalues(lambda x: h @ x, v[:, 0].copy(), n, 200)
    assert e0 == pytest.approx(w[0] + shift, abs=1e-8)
    assert e1 == pytest.approx(w[1] + shift, abs=1e-8)
    assert e1 - e0 >= 0.0


def test_gap_is_invariant_under_a_constant_shift_of_the_operator():
    """The whole point: a constant added to a Hermitian operator moves both
    eigenvalues together, so the gap must not move. The old code's answer
    moved by exactly the shift once the ground energy crossed zero."""
    n = 150
    h0 = _hermitian(n, seed=3)
    w, v = np.linalg.eigh(h0)
    gaps = []
    for ground in (-2.0, 0.0, 2.0, 20.0):
        h = h0 + (ground - w[0]) * np.eye(n)
        e0, e1 = _lowest_two_eigenvalues(lambda x: h @ x, v[:, 0].copy(), n, 200)
        gaps.append(e1 - e0)
    assert gaps == pytest.approx([w[1] - w[0]] * len(gaps), abs=1e-8)


def test_recovers_the_true_ground_state_when_handed_an_excited_one():
    """`psi0_stored` is only a candidate: if the growing algorithm's own
    local solve stopped on an excited eigenvalue (its convergence test
    bounds the distance to *some* eigenvalue, not to the lowest), the
    random-started solve must still find the real ground state, and the
    returned pair must still be the two lowest."""
    n = 120
    h = _hermitian(n, seed=7)
    w, v = np.linalg.eigh(h)
    e0, e1 = _lowest_two_eigenvalues(lambda x: h @ x, v[:, 5].copy(), n, 200)
    assert e0 == pytest.approx(w[0], abs=1e-8)
    assert e1 == pytest.approx(w[1], abs=1e-8)
