"""Standalone regression coverage for pyitensor's KPM energy truncation
(src/dmrgpy/pyitensor/kpm_energy_truncation.py, Holzner et al., PRB 83,
195115 (2011), Sec. III-B), independent of the Chebyshev recursion it
plugs into (see test_kpm_energy_truncation_accuracy.py for the
end-to-end KPM-vs-ED coverage of that wiring).

Two tests below build a random MPS via Many_Body_Chain/pyitensor's own
_default_mps() (randomMPS()), which draws from the *global* np.random
state (see mpsalgebra.py's randomMPS() docstring/body -- no seed
parameter of its own). np.random.seed(0) is called as the very first
statement of each such test, before even building the chain, since
sc.get_gs() itself already consumes global random draws for DMRG's own
random starting state (pyitensor has no product-state start, see
CLAUDE.md's mpscpp3-specific ConserveQNs=false note) -- seeding any later
than that would still be deterministic *given* an unchanged call
sequence, but reseeding right before the chain even exists removes that
fragility entirely.
"""
import numpy as np

from dmrgpy import spinchain
from dmrgpy.pyitensor.autompo import AutoMPO
from dmrgpy.pyitensor.dmrg import _all_right_environments, _extend_left, one_site_heff
from dmrgpy.pyitensor.kpm_energy_truncation import energy_truncate
from dmrgpy.pyitensor.mpobuilder import to_mpo

_BUILD_CUTOFF = 1e-12

# kpm_scale < 0.5 is what actually shrinks Ws below the full bandwidth W in
# this codebase's convention (chain.py's _scaled_hamiltonian: scale =
# 1/((emax-emin)*kpm_scale), so Ws = (emax-emin)*kpm_scale -- *smaller*
# kpm_scale means a *smaller* effective window, the opposite of what the
# name might suggest). At 0.3, the true band edges rescale to about
# +-1.67, comfortably outside [-1,1].
_NARROW_KPM_SCALE = 0.3


def _heisenberg_chain(n=6):
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    sc.setup_python()
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    sc.get_gs()
    return sc


def _excited_state(sc):
    """A genuinely broad-in-energy state: Sz on site 0 applied to the
    ground state, the same construction kpm_dynamical_correlator() uses
    to build its two Chebyshev seed vectors."""
    session = sc._session
    terms = sc.Sz[0].to_terms()
    m1 = to_mpo(AutoMPO.from_terms(session.sites, terms), cutoff=_BUILD_CUTOFF,
                maxdim=session.mpomaxm)
    return session._apply_mpo_with(m1, session.wf0, session.cutoff, session.maxm)


def test_energy_truncate_noop_when_window_is_wide():
    """With Ws~W (today's safe default kpm_scale), every site's local
    Krylov spectrum should already sit inside the [-1,1] threshold, so
    truncation must leave the state numerically unchanged."""
    np.random.seed(0)
    sc = _heisenberg_chain(6)
    session = sc._session
    scaled_H, emin, emax, scale = session._scaled_hamiltonian(kpm_scale=0.7)
    psi = _excited_state(sc)
    avg_w, dchange = energy_truncate(psi, scaled_H, dK=10, n_sweeps=2, threshold=1.0)
    assert avg_w < 1e-6
    assert dchange < 1e-8


def test_energy_truncate_removes_high_energy_weight():
    """Rescale with a deliberately narrow effective window so H' has
    excited-state eigenvalues far outside [-1,1], then confirm truncation
    actually changes the state substantially, and that re-truncating its
    own output only changes it a little further -- most of the
    high-energy weight is gone after the first call. (The residual isn't
    expected to hit exact zero: the paper's own Sec. V.C notes energy
    truncation is "a precautionary measure, not a variational procedure"
    -- with a finite Krylov dimension dK, the projected spectrum is only
    a Ritz approximation to the site's true local spectrum, so a small,
    stable residual above threshold is expected and does not vanish with
    further truncation calls; see test_energy_truncate_bounds_krylov_
    spectrum below for the actual size of that residual.)"""
    np.random.seed(0)
    sc = _heisenberg_chain(6)
    session = sc._session
    scaled_H, emin, emax, scale = session._scaled_hamiltonian(kpm_scale=_NARROW_KPM_SCALE)
    # A random MPS (not the ground state or a mild local excitation of it)
    # has support spread across the *entire* spectrum, so shrinking the
    # window is guaranteed to leave sizeable weight above threshold --
    # unlike a single Sz insertion on the ground state, which stays too
    # close to the low-energy end to actually exercise the truncation.
    psi = session._default_mps()

    _, dchange = energy_truncate(psi, scaled_H, dK=30, n_sweeps=10, threshold=1.0)
    assert dchange > 1e-3  # truncation actually did something

    _, dchange2 = energy_truncate(psi, scaled_H, dK=30, n_sweeps=10, threshold=1.0)
    assert dchange2 < 1e-2 * dchange  # bulk of the work was the first call


def test_energy_truncate_bounds_krylov_spectrum():
    """After truncation, every site's local Krylov energy spectrum
    (recomputed independently of energy_truncate's own bookkeeping, via a
    *full* dense diagonalization of that site's local effective
    Hamiltonian rather than energy_truncate's own dK-dimensional Krylov
    approximation of it) must have no more than a couple percent of its
    weight above threshold -- a direct check of the postcondition Eqs.
    (37)-(39) are meant to establish, not just of the diagnostic counters
    energy_truncate happens to report. Not an exact-zero check: see
    test_energy_truncate_removes_high_energy_weight's docstring for why a
    small residual is expected with a finite Krylov dimension dK."""
    np.random.seed(0)
    sc = _heisenberg_chain(6)
    session = sc._session
    scaled_H, emin, emax, scale = session._scaled_hamiltonian(kpm_scale=_NARROW_KPM_SCALE)
    psi = session._default_mps()
    energy_truncate(psi, scaled_H, dK=30, n_sweeps=10, threshold=1.0)

    psi.position(1)
    env = _all_right_environments(scaled_H, psi)
    side_env = {0: (None, None)}
    n = psi.length()
    for i in range(1, n + 1):
        L, Lbra = side_env[i - 1]
        R, Rbra = env[i + 1]
        matvec, order_in, shape, x0 = one_site_heff(L, Lbra, scaled_H, psi, i, R, Rbra)
        dim = x0.size
        basis = np.eye(dim, dtype=complex)
        Hmat = np.column_stack([matvec(basis[:, k]) for k in range(dim)])
        Hmat = (Hmat + Hmat.conj().T) / 2
        evals, evecs = np.linalg.eigh(Hmat)
        coeff = evecs.conj().T @ x0
        weight_above_threshold = np.sum(np.abs(coeff[np.abs(evals) >= 1.0]) ** 2)
        total = np.vdot(x0, x0).real
        assert weight_above_threshold < 0.02 * total
        if i < n:
            psi.position(i + 1)
            side_env[i] = _extend_left(L, Lbra, scaled_H, psi, i)
