"""Exact free-fermion reference for the connected Sz-Sz real-time
dynamical correlator of the XX spin chain
H = J*sum_i(Sx_i Sx_{i+1} + Sy_i Sy_{i+1}) [+ dimerized J2 on alternating
bonds], via Jordan-Wigner + brute-force diagonalization of the
single-particle hopping matrix.

Why this is useful as a test oracle: many-body ED is exponential in system
size, but a *non-interacting* system's exact solution is polynomial (one
N x N diagonalization), so this stays exact for chains far too large for
ED -- the only oracle that can validate a genuinely large-N real-time
calculation like iDMRG's IBC-window TDVP (idmrg_window.py) at all.

Derivation. Under Jordan-Wigner, nearest-neighbor terms have no string:
S+_i S-_{i+1} + S-_i S+_{i+1} = c_i^dag c_{i+1} + c_{i+1}^dag c_i, and
Sx_iSx_{i+1}+Sy_iSy_{i+1} = (1/2)(S+_iS-_{i+1}+S-_iS+_{i+1}), so H maps to
a plain tight-binding chain with hopping J/2 (J2/2 on alternating bonds
for the dimerized case) -- diagonalize once (`np.linalg.eigh`), fill the
lowest N/2 single-particle levels (half filling, Sz=0 sector -- matches
dmrgpy's own zero-field XX ground state), and Sz_i = n_i - 1/2 maps
*locally*, no string.

dmrgpy's own "TD" submode convention (see idmrg_window.py's
`window_tdvp_step`/`snapshot_correlator` docstrings, and
`mpscpp3::quench_tdvp`'s `Hshift = H - EGS*Id`) evolves under the
ground-state-energy-shifted Hamiltonian, `S(x,y,t) =
<psi|Sz_x e^{-i(H-EGS)t} Sz_y|psi>`. Writing c_i = sum_k phi_k(i) a_k
(phi_k real eigenvectors) and using Wick's theorem for the resulting
4-fermion-operator matrix element (no anomalous <cc>/<c^dag c^dag> terms,
since H conserves particle number), the *connected* correlator collapses
to a clean, EGS-independent closed form (the EGS-dependent global phase
cancels exactly against Hshift|psi>=0):

    S_conn(x,y,t) = Gocc(x,y,t) * Gunocc(x,y,t)
    Gocc(x,y,t)   = sum_{k occupied}   phi_k(x) phi_k(y) exp(+i*e_k*t)
    Gunocc(x,y,t) = sum_{k unoccupied} phi_k(x) phi_k(y) exp(-i*e_k*t)

independent of any boundary/filling asymmetry (Friedel oscillations near
an open chain's edges cancel out of the *connected* part exactly, even
though they affect the raw <Sz_x><Sz_y> background -- see `sz_sz_raw`).

Validated directly against exact many-body ED (`dmrgpy.timedependent.
evolution_DC(mode="ED")`) on an 8-site open chain: matches to ~1e-12
(machine precision) at every (x,t) checked, *after* accounting for two
independent, pre-existing conventions in the ED/C++ backends that this
class's own formula does not need:

  1. `edtk/timedependent.py::evolution_DC` (and identically
     `mpscpp3::quench_tdvp`/`quench`) evolve *before* measuring on every
     loop iteration (`wf = evolve(wf,ht,dt); ...; cs.append(...)`), so
     `cs[0]` is actually the value at `t=dt`, not `t=0`, despite `ts =
     dt*arange(nt)` labeling it `t=0` -- a genuine, separate, pre-existing
     off-by-one in that convention's own `ts`/`cs` correspondence
     (confirmed directly: this class's own `t=0` value exactly matches
     `cs[0]` only after shifting by `+dt`). `idmrg_window.py`'s own
     `dynamical_correlator_td` does not have this bug (it measures
     *before* evolving each step -- confirmed by its own
     `test_dynamical_correlator_td_matches_exact_static_correlator_at_t0`)
     -- this class's `t=0` matches *that* convention's `ts[0]` directly,
     with no shift needed.
  2. The ED backend's own `S(x,t)` comes out as the *complex conjugate*
     of this class's own convention (a sign convention internal to
     `edtk/tdtk.py`'s `evolve()`, unrelated to the off-by-one above).
     `idmrg_window.py`'s own convention (forward evolution via
     `_lanczos_expm_multiply(matvec,x0,-1j*tau,...)`, the standard
     Schrodinger-picture sign) matches this class directly, with no
     conjugation needed -- confirmed directly via a dense-matrix exact
     cross-check of `window_tdvp_step` itself (see
     `test_idmrg_window_free_fermion.py`).

So: compare directly (no shift, no conjugate) against
`idmrg_window.dynamical_correlator_td`'s own `S(x,t)`; shift `t -> t+dt`
and conjugate before comparing against `evolution_DC(mode="ED")` or any
`itensor_version in (2,3)` DMRG "TD" submode result.
"""
import numpy as np


def build_hopping(N, J=1.0, J2=None):
    """N x N open-boundary tight-binding hopping matrix: amplitude J/2 on
    odd bonds (0-1, 2-3, ...), J2/2 on even bonds (1-2, 3-4, ...).
    J2=None (or J2=J) gives the uniform (gapless) XX chain; J2!=J gives a
    dimerized (SSH-like) chain -- still pure hopping, exactly solvable
    the same way, but with a dimerization gap that makes iDMRG
    convergence far better-behaved than the gapless uniform case (see
    `dynamical_correlator_td_reproducible_and_matches_dimerized_free_fermion`'s
    own docstring)."""
    if J2 is None:
        J2 = J
    h = np.zeros((N, N))
    for i in range(N - 1):
        t = (J if i % 2 == 0 else J2) / 2.0
        h[i, i + 1] = h[i + 1, i] = t
    return h


class FreeFermionXX:
    """Exact ground-state/real-time machinery for the half-filled XX (or
    dimerized XX) chain, via a single N x N diagonalization (brute force,
    polynomial in N -- valid for N far beyond many-body ED's exponential
    reach)."""

    def __init__(self, N, J=1.0, J2=None):
        if N % 2 != 0:
            raise ValueError(
                "FreeFermionXX: N must be even for exact half filling, got {}".format(N))
        self.N = N
        self.J = J
        h = build_hopping(N, J, J2)
        evals, evecs = np.linalg.eigh(h)  # ascending order
        self.evals = evals
        self.evecs = evecs
        self.n_occ = N // 2
        self.occ = slice(0, self.n_occ)
        self.unocc = slice(self.n_occ, N)

    def density(self, x):
        """<n_x> = sum_{k occ} |phi_k(x)|^2."""
        return np.sum(self.evecs[x, self.occ] ** 2)

    def sz(self, x):
        return self.density(x) - 0.5

    def _g(self, x, y, t, sl, sign):
        phi = self.evecs[x, sl] * self.evecs[y, sl]
        return np.sum(phi * np.exp(sign * 1j * self.evals[sl] * t))

    def g_occ(self, x, y, t):
        return self._g(x, y, t, self.occ, +1.0)

    def g_unocc(self, x, y, t):
        return self._g(x, y, t, self.unocc, -1.0)

    def sz_sz_connected(self, x, y, t):
        """S_conn(x,y,t) = <psi|Sz_x e^{-i(H-EGS)t} Sz_y|psi>_connected,
        for scalar or array t -- see this module's own docstring for the
        derivation and the sign/timing convention this matches."""
        return self.g_occ(x, y, t) * self.g_unocc(x, y, t)

    def sz_sz_raw(self, x, y, t):
        return self.sz_sz_connected(x, y, t) + self.sz(x) * self.sz(y)
