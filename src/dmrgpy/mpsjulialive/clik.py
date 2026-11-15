"""The MPS and MPO operations of CLIK-MPS (get_krylov_spectrum) on
itensor_version="julia_live".

cliksession.py runs pyitensor/clik.py's algorithm through an operations
object; its Session is the v3 one and this is the julia_live one, with the
same methods and the same meaning. A state is a bare Julia MPS and an
operator a bare Julia MPO: the driver only hands them back to these
methods. Everything is an existing function of this folder except the two
of clik.jl, the seeded starts and the expansion.

Where the Julia operation is not the v3 one, and a number can move:

- the random starts are ITensorMPS's random circuit states from a seeded
  MersenneTwister, not Gaussian draws, so the unconverged levels differ by
  draw, as they already do between v3 and "python";
- the local exponentials are KrylovKit's (tdvp.jl's tdvp_updater_kwargs:
  tol=1e-12, krylovdim=30 with restarts), so `niter` is not used;
- the expansion keeps a new direction at a bond when its projected weight
  is above 1e-12 of the trace-normalized density matrix of the Krylov
  vectors, where v3 and "python" compare against its norm;
- the MPO of G^dag G for exact_residuals is cut at the chain's cutoff and
  at max(maxm, mpomaxm), as on v3, where "python" builds it exactly.
"""

import numpy as np

from .juliasession import Main
from ..pyitensor import clik as _clik


def _number(z):
    """A Python float, or a Python complex when z has an imaginary part,
    so that a real state stays real under a real factor"""
    z = complex(z)
    return z.real if z.imag == 0 else z


class Session:
    """The MPS/MPO operations the driver needs, on a julia_live chain."""

    def __init__(self, chain, hermitian=True):
        self.chain = chain
        self.hermitian = bool(hermitian)
        self.maxm = int(chain.maxm)
        self.cutoff = float(chain.cutoff)
        self.mpomaxm = int(max(chain.maxm, chain.mpomaxm))

    def mpo(self, op):
        return self.chain.toMPO(op).jlmpo

    def inner(self, a, b):
        return complex(Main.overlap(a, b))

    def sandwich(self, a, A, b):
        return complex(Main.aMb_mpo(a, A, b))   # exact, nothing truncated

    def norm(self, psi):
        return float(np.sqrt(abs(self.inner(psi, psi))))

    def scale(self, psi, z):
        return Main.mpstimesscalar(_number(z), psi)

    def normalized(self, psi):
        return self.scale(psi, 1.0 / self.norm(psi))

    def add(self, a, b):
        return Main.summps(a, b, self.maxm, self.cutoff)

    def step(self, G, psi, z, niter):
        # tdvp_step applies exp(-i dt G) without normalizing, so dt = i z
        # gives exp(z G). herm is always passed: Lanczos for a Hermitian
        # G, Arnoldi otherwise, as tdvp_step and tdvp_arnoldi_step on v3
        return Main.tdvp_step(G, psi, complex(1j * z), self.cutoff,
                              self.maxm, herm=self.hermitian)

    def expand(self, G, psi):
        return Main.clik_expand(G, psi, int(_clik._GSE_ORDER), self.cutoff,
                                self.maxm)

    def dagger(self, G):
        return Main.dagger_mpo(G)

    def product(self, first, then):
        """The MPO of `then` after `first` (first acts first): mult_mpo
        takes the operator product A*B with B acting first."""
        return Main.mult_mpo(then, first, self.mpomaxm, self.cutoff)

    def apply(self, A, psi):
        return Main.gs_extra_apply(A, psi, self.maxm, self.cutoff)

    def random_starts(self, m, seed):
        p, u = Main.clik_random_starts(self.chain.jlsites, int(m), int(seed))
        return p, u

    def _not_ported(self, what):
        raise NotImplementedError(
            "get_krylov_spectrum on julia_live has random starts only; "
            + what + " is on itensor_version=\"python\" and 3")

    def product_state(self, vecs):
        self._not_ported("a product start")

    def identity_mps(self):
        self._not_ported("the traceless branch of a Lindbladian")

    def sector_identity_mps(self, q):
        self._not_ported("the traceless branch of a Lindbladian")
