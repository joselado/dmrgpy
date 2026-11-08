"""cvm_solver="variational" (Jeckelmann's dynamical DMRG) on
itensor_version="julia_live".

The sweep itself is ddmrg.jl, a port of mpscpp3's
Chain::ddmrg_correction_vector. This file builds its inputs in the units
of the solver scale rule and turns the functional into the spectral
weight, exactly as v3 does.

Scale rule: s is the power of two that brings the largest coefficient of
H below 1 into [1,2), and 1 otherwise (nhdmrg._unit_scale_up). The solver
sees K_s = s*(H-E0-omega), eta_s = s*eta and b_s = -eta_s*B|GS>, so
M_s = K_s^2 + eta_s^2 = s^2 M, b_s = s b, the minimizer is x_s = x/s and
W_s = W, which gives C = -W_s/(pi*eta) at any s. Because s is a power of
two, every product with it is exact and at s = 1 nothing changes.
"""
import numpy as np

from .juliasession import Main as Mainjl


def solver_scale(self):
    """The unit-scale factor of the chain's Hamiltonian."""
    from .. import nhdmrg
    terms = self.hamiltonian.to_terms()
    return nhdmrg._unit_scale_up(nhdmrg._max_abs_coef(terms))


def correction_vector(self, B, omega, eta, x0=None, s=None, wf0=None,
                      e0=None):
    """(value, x) at one frequency: the spectral weight -W_min/(pi*eta)
    of <GS|B^dagger (z-H)^-1 B|GS> in dynamics.py's convention, and the
    correction vector (a Julia MPS, in solver units) to warm start the
    next frequency. The caller has checked A = B^dagger. wf0 and e0 are
    the ground state and its energy; without them the chain's own are
    taken (get_gs() also sets self.e0)."""
    if s is None: s = solver_scale(self)
    if wf0 is None: wf0 = self.get_gs()
    if e0 is None: e0 = self.e0
    K = self.toMPO(s*(self.hamiltonian - (e0 + omega)))
    S2 = self.toMPO(B)
    eta_s = s*eta
    W, x = Mainjl.ddmrg_correction_vector(K.jlmpo, S2.jlmpo, wf0.jlmps,
            float(eta_s), float(-eta_s), int(self.cvm_maxm),
            int(self.cvm_nsweeps), float(self.cutoff), int(self.mpomaxm), x0)
    return -float(W)/(np.pi*eta), x


def dynamical_correlator_ddmrg(self, A, B, es, delta, wf0=None, e0=None):
    """Frequency sweep, each point warm started from the previous one's
    correction vector, as cvm.dynamical_correlator_ddmrg does on the
    session backends. The ground state is taken once, before the loop,
    and the same wf0 and e0 serve every frequency."""
    s = solver_scale(self)
    if wf0 is None: wf0 = self.get_gs()
    if e0 is None: e0 = self.e0
    out, x = [], None
    for e in es:
        o, x = correction_vector(self, B, e, delta, x0=x, s=s, wf0=wf0,
                                 e0=e0)
        print("DDMRG in E = ", e, " value = ", o)
        out.append(o)
    return (es, np.array(out))
