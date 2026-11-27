"""The norm-growth gate of CLIK-MPS (pyitensor/clik.py, cliksession.py).

On the real axis exp(L t) is a quantum channel, so it contracts the trace
norm of any X, and the log-norm of a normalized vector gains at most
(1/2) ln D, D = <<1|1>>, per step and over a whole branch. The old gate
rejected any growth above 1, so a fast-purifying Lindbladian ended in
"reduce dt" at every sub-step count. On the complex contour there is no
such bound: a growth above 1 is redone at twice the sub-steps and kept
when both agree, and a branch taken over by one growing mode warns.
"""
import numpy as np
import pytest

from dmrgpy import spinchain, bosonchain
from dmrgpy.lindblad import Lindblad_Chain
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.mpsalgebra import inner
from dmrgpy.pyitensor.vectorize import identity_mps


def test_growth_verdicts():
    bound = 0.5*np.log(16)
    assert clik._growth(0.5, 0.0, bound) is None
    assert clik._growth(1.2, 0.0, bound) == "redo"
    assert clik._growth(1.5, 0.0, bound) == "reject"       # above the bound
    assert clik._growth(0.5, 1.0, bound) == "reject"       # the sum is
    assert clik._growth(2.7, 0.0, None) == "redo"          # complex contour
    assert clik._growth(np.inf, 0.0, None) == "reject"


def test_identity_under_strong_decay_is_accepted():
    """4 spins decaying at gamma=4 from the identity gain 1.117 per step,
    below (1/2) ln 16 = 1.386; the old gate rejected it at every count"""
    n, gamma = 4, 4.0
    sc = spinchain.Spin_Chain([2]*n, itensor_version="python")
    sc.maxm = 64
    sc.set_hamiltonian(sum(0.1*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1])
                           for i in range(n-1)))
    lc = Lindblad_Chain(sc, [np.sqrt(gamma)*(sc.Sx[i] - 1j*sc.Sy[i])
                             for i in range(n)])
    S = lc.doubled._session
    G = S._mpo(lc.L.to_terms())
    one = identity_mps(S.sites)
    bound = 0.5*np.log(abs(complex(inner(one, one))))
    assert bound == pytest.approx(0.5*np.log(16))
    states, sub = clik.evolve_branch(G, one, 0.5, 4, 1e-10, 64, meter=one,
                                     bounded=True, growth_bound=bound)
    assert len(states) == 4 and sub == 1


def test_complex_branch_taken_over_by_one_mode_warns():
    """seeded 3-site Bose-Hubbard at the default dt=2: an eigenmode outside
    the sector with Re(lambda z) = 3.98 fills the later snapshots"""
    L, N = 3, 2
    bc = bosonchain.Bosonic_Chain(L, maxnb=[N+1]*L, itensor_version="python")
    h = sum(-(bc.Adag[j+1]*bc.A[j] + bc.Adag[j]*bc.A[j+1]) for j in range(L-1))
    bc.set_hamiltonian(h + 0.5*sum(bc.N[j]*(bc.N[j]-1) for j in range(L)))
    bc.maxm = 64
    lc = Lindblad_Chain(bc, [(bc.Adag[j+1]+bc.Adag[j])*(bc.A[j+1]-bc.A[j])
                             for j in range(L-1)])
    start = ([1]*N + [0]*(L-N), [0]*(L-N) + [1]*N)
    with pytest.warns(RuntimeWarning, match="takes over"):
        lam = lc.get_spectrum(start=start)[0]
    lamE = lc.get_spectrum(mode="ED")[0]
    # the steady state and the slowest levels are still right
    assert np.min(np.abs(lam - lamE[0])) < 1e-6
