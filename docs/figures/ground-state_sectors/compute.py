"""Figure 3, ground-state_sectors: the magnetization curve M(h) of the XX
and of the Heisenberg chain from the sector energies E0(Sz).

Run with  bash docs/figures/run.sh ground-state_sectors

In a field, H - h Sz_tot, the ground state is in the sector M that
minimizes E0(M) - h M, so the magnetization steps from M to M+1 at
h_M = E0(M+1) - E0(M). One sector solve per M gives the whole curve.

References (computed here, stored below):
  XX chain, J (SxSx + SySy): free fermions with hopping J/2 on an open
  chain, single-particle levels J cos(pi k/(n+1)), k=1..n. E0(M) is the
  sum of the lowest n/2+M levels. This is exact at every M, and it is the
  check. The curve of the infinite chain, m = arcsin(h/J)/pi, is drawn too.
  Heisenberg chain: the last step is one magnon below the saturated
  state, h_s = E0(n/2) - E0(n/2-1) = J (1 + cos(pi/n)) on an open chain,
  which tends to 2J.

Arrays in data.npz:
    M         (k,)   total Sz, 0..n/2
    e_xx      (k,)   DMRG E0(M), XX chain
    e_heis    (k,)   DMRG E0(M), Heisenberg chain
    e_xx_ff   (k,)   free-fermion E0(M), XX chain
    hs_exact  ()     J (1 + cos(pi/n)), Heisenberg saturation field, open chain
    params           JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain

params = common.start(dict(
    model="open spin-1/2 chains, J=1: XX (SxSx+SySy) and Heisenberg",
    n=32,
    maxm=60,
    nsweeps=12,
    cutoff=1e-12,
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
n = params["n"]
Ms = np.arange(0, n//2 + 1)


def sector_energy(model, M):
    sc = spinchain.Spin_Chain(["S=1/2"]*n,
                              itensor_version=params["itensor_version"])
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
        if model == "heis":
            h = h + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = params["maxm"]
    sc.nsweeps = params["nsweeps"]
    sc.cutoff = params["cutoff"]
    sc.set_conserved_sector(Sz=2*int(M))    # Sz in units of 1/2
    return sc.gs_energy()


eps = np.sort(np.cos(np.pi*np.arange(1, n+1)/(n+1)))
e_xx_ff = np.array([np.sum(eps[:n//2 + M]) for M in Ms])
hs_exact = 1.0 + np.cos(np.pi/n)

e_xx, e_heis = [], []
for M in Ms:
    t0 = time.time()
    ex = sector_energy("xx", M)
    eh = sector_energy("heis", M)
    e_xx.append(ex)
    e_heis.append(eh)
    print("M=%2d  XX: E0=%.10f  FF=%.10f  diff=%.1e   Heisenberg: E0=%.10f"
          "  (%.1f s)" % (M, ex, e_xx_ff[M], ex - e_xx_ff[M], eh,
                          time.time()-t0), flush=True)
e_xx, e_heis = np.array(e_xx), np.array(e_heis)
print("XX: max|E0_DMRG - E0_FF| = %.2e" % np.max(np.abs(e_xx - e_xx_ff)))
print("XX: max|h_M DMRG - h_M FF| = %.2e"
      % np.max(np.abs(np.diff(e_xx) - np.diff(e_xx_ff))))
print("Heisenberg: h_s = E0(n/2)-E0(n/2-1) = %.8f, exact J(1+cos(pi/n)) = %.8f"
      % (e_heis[-1] - e_heis[-2], hs_exact))
print("Heisenberg steps h_M:", np.round(np.diff(e_heis), 4))

common.save(__file__, params, M=Ms, e_xx=e_xx, e_heis=e_heis,
            e_xx_ff=e_xx_ff, hs_exact=np.array(hs_exact))
