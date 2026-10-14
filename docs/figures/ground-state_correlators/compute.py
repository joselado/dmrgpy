"""Figure 2, ground-state_correlators: decay of <Sz_i Sz_{i+r}> in the
spin-1/2 and in the spin-1 Heisenberg chain, power law against
exponential.

Run with  bash docs/figures/run.sh ground-state_correlators

Both chains are open, n sites, solved in the sector Sz_tot=0, and the pairs
go from the site i0=n/2-16 outwards, r=1 to 32, so both sites stay 16 sites
or more away from the ends.

The spin-1 chain has a spin 1/2 at each end. A plain open spin-1 chain has
two free spin-1/2 edge states, and its four lowest states are split only by
about exp(-n/xi); DMRG lands on any mix of them, and the edge-edge part of
the correlator is of the same size as the bulk one at r=32. The spin 1/2 at
each end binds the edge spin into a singlet, and the ground state is unique
with the Haldane gap above it. The bulk is unchanged.

References (fitted in compute.py, stored below):
  spin 1/2: a sqrt(ln r)/r, the Luttinger-liquid form with its log
            correction, a fitted at r=8 to 32.
  spin 1:   b exp(-r/xi)/sqrt(r), the Ornstein-Zernike form of a gapped
            chain in 1D, with xi=6.03 fixed and b fitted at r=8 to 32.
            xi is also fitted freely, with this form and with a pure
            exponential, and both values are printed and stored.

Arrays in data.npz:
    r          (m,)   distances 1..rmax
    czz_half   (m,)   <Sz_i0 Sz_i0+r>, spin 1/2
    czz_one    (m,)   <Sz_i0 Sz_i0+r>, spin 1
    ref_half   (m,)   a sqrt(ln r)/r (zero at r=1)
    ref_one    (m,)   b exp(-r/xi)/sqrt(r), xi=6.03
    xi_fit_oz  ()     xi fitted freely with the Ornstein-Zernike form
    xi_fit_exp ()     xi fitted freely with a pure exponential
    e_half, e_one ()  ground-state energies
    params            JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain

params = common.start(dict(
    model="open Heisenberg chains, J=1: spin 1/2, and spin 1 with a spin 1/2"
          " at each end",
    n=64,
    rmax=32,
    maxm_half=120,
    maxm_one=100,
    nsweeps=20,
    cutoff=1e-12,
    sector="Sz=0",
    fit_window=[8, 32],
    xi_haldane=6.03,
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
n, rmax = params["n"], params["rmax"]
i0 = n//2 - 16
rs = np.arange(1, rmax+1)


def correlator(sites, maxm):
    t0 = time.time()
    sc = spinchain.Spin_Chain(sites, itensor_version=params["itensor_version"])
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    sc.set_hamiltonian(h)
    sc.maxm = maxm
    sc.nsweeps = params["nsweeps"]
    sc.cutoff = params["cutoff"]
    sc.set_conserved_sector(Sz=0)
    e = sc.gs_energy()
    c = np.array([sc.vev(sc.Sz[i0]*sc.Sz[i0+r]).real for r in rs])
    print("%s: E0=%.10f  (%.1f s)" % (sites[1], e, time.time()-t0), flush=True)
    for r, v in zip(rs, c):
        print("   r=%2d  <SzSz>=% .6e" % (r, v))
    return e, c


e_half, c_half = correlator(["S=1/2"]*n, params["maxm_half"])
e_one, c_one = correlator(["S=1/2"] + ["S=1"]*(n-2) + ["S=1/2"],
                          params["maxm_one"])

lo, hi = params["fit_window"]
win = (rs >= lo) & (rs <= hi)
# spin 1/2: a sqrt(ln r)/r, a from a log-space mean over the window
f_half = np.sqrt(np.log(rs))/rs
a = np.exp(np.mean(np.log(np.abs(c_half[win])) - np.log(f_half[win])))
ref_half = a*f_half
print("spin 1/2: a=%.4f, |C|/(a sqrt(ln r)/r) over r=%d..%d: min %.3f max %.3f"
      % (a, lo, hi, np.min(np.abs(c_half[win])/ref_half[win]),
         np.max(np.abs(c_half[win])/ref_half[win])))
# spin 1: b exp(-r/xi)/sqrt(r) with xi fixed, and xi fitted freely
xi = params["xi_haldane"]
f_one = np.exp(-rs/xi)/np.sqrt(rs)
b = np.exp(np.mean(np.log(np.abs(c_one[win])) - np.log(f_one[win])))
ref_one = b*f_one
y = np.log(np.abs(c_one[win]))
xi_fit_oz = -1.0/np.polyfit(rs[win], y + 0.5*np.log(rs[win]), 1)[0]
xi_fit_exp = -1.0/np.polyfit(rs[win], y, 1)[0]
print("spin 1: b=%.4f, xi fitted with exp(-r/xi)/sqrt(r): %.3f, with "
      "exp(-r/xi): %.3f (reference 6.03)" % (b, xi_fit_oz, xi_fit_exp))
print("spin 1: |C|/ref over r=%d..%d: min %.3f max %.3f"
      % (lo, hi, np.min(np.abs(c_one[win])/ref_one[win]),
         np.max(np.abs(c_one[win])/ref_one[win])))

common.save(__file__, params, r=rs, czz_half=c_half, czz_one=c_one,
            ref_half=ref_half, ref_one=ref_one,
            xi_fit_oz=np.array(xi_fit_oz), xi_fit_exp=np.array(xi_fit_exp),
            e_half=np.array(e_half), e_one=np.array(e_one))
