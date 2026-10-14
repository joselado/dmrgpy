"""Figure 4, excited-states_haldane: edge states and bulk gap of the open
spin-1 Heisenberg chain, from the sector ground states.

Run with  bash docs/figures/run.sh excited-states_haldane

On an open spin-1 chain of even length the two spin-1/2 edge states form a
singlet (Sz_tot=0) and a triplet (Sz_tot=1), split by about exp(-n/xi).
The bulk gap is one magnon on top of the edge triplet:
    edge splitting  E0(Sz=1) - E0(Sz=0)
    gap             E0(Sz=2) - E0(Sz=1)
(set_conserved_sector takes Sz in units of 1/2, so Sz=0, 2 and 4.)
The magnon of an open chain has momentum pi - pi/n, so the gap converges
as 1/n^2 to the gap of the infinite chain, 0.41048 J (White and Huse).
The extrapolation is extrapolate.size_extrapolation(ns**2, gaps, n=1), a
fit of a + b/n^2 weighted by n^2, on the sizes n >= fit_min_n. The next
correction is of order 1/n^3 (the effective length of the box is n plus a
constant), and with all sizes from n=16 the fit lands about 0.002 above
the reference; from n=48 on it is below 0.001.

Arrays in data.npz:
    ns          (k,)   chain lengths
    e0, e1, e2  (k,)   E0 in the sectors Sz_tot = 0, 1, 2
    split       (k,)   edge splitting E0(1) - E0(0)
    gap         (k,)   bulk gap E0(2) - E0(1)
    gap_inf     ()     extrapolated gap, a of a + b/n^2, sizes n >= fit_min_n
    gap_all     ()     the same fit with every size
    gap_slope   ()     b of a + b/n^2
    gap_ref     ()     0.41048, White and Huse
    xi_ref      ()     6.03, correlation length of the spin-1 chain
    profile     (n_p,) <Sz_i> in the Sz_tot=1 sector, n_p sites
    params             JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain
from dmrgpy import extrapolate

params = common.start(dict(
    model="open spin-1 Heisenberg chain, J=1",
    ns=[16, 24, 32, 48, 64, 96, 128],
    fit_min_n=48,
    n_profile=48,
    maxm=150,
    nsweeps=30,
    cutoff=1e-12,
    itensor_version=3,
    gap_ref=0.41048,
    xi_ref=6.03,
    seeds="none (ITensor v3 random initial MPS)",
))


def sector(n, sz2):
    sc = spinchain.Spin_Chain(["S=1"]*n,
                              itensor_version=params["itensor_version"])
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    sc.set_hamiltonian(h)
    sc.maxm = params["maxm"]
    sc.nsweeps = params["nsweeps"]
    sc.cutoff = params["cutoff"]
    sc.set_conserved_sector(Sz=sz2)
    e = sc.gs_energy()
    return sc, e


ns = np.array(params["ns"])
e0, e1, e2 = [], [], []
profile = None
for n in ns:
    t0 = time.time()
    _, a = sector(n, 0)
    sc1, b = sector(n, 2)
    _, c = sector(n, 4)
    e0.append(a)
    e1.append(b)
    e2.append(c)
    if n == params["n_profile"]:
        profile = np.array([sc1.vev(sc1.Sz[i]).real for i in range(n)])
    print("n=%3d  E0(0)=%.10f  E0(1)=%.10f  E0(2)=%.10f  split=%.3e  "
          "gap=%.6f  (%.1f s)" % (n, a, b, c, b-a, c-b, time.time()-t0),
          flush=True)
e0, e1, e2 = np.array(e0), np.array(e1), np.array(e2)
split = e1 - e0
gap = e2 - e1
fit = ns >= params["fit_min_n"]
gap_inf = extrapolate.size_extrapolation(ns[fit]**2, gap[fit], n=1)
gap_all = extrapolate.size_extrapolation(ns**2, gap, n=1)
print("extrapolated gap with every size: %.5f" % gap_all)
# the slope b of the same weighted fit, for the line on the plot
from numpy.polynomial.polynomial import polyfit
coef = polyfit(1.0/ns[fit]**2, gap[fit], 1, w=ns[fit]**2)
print("extrapolated gap (size_extrapolation) = %.5f, a+b/n^2 with a=%.5f "
      "b=%.3f; reference %.5f" % (gap_inf, coef[0], coef[1], params["gap_ref"]))
for n, g in zip(ns, gap):
    print("   n=%3d  gap=%.6f  fit=%.6f" % (n, g, coef[0] + coef[1]/n**2))
ok = split > 1e-8
xi = -1.0/np.polyfit(ns[ok], np.log(np.abs(split[ok])), 1)[0]
print("edge splitting decay length from the sizes with split > 1e-8: %.3f "
      "(reference %.2f)" % (xi, params["xi_ref"]))
if profile is not None:
    nh = params["n_profile"]//2
    print("profile: sum over the left half %.6f, right half %.6f, "
          "max |<Sz_i>| in the middle 8 sites %.2e"
          % (profile[:nh].sum(), profile[nh:].sum(),
             np.max(np.abs(profile[nh-4:nh+4]))))

common.save(__file__, params, ns=ns, e0=e0, e1=e1, e2=e2, split=split,
            gap=gap, gap_inf=np.array(gap_inf), gap_all=np.array(gap_all), gap_slope=np.array(coef[1]),
            gap_ref=np.array(params["gap_ref"]),
            xi_ref=np.array(params["xi_ref"]), profile=profile)
