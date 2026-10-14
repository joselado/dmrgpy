"""Figure 13, dynamical-correlators_sqw: dynamical spin structure factor
S^zz(q,omega) of the open spin-1/2 Heisenberg chain, from one KPM run per
momentum with the standing-wave operator of the open chain,

    S^z_q = sqrt(2/(n+1)) sum_i sin(q (i+1)) S^z_i,   q_k = pi k/(n+1).

Run with  bash docs/figures/run.sh dynamical-correlators_sqw

Arrays in data.npz:
    qs        (nq,)       momenta q_k, k=1..n
    es        (ne,)       frequency grid of the map, in J
    sqw       (nq, ne)    S^zz(q,omega) on es, Jackson kernel, real part
    sumrule   (nq,)       int S(q,omega) domega over the full spectrum
    sq        (nq,)       <GS|S^z_q S^z_q|GS> from vev, the sum-rule reference
    npol      (nq,)       number of Chebyshev moments used per q
    e0        ()          DMRG ground-state energy, total, in J
    calib_kpmmaxm (c,)    kpmmaxm values of the calibration at q_n
    calib_sumrule (c,)    sum rule at q_n for each of them
    calib_diff    (c,)    max |S - S_previous| on es, against the previous kpmmaxm
    params                JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain, kpmdmrg

params = common.start(dict(
    model="spin-1/2 Heisenberg chain, open, J=1",
    n=int(os.environ.get("FIG_N", 32)),
    maxm=80,
    nsweeps=14,
    cutoff=1e-12,
    delta=0.15,
    kernel="jackson",
    calib_kpmmaxm=[40, 80, 120],
    calib_tol=0.01,
    es=[0.0, 4.0, 400],
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
n = params["n"]
delta = params["delta"]
es = np.linspace(*params["es"][:2], params["es"][2])

sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=params["itensor_version"])
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
sc.set_hamiltonian(h)
sc.maxm = params["maxm"]
sc.nsweeps = params["nsweeps"]
sc.cutoff = params["cutoff"]
t0 = time.time()
e0 = sc.gs_energy()
print("E0 = %.10f  E0/n = %.8f  (%.1f s)" % (e0, e0/n, time.time()-t0), flush=True)

qs = np.pi*np.arange(1, n+1)/(n+1)


def sq_operator(q):
    amp = np.sqrt(2.0/(n+1))*np.sin(q*np.arange(1, n+1))
    return sum(amp[i]*sc.Sz[i] for i in range(n))


def spectrum(Sq):
    """S(q,omega) on es, and its integral over the whole spectrum."""
    mus, emin, emax, scale, npol, d = sc.get_dynamical_correlator_moments(
        name=(Sq, Sq), delta=delta)
    _, y = kpmdmrg.dynamical_correlator_from_moments(
        mus, emin, emax, scale, npol, es, kernel=params["kernel"],
        delta=d, origin=sc.e0)
    ew = np.linspace(-1.0, emax - sc.e0 + 1.0, 20000)
    _, yw = kpmdmrg.dynamical_correlator_from_moments(
        mus, emin, emax, scale, npol, ew, kernel=params["kernel"],
        delta=d, origin=sc.e0)
    return y.real, np.trapezoid(yw.real, ew), npol


# calibration of kpmmaxm at q close to pi, where the weight piles up
Sq = sq_operator(qs[-1])
ref = sc.vev(Sq*Sq).real
cal_m, cal_s, cal_d = [], [], []
prev = None
kpmmaxm = params["calib_kpmmaxm"][-1]
for m in params["calib_kpmmaxm"]:
    t0 = time.time()
    sc.kpmmaxm = m
    y, s, npol = spectrum(Sq)
    diff = np.nan if prev is None else np.max(np.abs(y - prev))
    cal_m.append(m); cal_s.append(s); cal_d.append(diff)
    print("calibration q=%.4f pi kpmmaxm=%3d  sum=%.5f  <SqSq>=%.5f  "
          "rel=%.2e  max|dS|=%.2e  npol=%d  (%.1f s)"
          % (qs[-1]/np.pi, m, s, ref, abs(s-ref)/ref, diff, npol,
             time.time()-t0), flush=True)
    if prev is not None and abs(s-ref)/ref < params["calib_tol"] \
            and diff < params["calib_tol"]*np.max(np.abs(y)):
        kpmmaxm = cal_m[-2]   # the previous one already agreed with this one
        break
    prev = y
print("chosen kpmmaxm =", kpmmaxm, flush=True)
params["kpmmaxm"] = kpmmaxm
sc.kpmmaxm = kpmmaxm

sqw, sumrule, sq, npols = [], [], [], []
for k, q in enumerate(qs):
    t0 = time.time()
    Sq = sq_operator(q)
    y, s, npol = spectrum(Sq)
    r = sc.vev(Sq*Sq).real
    sqw.append(y); sumrule.append(s); sq.append(r); npols.append(npol)
    print("k=%2d q=%.4f pi  sum=%.5f  <SqSq>=%.5f  rel=%.2e  npol=%d  (%.1f s)"
          % (k+1, q/np.pi, s, r, abs(s-r)/r, npol, time.time()-t0), flush=True)

common.save(__file__, params, qs=qs, es=es, sqw=np.array(sqw),
            sumrule=np.array(sumrule), sq=np.array(sq), npol=np.array(npols),
            e0=np.array(e0), calib_kpmmaxm=np.array(cal_m),
            calib_sumrule=np.array(cal_s), calib_diff=np.array(cal_d))
