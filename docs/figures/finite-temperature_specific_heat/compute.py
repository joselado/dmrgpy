"""Figure 19, finite-temperature_specific_heat: thermal energy and specific
heat of the spin-1/2 Heisenberg chain from purification and from METTS,
against the exact result from the full spectrum.

Run with  bash docs/figures/run.sh finite-temperature_specific_heat

Arrays in data.npz (energies total, in units of J; divide by n per site):
    T          (k,)  temperatures of the purification points
    E_pur      (k,)  purification <H>_T
    C_pur      (k,)  purification C = (<H^2>_T-<H>_T^2)/T^2, with the
                     variance from gs_energy_fluctuation() (||(H-<H>)psi||^2)
    C_pur_npow (k,)  the same from vev(h,npow=2)-vev(h)^2, as a cross-check
    T_metts    (m,)  temperatures of the METTS points
    E_metts    (m,)  METTS <H>_T
    err_metts  (m,)  METTS statistical error of <H>_T
    T_exact    (K,)  fine temperature grid of the exact curves
    E_exact    (K,)  exact <H>_T from eigvalsh of the ED Hamiltonian
    C_exact    (K,)  exact C(T)
    params           JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain, thermal

params = common.start(dict(
    model="spin-1/2 Heisenberg chain, open, J=1",
    n=10,
    T_purification="np.geomspace(0.1, 5, 20)",
    purification_maxm=100,
    T_metts=[0.25, 0.5, 1.0, 2.0, 4.0],
    metts_nsamples=300,
    metts_nwarmup=30,
    metts_seed=1,
    T_exact="np.geomspace(0.05, 10, 300)",
    itensor_version=3,
))
n = params["n"]
spins = ["S=1/2"] * n


def heisenberg(c):
    return sum(c.Sx[i]*c.Sx[i+1] + c.Sy[i]*c.Sy[i+1] + c.Sz[i]*c.Sz[i+1]
               for i in range(n-1))


# exact reference from the full spectrum (1024 states)
sc = spinchain.Spin_Chain(spins, itensor_version=params["itensor_version"])
h = heisenberg(sc)
sc.set_hamiltonian(h)
hm = sc.get_ED_obj().get_hamiltonian()
hm = hm.toarray() if hasattr(hm, "toarray") else np.asarray(hm)
ev = np.linalg.eigvalsh(hm)
print("ED spectrum: %d states, E0 = %.10f" % (len(ev), ev[0]), flush=True)


def exact(T):
    w = np.exp(-(ev - ev.min())/T)
    z = w.sum()
    e = (ev*w).sum()/z
    e2 = (ev**2*w).sum()/z
    return e, (e2 - e**2)/T**2


T_exact = np.geomspace(0.05, 10, 300)
E_exact, C_exact = np.array([exact(T) for T in T_exact]).T

# purification: a fresh Thermal_Spin_Chain per temperature
T_pur = np.geomspace(0.1, 5, 20)
E_pur, C_pur, C_npow = [], [], []
for T in T_pur:
    t0 = time.time()
    tc = thermal.Thermal_Spin_Chain(spins, T=T,
                                    itensor_version=params["itensor_version"])
    tc.MBChain.maxm = params["purification_maxm"]
    ht = heisenberg(tc)
    tc.set_hamiltonian(ht)
    tc.get_gs()
    e = tc.MBChain.vev(ht).real
    de = tc.MBChain.gs_energy_fluctuation()
    e2 = tc.MBChain.vev(ht, npow=2).real
    E_pur.append(e)
    C_pur.append(de**2/T**2)
    C_npow.append((e2 - e**2)/T**2)
    ee, ce = exact(T)
    print("purification T=%.4f  E=%.8f (exact %.8f)  C=%.6f npow %.6f "
          "(exact %.6f)  (%.1f s)" % (T, e, ee, C_pur[-1], C_npow[-1], ce,
                                       time.time() - t0), flush=True)

# METTS on the physical chain
E_m, err_m = [], []
for T in params["T_metts"]:
    t0 = time.time()
    m, err = sc.metts_vev(h, T, nsamples=params["metts_nsamples"],
                          nwarmup=params["metts_nwarmup"],
                          seed=params["metts_seed"])
    E_m.append(m.real)
    err_m.append(err)
    ee, _ = exact(T)
    print("METTS T=%.2f  E=%.6f +- %.6f (exact %.6f, %.1f error bars)  "
          "(%.1f s)" % (T, m.real, err, ee, abs(m.real - ee)/err,
                        time.time() - t0), flush=True)

common.save(__file__, params, T=T_pur, E_pur=np.array(E_pur),
            C_pur=np.array(C_pur), C_pur_npow=np.array(C_npow),
            T_metts=np.array(params["T_metts"]), E_metts=np.array(E_m),
            err_metts=np.array(err_m), T_exact=T_exact, E_exact=E_exact,
            C_exact=C_exact)
