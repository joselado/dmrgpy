# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Timing comparison for the third-order Kondo two-time-correlator
# construction (kondospectrumtk/dmrgtwotime.py) across three ways of
# computing the exact same quantity: the ED excited-state-sum reference
# (conductance.third_order_kondo_dIdV), DMRG via the compiled
# itensor_version=3 extension (real TDVP time evolution), and DMRG via
# the pure-Python pyitensor backend (itensor_version="python", no
# compiled extension needed at all).
#
# dmrgtwotime.py never hardcodes which DMRG backend it's driving -- every
# call goes through chain._session/chain.toMPO/MPS's cpp_handle-based
# dispatch, the same generic mechanism every other itensor_version in
# (2,3,"python") call site in this codebase uses (see CLAUDE.md's
# "In-process pybind11 extension" section) -- so it was tried directly
# against a pyitensor chain and found to work completely unmodified. At
# matching (t2,tau) grid resolution the two DMRG backends' own G(t2,tau)
# constructions agree with the ED-eigenbasis G(t2,tau) construction to
# ~1e-9-1e-14 pointwise (see tests/test_kondo_spectrum_dmrgtwotime.py);
# what this example instead measures is convergence of the *physical*
# answer against the exact excited-state-sum reference as a function of
# (t2,tau) grid resolution, which is a separate, much coarser-converging
# question -- see the Gamma0/grid discussion below and
# examples/kondo_third_order_dmrg_VS_ED's module docstring for a
# from-scratch account of a wrong turn taken here (an under-resolved grid
# that made both DMRG backends "agree with each other" while both were
# equally wrong relative to the exact answer).
#
# Gamma0 is set much larger here (2e-3) than kondospectrumtk's own
# default (5e-6): the third-order Kondo term's t2/tau grid needs spacing
# finer than 1/omega0 and a *range* wider than several/Gamma0 (see
# dmrgtwotime.py's/twotime.py's module docstrings). At the default
# Gamma0, resolving that range needs on the order of 1e5-1e6 t2
# checkpoints -- infeasible for a live demo, each one its own real time
# evolution. Raising Gamma0 shrinks the required range proportionally,
# letting a modest grid (n_t2_half=n_tau_half=20 below) converge to
# within single-digit percent of the exact answer at practical cost --
# this is a legitimate physical choice (Gamma0 is just a lifetime-
# broadening parameter), not a numerical shortcut that changes what's
# being computed. See examples/kondo_third_order_dmrg_VS_ED for a plot of
# the resulting spectra (v3 only); this example is about wall-clock cost,
# not the physics itself.
import time
import numpy as np
from dmrgpy import spinchain, cppext
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_kondo_dIdV
from dmrgpy.kondospectrumtk.dmrgtwotime import two_time_kondo_term_dmrg

G = 2.0
MUB = 5.7883818066e-5 # eV/T


def build_chain(itensor_version):
    """Same 3-site chain as examples/kondo_third_order_dmrg_VS_ED and
    tests/test_kondo_spectrum_dmrgtwotime.py (ITensor v3's two-site
    DMRG/TDVP needs at least 3 sites; pyitensor has no such limitation,
    but the same chain is used for all three backends for a fair
    comparison)."""
    sc = spinchain.Spin_Chain(["1/2", "1/2", "1/2"], itensor_version=itensor_version)
    h = G*MUB*10.0*sc.Sz[0]
    for i in range(2):
        h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    return sc


eVs = np.linspace(-2e-3, 2e-3, 21)
omega0, Gamma0 = 2e-3, 2e-3 # see the module docstring above
n_t2_half, n_tau_half = 20, 20
dt2, dtau = 150.0, 150.0
T0, Jrho_s = 1.0, 1.0 # third_order_kondo_dIdV's own prefactor convention
                       # (4*pi*T0**2*Jrho_s*total) -- applied by hand
                       # below to two_time_kondo_term_dmrg's raw output,
                       # matching Spin_Chain._get_kondo_spectrum_dmrg

timings = {}
results = {}

# --- ED: exact eigenstate-sum construction --------------------------------
print("ED (exact eigenstate sum) ...")
sc_ed = build_chain(2) # itensor_version is irrelevant here: KondoSpectrum
                        # always does its own full diagonalization via
                        # get_ED_obj(), independent of the chain's DMRG
                        # backend setting
t0 = time.time()
ks = KondoSpectrum(sc_ed, site=0, T=0.0)
results["ED"] = third_order_kondo_dIdV(ks, eVs, Jrho_s, T0=T0, omega0=omega0,
                                        Gamma0=Gamma0)
timings["ED"] = time.time()-t0

# --- DMRG via the compiled itensor_version=3 extension ---------------------
if cppext.available(3):
    print("DMRG (itensor_version=3, real TDVP) ...")
    sc3 = build_chain(3)
    t0 = time.time()
    sc3.get_gs()
    term3 = two_time_kondo_term_dmrg(sc3, 0, eVs, omega0=omega0, Gamma0=Gamma0,
                                      dt2=dt2, n_t2_half=n_t2_half, dtau=dtau,
                                      n_tau_half=n_tau_half)
    results["DMRG (v3)"] = 4*np.pi*T0**2*Jrho_s*term3
    timings["DMRG (v3)"] = time.time()-t0
else:
    print("itensor_version=3 not compiled -- skipping that backend")

# --- DMRG via the pure-Python pyitensor backend -----------------------------
print("DMRG (itensor_version=\"python\", pyitensor) ...")
scpy = build_chain("python")
t0 = time.time()
scpy.get_gs()
termpy = two_time_kondo_term_dmrg(scpy, 0, eVs, omega0=omega0, Gamma0=Gamma0,
                                   dt2=dt2, n_t2_half=n_t2_half, dtau=dtau,
                                   n_tau_half=n_tau_half)
results["DMRG (pyitensor)"] = 4*np.pi*T0**2*Jrho_s*termpy
timings["DMRG (pyitensor)"] = time.time()-t0

# --- report -----------------------------------------------------------------
print()
for name, dt in timings.items():
    print("%-18s %8.2f s"%(name, dt))

ref = results["ED"]
diffs = {}
for name, val in results.items():
    if name=="ED": continue
    diffs[name] = np.max(np.abs(val-ref))
    print("%-18s max|diff vs ED| = %.4f (max|ED| = %.4f)"%(name, diffs[name], np.max(np.abs(ref))))
for name, diff in diffs.items():
    # a couple tens of percent (t2/tau-grid discretization error at this
    # resolution), not machine precision -- see
    # test_kondo_spectrum_dmrgtwotime.py for the much tighter
    # DMRG-vs-ED-at-matching-grid check (that one isolates the DMRG
    # construction itself from grid-resolution error, independent of how
    # well that grid resolution itself approximates the true answer)
    assert diff < 0.25*np.max(np.abs(ref)), name

import matplotlib.pyplot as plt

plt.bar(list(timings.keys()), list(timings.values()), color=["C0", "C1", "C2"][:len(timings)])
plt.ylabel("wall-clock time (s)")
plt.title("Third-order Kondo term: ED vs DMRG (v3, pyitensor)")
plt.show()
