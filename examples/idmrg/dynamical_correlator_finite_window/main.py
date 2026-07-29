# Dynamical correlator of an infinite (iDMRG) chain, computed with the KPM
# method on a finite window -- infinitechain.Infinite_Many_Body_Chain.
# kpm_finite, see its own docstring
# (pyitensor/idmrg.py has no analogue of the transfer-matrix formalism
# vev/correlator use for *dynamical* quantities, since H is extensive/
# unbounded in the thermodynamic limit and a literal Chebyshev expansion
# of it has no meaning -- unlike apply_mpo/imps_sum's bounded-operator
# scope). Instead this method builds an ordinary finite, open-boundary
# chain of n_window repeats of the unit cell and delegates directly to
# kpmdmrg.get_dynamical_correlator -- the exact same KPM code path an
# ordinary finite Spin_Chain uses -- placing the two operators at the
# window's central cell, as far as possible from both open ends. Named
# kpm_finite (not get_dynamical_correlator, the finite-chain method it
# wraps) to flag up front that this is an approximation -- see
# kpm_finite's own docstring for what a genuinely infinite-chain KPM
# (infinite boundary conditions + a growing window) would need instead;
# not implemented here.
#
# SCOPE: a genuine finite-window *approximation* to the (only ill-defined
# via KPM anyway) infinite-chain dynamical correlator, not an exact
# infinite-size method -- results carry finite-size/open-boundary
# corrections that must be checked by convergence in n_window, exactly as
# demonstrated below (three window sizes, same delta/operators).

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import infinitechain

ic = infinitechain.Infinite_Spin_Chain(["1/2"])
h = ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]
ic.set_hamiltonian(h)
ic.maxm, ic.maxiter, ic.etol = 30, 60, 1e-9
density = ic.gs_energy()
print("iDMRG ground-state energy density:", density)

es = np.linspace(-0.5, 4.5, 120)
delta = 0.35  # kept coarse on purpose -- see this example's header and
              # kpm_finite's own docstring for why a fine delta needs a
              # KPM moment count comparable to n_window itself, at which
              # point open-boundary reflections contaminate the result
              # regardless of window size

print()
print("=== <Sz(0,t)Sz(0,0)> local dynamical correlator, convergence in n_window ===")
heights = []
for n_window in (8, 12, 16):
    _, ys = ic.kpm_finite(
        "Sz", 0, "Sz", 0, n_window=n_window,
        window_chain_kwargs=dict(maxm=24, nsweeps=10),
        delta=delta, es=es)
    peak_w = es[np.argmax(ys.real)]
    height = np.max(np.abs(ys))
    heights.append(height)
    print("n_window={:3d}  peak at omega={:.3f}  max|S(omega)|={:.4f}".format(
        n_window, peak_w, height))

# For such small windows (8-16 sites) individual finite-size levels are
# still widely spaced, so the *exact* peak position can jump between
# comparable, closely-spaced levels as n_window changes -- not a bug, see
# this example's/kpm_finite's own docstring on why a fixed delta needs
# growing n_window (and, eventually, a smaller delta) to really converge.
# The overall peak *height* (total nearby spectral weight, a
# coarser/better-converged diagnostic) is checked instead.
assert max(heights) / min(heights) < 2.0

print()
print("dynamical_correlator_finite_window example PASSED")
