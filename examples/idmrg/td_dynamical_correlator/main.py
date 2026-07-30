# Real-time dynamical correlator of an infinite (iDMRG) chain, via
# infinite-boundary-condition (IBC) window TDVP evolution --
# infinitechain.Infinite_Many_Body_Chain.td_dynamical_correlator, see its
# own docstring and pyitensor/idmrg_window.py's module docstring.
#
# Unlike kpm_finite (an ordinary finite, *open-boundary* window + KPM/
# Chebyshev -- see dynamical_correlator_finite_window/main.py in this same
# directory), this method caps the window's two ends with the *converged*
# iDMRG growth environment (idmrg_ground_state's own HL/HR, exposed on
# IDMRGResult -- see pyitensor/idmrg.py's own env_window_boundary comment)
# instead of plain open boundaries, following Milsted/Vanderstraeten et
# al., "Infinite boundary conditions for response functions and limit
# cycles in iDMRG" (arXiv:1804.09163). This removes the open-boundary
# artifacts (e.g. Friedel-oscillation-like contamination of the window's
# own central region) that no amount of n_window alone can fix for
# kpm_finite -- so a much smaller n_window margin should suffice here for
# comparable accuracy.
#
# The dynamics itself is real-time two-site TDVP (pyitensor/tdvp.py's own
# Krylov propagator, reused via this module's own window-aware sweep --
# see idmrg_window.py's module docstring for why tdvp.py's own sweep
# functions cannot be reused directly on a capped window), evolving a
# single perturbed window and reading off every distance x (and, via a
# spatial Fourier transform, every momentum k) from that one run -- the
# paper's own headline efficiency result, simplified here to t1=0 (see
# dynamical_correlator_td's own docstring for the full two-branch Eq. 7
# extension this does not implement).
#
# SCOPE: still an approximation (a finite window, however well-capped, is
# not literally the infinite chain) -- convergence in n_window (and in the
# TDVP truncation/Krylov parameters maxdim/cutoff/niter) should always be
# checked for quantitative work, demonstrated below.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import infinitechain

ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
h = 1.4*ic.SzC[0] + ic.SxC[0]*ic.SxR[0]  # gapped, transverse-field-like model
ic.set_hamiltonian(h)
ic.maxm, ic.maxiter, ic.etol, ic.niter = 20, 300, 1e-12, 150
density = ic.gs_energy()
print("iDMRG ground-state energy density:", density)
print("state_overlap (self-consistency diagnostic):", ic._result.state_overlap)

print()
print("=== S(k,omega) convergence in n_window ===")
# x_values is kept FIXED across the sweep (not scaled with n_window): the
# point of this check is whether *more environment margin* around a
# fixed observation region changes the result, not whether summing over
# more x itself does (a slowly, but still exponentially, decaying
# connected correlator can keep adding non-negligible contributions as
# the x-range grows, which is a separate convergence question -- see
# dynamical_correlator_td's own docstring on choosing x_values/n_window
# together for a fixed total simulated time).
ks = np.linspace(-np.pi, np.pi, 9)
es = np.linspace(-1, 6, 200)
xs_fixed = range(-4, 5)
peaks = []
for n_window in (10, 14, 18):
    _, _, Skw = ic.td_dynamical_correlator(
        "Sz", 0, "Sz", n_window=n_window, dt=0.05, nt=40,
        maxdim=60, cutoff=1e-10, niter=50, x_values=xs_fixed,
        ks=ks, es=es, delta=0.15, window=[-1, 6])
    peak = np.max(np.abs(Skw))
    peaks.append(peak)
    print("n_window={:3d}  max|S(k,omega)|={:.4f}".format(n_window, peak))
assert max(peaks) / min(peaks) < 1.5  # should stabilize as n_window grows

print()
print("=== cross-check vs kpm_finite (independent approximation) at r=1 ===")
# Same physical correlator (<Sz(0,t)Sz(1,0)>), same delta/window, via the
# two independent methods -- an exact match isn't expected (different
# approximation schemes, different systematic errors -- see
# examples/dynamical_correlator/dynamical_correlator_time_evolution/
# main.py for the same "compare, don't expect exact agreement" spirit
# between KPM and TD submodes on an ordinary *finite* chain), but both
# should agree on *where* the dominant spectral weight sits.
es_kpm, y_kpm = ic.kpm_finite("Sz", 0, "Sz", 1, n_window=20,
                                window_chain_kwargs=dict(maxm=30, nsweeps=10),
                                delta=0.3, es=np.linspace(-1, 6, 100))

from dmrgpy.pyitensor import idmrg_window
from dmrgpy.timedependent import _fourier_transform_correlator
ts, xs, S = idmrg_window.dynamical_correlator_td(
    ic._result, n_window=16, opname_A="Sz", opname_B="Sz", dt=0.05, nt=60,
    cutoff=1e-10, maxdim=60, niter=50, x_values=[1])
es_td, g_td = _fourier_transform_correlator(ts, S[:, 0], 0.05,
                                              es=np.linspace(-1, 6, 100),
                                              delta=0.3, window=[-1, 6])

peak_kpm = es_kpm[np.argmax(np.abs(y_kpm))]
peak_td = es_td[np.argmax(np.abs(g_td))]
print("kpm_finite peak at omega={:.3f}   td_dynamical_correlator peak at omega={:.3f}".format(
    peak_kpm, peak_td))
assert abs(peak_kpm - peak_td) < 1.5

print()
print("td_dynamical_correlator example PASSED")
