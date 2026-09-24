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
import matplotlib.pyplot as plt
from dmrgpy import infinitechain

ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
ic.gs_method = "idmrg"
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
print("=== plotting S(k,omega) heatmap (most-converged n_window from above) ===")
# Reuse ks/es/Skw as left over from the loop's last (most-converged, n_window=18)
# iteration above -- no extra DMRG work needed. Skw is complex and shaped
# (len(ks), len(es)); plot np.abs(Skw), matching this same script's own
# convergence check above (np.max(np.abs(Skw))) and the established
# convention for the TD/TDVP-family dynamical correlator elsewhere in this
# codebase (examples/dynamical_correlator/dynamical_correlator_time_evolution/
# main.py plots its own "TD" submode with np.abs(...), while the separate KPM
# submode there uses .real instead -- the two submodes are not interchangeable
# conventions, they reflect what's actually real/physical for each method).
plt.figure(figsize=(6, 5))
plt.pcolormesh(ks, es, np.abs(Skw).T, shading="auto", cmap="viridis")
plt.colorbar(label=r"|S(k,$\omega$)|")
plt.xlabel("k")
plt.ylabel(r"$\omega$")
plt.title("iDMRG IBC-window TDVP dynamical correlator\n"
          "S(k,$\\omega$) = <Sz(0,t) Sz(x,0)> ($n_\\mathrm{window}$=%d)" % n_window)
plt.tight_layout()
plt.savefig("td_dynamical_correlator_Skw.png", dpi=150)
print("Plot saved to td_dynamical_correlator_Skw.png")
plt.show()

print()
print("=== cross-check vs kpm_finite (independent approximation), local Sx,Sx ===")
# The two independent methods on a correlator whose answer is known: Sx
# creates exactly one quasiparticle of this paramagnet, eps(k) =
# 2*sqrt(0.5525 + 0.35*cos k), so the LOCAL Sx,Sx spectrum (TD at x=0,
# kpm_finite at r=0) lives in the magnon band [0.9, 1.9], and at positive
# frequency only, D_n = E_n - E_0 > 0. An exact match isn't expected
# (different approximation schemes, different systematic errors -- see
# examples/dynamical_correlator/dynamical_correlator_time_evolution/
# main.py for the same "compare, don't expect exact agreement" spirit
# between KPM and TD submodes on an ordinary *finite* chain), but both
# must put the line inside the band. The TD side is the public route at
# x_values=[0], so it goes through sxt_to_skomega's own conjugation, and
# the grid is two-signed so that a mirrored spectrum would show.
#
# This block used to compare Sz,Sz at r=1, where kpm_finite (which
# subtracts no <A><B>) peaks on the elastic line at omega=0 and the TD
# route on a connected continuum, two different quantities; it passed
# only because every infinite-chain S(k,omega) was mirrored in omega at
# the time, which put the TD line next to the elastic one (2026-09-24b
# audit, finding 8).
es_loc = np.linspace(-3, 3, 301)
es_kpm, y_kpm = ic.kpm_finite("Sx", 0, "Sx", 0, n_window=20,
                                window_chain_kwargs=dict(maxm=30, nsweeps=10),
                                delta=0.3, es=es_loc)
y_kpm = np.asarray(y_kpm).real
_k, es_td, Sloc = ic.td_dynamical_correlator(
    "Sx", 0, "Sx", n_window=12, dt=0.1, nt=120, x_values=[0],
    maxdim=60, cutoff=1e-10, niter=50, ks=[0.0], es=es_loc, delta=0.3)
y_td = Sloc[0].real

peak_kpm = es_kpm[np.argmax(y_kpm)]
peak_td = es_td[np.argmax(y_td)]
a = np.abs(y_td)
w_neg = np.trapezoid(a[es_td < 0], es_td[es_td < 0]) / np.trapezoid(a, es_td)
print("kpm_finite peak at omega={:.3f}   td_dynamical_correlator peak at omega={:.3f}".format(
    peak_kpm, peak_td))
print("TD weight on omega<0: {:.3f}   shape correlation TD vs KPM: {:+.3f}".format(
    w_neg, np.corrcoef(y_td, y_kpm)[0, 1]))
assert 0.9 <= peak_kpm <= 1.9 and 0.9 <= peak_td <= 1.9
assert w_neg < 0.15

plt.figure(figsize=(6, 4))
plt.plot(es_kpm, y_kpm, label="kpm_finite, r=0")
plt.plot(es_td, y_td, "--", label="td_dynamical_correlator, x=0 (Re)")
plt.axvspan(0.9, 1.9, color="grey", alpha=0.15, label="magnon band")
plt.xlabel(r"$\omega$")
plt.ylabel(r"local $S^{xx}(\omega)$")
plt.title("local Sx,Sx spectrum, two independent methods")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("td_dynamical_correlator_local_crosscheck.png", dpi=150)
print("Plot saved to td_dynamical_correlator_local_crosscheck.png")
plt.show()

print()
print("td_dynamical_correlator example PASSED")
