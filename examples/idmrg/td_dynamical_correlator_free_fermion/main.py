# Cross-check iDMRG's real-time IBC-window dynamical correlator
# (pyitensor/idmrg_window.py's dynamical_correlator_td, wired up publicly
# as Infinite_Many_Body_Chain.td_dynamical_correlator) against an EXACT
# non-interacting (free-fermion) reference, on a system far larger than
# many-body ED could ever reach.
#
# WHY THIS WORKS: H = J*sum(Sx_i Sx_{i+1} + Sy_i Sy_{i+1}) [+ dimerized J2
# on alternating bonds] maps *exactly* to a free-fermion tight-binding
# chain under Jordan-Wigner (nearest-neighbor terms have no string), so
# its connected Sz-Sz real-time correlator is computable via a single
# N x N diagonalization -- polynomial cost, unlike many-body ED's
# exponential one -- for N in the thousands. See tests/
# _free_fermion_reference.py's own module docstring for the full
# derivation and dmrgpy's own convention this matches (a ground-state-
# energy-shifted Schrodinger-picture correlator, e^{-i(H-EGS)t}).
#
# This exact comparison is what found a real, previously-undiscovered bug
# in idmrg_window.py: window.env_HL/env_HR are not energy-baseline-
# subtracted, so window_tdvp_step's TDVP evolution used to carry an
# arbitrary, macro-iteration-count-dependent global phase -- see
# window_tdvp_step's own docstring for the full story and the fix
# (eshift, restoring consistency with mpscpp3::quench_tdvp's own
# established Hshift=H-EGS convention).
#
# MODEL CHOICE: the *dimerized* (SSH-like) XX chain, not the plain
# gapless uniform one -- still exactly free-fermion-solvable the same
# way, but the dimerization gap makes iDMRG converge far better
# (state_overlap ~0.98-0.999 here vs ~0.87-0.96 for the gapless chain at
# the same bond dimension) -- confirmed directly. The remaining
# real-time discrepancy against the exact answer reflects iDMRG's own
# ground-state convergence quality, not a bug in window_tdvp_step itself
# -- verified separately and exactly by a dense-matrix cross-check (see
# tests/test_idmrg_window_free_fermion.py's own
# test_window_tdvp_step_eshift_matches_exact_dense_evolution).

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')
sys.path.append(os.getcwd()+'/../../../tests')  # for _free_fermion_reference

import numpy as np
from dmrgpy.pyitensor import idmrg, idmrg_window
from _free_fermion_reference import FreeFermionXX

J1, J2 = 1.0, 0.2  # dimerization: J2 << J1 opens a clear gap

hi = [[J1, ["Sx", 0], ["Sx", 1]], [J1, ["Sy", 0], ["Sy", 1]]]
he = [[J2, ["Sx", 1], ["Sx", 2]], [J2, ["Sy", 1], ["Sy", 2]]]

n_window = 8

# iDMRG's own local two-site solve is unseeded (see idmrg.py's own
# _local_two_site_solve docstring) -- run a few attempts and keep the
# best-converged (highest state_overlap) one that also passes
# build_window's own periodic-wraparound self-consistency check.
candidates = []
for seed in range(25):
    np.random.seed(seed)
    result = idmrg.idmrg_ground_state([2, 2], hi, he, 2, maxm=24, cutoff=1e-12,
                                       maxiter=150, etol=1e-14, niter=100, verbose=False)
    candidates.append((result.state_overlap or 0, result))
    candidates.sort(key=lambda p: -p[0])
    if candidates[0][0] > 0.995:
        break

# A high state_overlap alone does not guarantee the periodic-wraparound
# bond dimensions build_window needs actually match (see build_window's
# own RuntimeError message) -- take the best-converged candidate that
# also actually builds a window.
best = None
for _overlap, result in candidates:
    try:
        idmrg_window.build_window(result, n_window)
        best = result
        break
    except RuntimeError:
        continue
if best is None:
    raise RuntimeError("no candidate iDMRG solve passed build_window's own "
                        "self-consistency check -- rerun (unseeded convergence)")

print("iDMRG ground-state energy density:", best.e0)
print("state_overlap (self-consistency diagnostic):", best.state_overlap)

ff = FreeFermionXX(2000, J=J1, J2=J2)
exact_e0 = np.sum(ff.evals[ff.occ]) / ff.N
print("exact free-fermion energy density:", exact_e0)
print()

dt, nt = 0.05, 6
ts, xs, S = idmrg_window.dynamical_correlator_td(
    best, n_window=n_window, opname_A="Sz", opname_B="Sz", dt=dt, nt=nt,
    cutoff=1e-10, maxdim=40, niter=50, x_values=[-1, 0, 1], connected=True)

x0 = 1000  # reference site, far from the free-fermion chain's own open edges
xs = list(xs)
print("connected Sz-Sz correlator S(x,t): iDMRG (window TDVP) vs exact free fermion")
print(f"{'t':>6} {'x':>3}  {'iDMRG':>20}  {'exact':>20}  {'|diff|':>8}")
for it, t in enumerate(ts):
    for ix, x in enumerate(xs):
        exact = ff.sz_sz_connected(x0 + int(x), x0, t)
        got = S[it, ix]
        print(f"{t:6.3f} {x:+3d}  {got.real:+.6f}{got.imag:+.6f}j  "
              f"{exact.real:+.6f}{exact.imag:+.6f}j  {abs(got - exact):.4f}")
