# Real-time IBC-window dynamical correlator (arXiv:1804.09163, Sec. V.1)
# on the uniform S=1/2 Heisenberg chain: pure-Python (pyitensor) backend
# vs the native ITensor v3 (mpscpp3) backend.
#
# THE CHECK THAT MATTERS. S(x, t=0) is not an approximation of anything:
# no evolution has happened yet, so it must equal the chain's own static
# two-point correlator <A_x B_0> exactly, on whatever state the backend
# converged to. That identity needs no exact solution, no second backend,
# and no tolerance argument -- it is the sharpest cheap oracle this
# machinery has. This example evaluates it per backend, against that
# backend's own `correlator`, and then compares the two backends' whole
# S(x,t) trajectories against each other.
#
# A HISTORY NOTE, because it is the transferable lesson. Between
# 2026-08-29 and its fix later the same day, the v3 half of this
# comparison was DISABLED: that backend's window tiled the raw
# per-micro-step `idmrg_U_` factors, whose two ends live in bases minted
# by *different* iDMRG micro-steps, instead of the gauge-consistent unit
# cell that this very backend's own idmrg_onsite_expectation/
# idmrg_two_point_correlator had already been moved onto for exactly that
# reason. The symptom was: exact at x=0, missing by up to ~1e-1 and worse
# with growing |x| -- the signature of a bond-basis (gauge) mismatch, not
# of a physics bug or of convergence. And what let it survive so long was
# this example's own earlier form: it compared the two backends at x=0
# only, with a 0.1 tolerance, and x=0 is the one point a gauge error
# cannot touch. A comparison is only as sharp as the point it is
# evaluated at.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg_window

MAXM, CUTOFF, MAXITER, ETOL, NITER, RESTARTS = 10, 1e-12, 40, 1e-11, 30, 1
N_WINDOW, DT, NT = 6, 0.05, 4
X_VALUES = [-2, -1, 0, 1, 2]


def build(itensor_version):
    # An explicit TWO-site unit cell for what is a uniform chain, on
    # purpose. With n_uc=1 the growing algorithm's extracted cell is still
    # two sites long (idmrg._theta_cell), and at finite bond dimension its
    # two bonds are not exactly equivalent -- so "the" static correlator at
    # an ODD separation depends on which of the two cell positions it
    # starts from, while `correlator`'s own p_i can only ever be 0 there.
    # Comparing the window's own centre (which lands on cell position 1)
    # against a p_i=0 reference then shows a ~3e-2 mismatch that is neither
    # backend's fault -- confirmed directly, it is the period-2 structure,
    # not the window. With n_uc=2 the sublattice is addressable and the
    # identity below is exact at every x.
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * 2,
                                            itensor_version=itensor_version)
    ic.gs_method = "idmrg"  # the window needs the growing algorithm's own
                            # converged-environment snapshot -- not built
                            # by gs_method="vumps" (the default)
    ic.maxm, ic.cutoff, ic.maxiter = MAXM, CUTOFF, MAXITER
    ic.etol, ic.niter, ic.restarts = ETOL, NITER, RESTARTS
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0]
         + ic.SzC[1] * ic.SzR[0])
    ic.set_hamiltonian(h)
    density = ic.gs_energy()
    return ic, density


print("=== ground-state energy density ===")
ic_py, e_py = build("python")
ic_v3, e_v3 = build(3)
print("python density={: .10f}".format(e_py))
print("v3     density={: .10f}".format(e_v3))
print("|diff|={:.2e}   <- the ENERGY is fine on both backends".format(abs(e_py - e_v3)))
assert abs(e_py - e_v3) < 1e-4, "python and v3 iDMRG energy densities disagree"

# -- S(x, t=0) on each backend, against that backend's own static correlator
ts_py, xs_py, S_py = idmrg_window.dynamical_correlator_td(
    ic_py._result, N_WINDOW, "Sz", "Sz", DT, NT,
    cutoff=1e-10, maxdim=20, niter=30, x_values=X_VALUES, connected=False)

# The v3 counterpart of the same S(x,t) array. Both backends' PUBLIC
# td_dynamical_correlator returns S(k,omega) (it Fourier-transforms this
# before handing it back), so the real-space quantity the t=0 identity
# lives in has to come from one level down on each side: the pyitensor
# module function above, and the session method here.
ts_v3, xs_v3, S_v3 = ic_v3._session3.td_dynamical_correlator_window(
    N_WINDOW, "Sz", "Sz", DT, NT, sorted(X_VALUES), 20, 1e-10, 30, False, 0)


def static(ic, x):
    """<Sz_x Sz_0> from the backend's own converged state. Sz is Hermitian
    and parity-even, so no operator-ordering or Jordan-Wigner subtlety
    enters here -- see tests/test_idmrg_window_fermionic.py for the
    fermionic version of the same identity, where both do. What does
    matter is the sublattice: the perturbation sits on position 0, so for
    x<0 the *other* operator is the one at position x%n_uc."""
    x = int(x)
    if x >= 0:
        return complex(ic.correlator("Sz", 0, "Sz", x)).real
    return complex(ic.correlator("Sz", x % 2, "Sz", -x)).real


exact_py = np.array([static(ic_py, x) for x in xs_py])
exact_v3 = np.array([static(ic_v3, x) for x in xs_v3])
err_py = np.abs(S_py[0].real - exact_py)
err_v3 = np.abs(np.array(S_v3)[0].real - exact_v3)

print()
print("=== S(x,t=0) against each backend's own static correlator ===")
print("  x    python S(x,0)   its <SzSz>   |err|      v3 S(x,0)    its <SzSz>   |err|")
for i, x in enumerate(xs_py):
    print("  %+d   %+ .8f   %+ .8f  %.1e   %+ .8f  %+ .8f  %.1e"
          % (x, S_py[0][i].real, exact_py[i], err_py[i],
             np.array(S_v3)[0][i].real, exact_v3[i], err_v3[i]))
print("worst |err|:  python %.2e   v3 %.2e" % (err_py.max(), err_v3.max()))

# Both identities are exact, and both are genuine regression checks.
assert err_py.max() < 1e-9, "the python window broke an exact identity"
assert err_v3.max() < 1e-9, (
    "the v3 window broke the exact S(x,0) == correlator(x) identity -- if "
    "the error is exact at x=0 and grows with |x|, this is the gauge "
    "defect fixed on 2026-08-29 coming back: the window must tile "
    "Chain::idmrg_cell_raw_, not the raw per-micro-step idmrg_U_ factors")

# -- and the two backends against each other, over the whole trajectory
S_py_a = np.array(S_py)
S_v3_a = np.array(S_v3)
print()
print("=== python vs v3, S(x,t) over the whole evolution ===")
print("   t      max_x |S_py - S_v3|")
for it, t in enumerate(ts_py):
    print("  %.2f    %.2e" % (t, np.abs(S_py_a[it] - S_v3_a[it]).max()))
worst = np.abs(S_py_a - S_v3_a).max()
worst_t0 = np.abs(S_py_a[0] - S_v3_a[0]).max()
print("worst at t=0:         %.2e   <- the two backends' own STATES differ"
      % worst_t0)
print("worst over all (x,t): %.2e" % worst)
# The floor here is not the window machinery -- each backend reproduces its
# OWN static correlator to ~1e-16 above. It is that the two iDMRG runs
# converge to slightly different states at this deliberately small MAXM
# (their static correlators already differ by ~1e-4 at t=0), and TDVP
# amplifies that over the evolution. Tighten MAXM/MAXITER and this shrinks;
# a regression in the window itself would blow straight past it.
assert worst < 5e-3, "python and v3 windows disagree on S(x,t)"

# The public API is what a user actually calls, and it returns S(k,omega).
# Check that it runs on v3 (it raised outright while the gauge was broken)
# and lands on the same spectrum as the python backend.
ks_py, es_py, Sk_py = ic_py.td_dynamical_correlator(
    "Sz", 0, "Sz", N_WINDOW, dt=DT, nt=NT, x_values=X_VALUES,
    maxdim=20, cutoff=1e-10, niter=30)
ks_v3, es_v3, Sk_v3 = ic_v3.td_dynamical_correlator(
    "Sz", 0, "Sz", N_WINDOW, dt=DT, nt=NT, x_values=X_VALUES,
    maxdim=20, cutoff=1e-10, niter=30)
scale = np.abs(Sk_py).max()
print()
print("=== public API, S(k,omega): max |python - v3| = %.2e  (peak %.3f) ="
      % (np.abs(Sk_py - Sk_v3).max(), scale))
assert np.abs(Sk_py - Sk_v3).max() < 1e-2 * scale

# -------------------------------------------------------------------- plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

ax1.plot(xs_py, exact_py, "k-", lw=2, label="static correlator (the answer)")
ax1.plot(xs_py, S_py[0].real, "o", color="tab:blue", label="window, python")
ax1.plot(xs_v3, S_v3_a[0].real, "s", mfc="none", color="tab:orange",
         label="window, v3")
ax1.set_xlabel("distance $x$") ; ax1.set_ylabel(r"$S(x,t=0)$")
ax1.set_title("t=0 must reproduce the static correlator exactly")
ax1.legend(fontsize=8)

for i, x in enumerate(xs_py):
    line, = ax2.plot(ts_py, S_py_a[:, i].real, "-", label="$x=%+d$" % x)
    ax2.plot(ts_v3, S_v3_a[:, i].real, "o", mfc="none",
             color=line.get_color())
ax2.set_xlabel("time $t$") ; ax2.set_ylabel(r"$\mathrm{Re}\,S(x,t)$")
ax2.set_title("lines: python    circles: v3")
ax2.legend(fontsize=8, ncol=2)

fig.suptitle("iDMRG IBC-window dynamical correlator, python vs ITensor v3")
fig.tight_layout()
plt.show()
