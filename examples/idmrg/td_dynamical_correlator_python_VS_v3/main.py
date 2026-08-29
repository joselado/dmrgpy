# Real-time IBC-window dynamical correlator (arXiv:1804.09163, Sec. V.1)
# on the uniform S=1/2 Heisenberg chain: pure-Python (pyitensor) backend
# vs the native ITensor v3 (mpscpp3) backend -- and, since 2026-08-29, a
# demonstration of WHY the v3 half of that comparison is currently
# disabled rather than a demonstration that the two agree.
#
# THE CHECK THAT MATTERS. S(x, t=0) is not an approximation of anything:
# no evolution has happened yet, so it must equal the chain's own static
# two-point correlator <A_x B_0> exactly, on whatever state the backend
# converged to. That identity needs no exact solution, no second backend,
# and no tolerance argument -- it is the sharpest cheap oracle this
# machinery has. This example evaluates it per backend, against that
# backend's own `correlator`:
#
#   itensor_version="python": exact to ~1e-15 at every x.
#   itensor_version=3:        misses by up to ~1e-1, growing with |x|,
#                             while staying exact at x=0.
#
# That last pattern -- exact at x=0, worse the further apart the two
# operators sit -- is the signature of a bond-basis (gauge) mismatch, not
# of a physics bug or of convergence: Chain::idmrg_build_window and
# Chain::idmrg_window_snapshot_correlator tile the raw per-micro-step
# `idmrg_U_` factors, whose two ends live in bases minted by *different*
# iDMRG micro-steps, instead of the gauge-consistent unit cell
# (idmrg_theta_cell/ic_canonicalize_cell) that this very backend's own
# idmrg_onsite_expectation/idmrg_two_point_correlator were moved onto for
# exactly this reason. See docs/known_issue_v3_window_gauge.md.
#
# WHAT THIS EXAMPLE USED TO DO, AND WHY IT MISSED IT: it compared the two
# backends' S(x,t=0) at x=0 only, with a 0.1 tolerance -- and x=0 is the
# one point the gauge error cannot touch. A comparison is only as sharp as
# the point it is evaluated at; an identity evaluated everywhere would
# have caught this years earlier. That is the transferable lesson here.
#
# Because the public API now raises for itensor_version=3 (a wrong number
# is worse than an error), the v3 side below is reached through the
# deliberate, tests-only opt-in Chain::set_allow_defective_window -- the
# same one tests/test_idmrg_window_fermionic.py's strict xfail uses.
# When the gauge is fixed, delete the opt-in, the raise, and this
# paragraph, and this example becomes an ordinary agreement check again.

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

# the v3 window refuses to run unless the caller acknowledges the defect
ic_v3._session3.set_allow_defective_window(True)
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

# The python backend's identity is exact and is a genuine regression check.
assert err_py.max() < 1e-9, "the python window broke an exact identity"
# The v3 backend's is not -- this asserts the defect is still what the
# known-issue document says it is, so this example fails loudly if the
# situation changes in either direction (fixed, or worse).
assert err_v3.max() > 1e-3, (
    "v3's S(x,0) now matches its own static correlator -- if the window "
    "gauge was fixed, delete the opt-in above, the raise in "
    "infinitechain.td_dynamical_correlator, the C++ refusal, and "
    "docs/known_issue_v3_window_gauge.md, and turn this example back into "
    "an agreement check")

print()
print("=== the public API refuses itensor_version=3 ===")
try:
    ic_v3.td_dynamical_correlator("Sz", 0, "Sz", n_window=N_WINDOW, dt=DT, nt=4)
    raise AssertionError("expected td_dynamical_correlator to refuse on v3")
except RuntimeError as e:
    print(str(e))

# -------------------------------------------------------------------- plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

ax1.plot(xs_py, exact_py, "k-", lw=2, label="static correlator (the answer)")
ax1.plot(xs_py, S_py[0].real, "o", color="tab:blue", label="window, python")
ax1.plot(xs_v3, np.array(S_v3)[0].real, "s", color="tab:orange",
         label="window, v3 (gauge-defective)")
ax1.set_xlabel("distance $x$") ; ax1.set_ylabel(r"$S(x,t=0)$")
ax1.set_title("t=0 must reproduce the static correlator exactly")
ax1.legend(fontsize=8)

ax2.semilogy(xs_py, np.maximum(err_py, 1e-17), "o-", color="tab:blue",
             label="python (~1e-15, i.e. exact)")
ax2.semilogy(xs_v3, np.maximum(err_v3, 1e-17), "s-", color="tab:orange",
             label="v3 (grows with |x|, exact at x=0)")
ax2.set_xlabel("distance $x$") ; ax2.set_ylabel("|window - static|")
ax2.set_title("the same identity, as an error")
ax2.legend(fontsize=8)

fig.suptitle("iDMRG IBC-window dynamical correlator: the t=0 identity, "
              "python vs ITensor v3")
fig.tight_layout()
plt.show()
