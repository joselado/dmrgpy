# Real-time IBC-window dynamical correlator (arXiv:1804.09163, Sec. V.1):
# pure-Python (pyitensor) backend vs the native ITensor v3 (mpscpp3)
# backend, on the uniform S=1/2 Heisenberg chain (n_uc=1) --
# infinitechain.Infinite_Many_Body_Chain.td_dynamical_correlator with
# itensor_version="python" vs itensor_version=3.
#
# Mirrors heisenberg_infinite_python_VS_v3/main.py's own python-vs-v3
# idmrg_ground_state comparison, extended to the time-dependent
# correlator: Chain::td_dynamical_correlator_window (mpscpp3/
# chain_session.h) reuses the *native* ITensorTDVP boundary-tensor tdvp()
# overload (TDVP/tdvp.h's tdvp(psi,H,t,LH,RH,sweeps,args), backed by
# LocalMPO(H,LH,RH,args)) directly against a tiled window MPS/MPO, rather
# than pyitensor/idmrg_window.py's own hand-rolled window-aware TDVP sweep
# (needed there only because pyitensor's generic tdvp.py infers a site's
# Link via a same-Index chain-neighbor lookup that cannot see a window's
# extra boundary legs -- see that module's own docstring).
#
# Model choice: pure isotropic Heisenberg (no onsite field term) --
# deliberately *not* the transverse-field model this same directory's
# own td_dynamical_correlator/main.py uses: confirmed directly (on an
# unmodified checkout, unrelated to the window/TDVP code here) that
# Chain::idmrg_ground_state's onsite ("field") term handling diverges on
# the v3 backend for that model. A known, pre-existing, orthogonal v3
# limitation -- flagged as a follow-up, not fixed here.
#
# Both backends solve their local two-site iDMRG problem from an unseeded
# random start, so two independent runs converge to independently-gauged
# realizations of the same physical ground state, not bit-identical ones
# -- the comparison below uses a loose tolerance for the same reason
# test_idmrg_window_v3.py's own cross-backend checks do.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import infinitechain

MAXM, CUTOFF, MAXITER, ETOL, NITER, RESTARTS = 10, 1e-12, 40, 1e-11, 30, 1
N_WINDOW, DT, NT = 6, 0.05, 4
X_VALUES = [-1, 0, 1]


def build(itensor_version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "idmrg"  # td_dynamical_correlator needs the growing
                            # algorithm's own converged-environment
                            # snapshot -- not built by gs_method="vumps"
                            # (the default since 2026-08-08)
    ic.maxm, ic.cutoff, ic.maxiter = MAXM, CUTOFF, MAXITER
    ic.etol, ic.niter, ic.restarts = ETOL, NITER, RESTARTS
    h = ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]
    ic.set_hamiltonian(h)
    density = ic.gs_energy()
    return ic, density


print("=== ground-state energy density ===")
ic_py, e_py = build("python")
ic_v3, e_v3 = build(3)
print("python density={: .10f}".format(e_py))
print("v3     density={: .10f}".format(e_v3))
print("|diff|={:.2e}".format(abs(e_py - e_v3)))
assert abs(e_py - e_v3) < 1e-4, "python and v3 iDMRG energy densities disagree"

print()
print("=== S(x,t) via IBC-window TDVP ===")
# td_dynamical_correlator's own public API returns (ks,es,Skw) -- reach
# the raw S(x,t) via the same lower-level functions each backend's own
# dispatch calls instead, so this example compares the
# *pre-Fourier-transform* correlator (the quantity
# td_dynamical_correlator_window/dynamical_correlator_td actually
# compute) rather than the doubly-transformed S(k,omega), which would
# additionally fold in the (independent, backend-agnostic)
# Lorentzian-broadening FFT convention.
from dmrgpy.pyitensor import idmrg_window
ts_py, xs_py, S_py = idmrg_window.dynamical_correlator_td(
    ic_py._result, N_WINDOW, "Sz", "Sz", DT, NT,
    cutoff=1e-10, maxdim=20, niter=30, x_values=X_VALUES, connected=True)
ts_v3, xs_v3, S_v3 = ic_v3._session3.td_dynamical_correlator_window(
    N_WINDOW, "Sz", "Sz", DT, NT, sorted(X_VALUES), 20, 1e-10, 30, True, 0)

print("python S(x,t=0)  :", np.round(S_py[0].real, 4))
print("v3     S(x,t=0)  :", np.round(S_v3[0].real, 4))
ix0 = list(xs_py).index(0)
print("|diff| at x=0,t=0: {:.4f}".format(abs(S_py[0][ix0] - S_v3[0][ix0])))
assert abs(S_py[0][ix0] - S_v3[0][ix0]) < 0.1
assert np.all(np.isfinite(S_v3))

print()
print("=== S(k,omega) via the public API (itensor_version=3) ===")
ks = np.linspace(-np.pi, np.pi, 5)
ks_out, es, Skw = ic_v3.td_dynamical_correlator(
    "Sz", 0, "Sz", n_window=N_WINDOW, dt=DT, nt=6, x_values=X_VALUES,
    maxdim=20, cutoff=1e-10, niter=30, ks=ks, delta=0.2)
assert np.all(np.isfinite(Skw))
print("max|S(k,omega)| (v3) =", np.max(np.abs(Skw)))

print()
print("iDMRG td_dynamical_correlator python-VS-v3 regression test PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

xs_py_list = list(xs_py)
xs_v3_list = list(xs_v3)
ax1.plot(xs_py_list, S_py[0].real, "o-", color="tab:blue", label="python")
ax1.plot(xs_v3_list, S_v3[0].real, "s--", color="tab:orange", label="v3")
ax1.set_xlabel("distance $x$")
ax1.set_ylabel(r"Re $S(x,t=0)$")
ax1.set_title("S(x,t=0): python vs v3")
ax1.legend()

pcm = ax2.pcolormesh(ks_out, es, np.abs(Skw).T, shading="auto", cmap="viridis")
fig.colorbar(pcm, ax=ax2, label=r"|S(k,$\omega$)|")
ax2.set_xlabel("k")
ax2.set_ylabel(r"$\omega$")
ax2.set_title("S(k,$\\omega$) (v3)")

fig.suptitle("iDMRG IBC-window TDVP dynamical correlator: python vs ITensor v3")
fig.tight_layout()
fig.savefig("td_dynamical_correlator_python_VS_v3.png", dpi=150)
print("saved plot to td_dynamical_correlator_python_VS_v3.png")
