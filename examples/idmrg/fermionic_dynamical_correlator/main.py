# Real-time dynamical correlator of a FERMIONIC infinite chain, i.e. one
# whose operators drag a Jordan-Wigner string:
#
#     S(x,t) = <Cdag_x(t) C_0>,   and its Fourier transform S(k,omega),
#
# computed directly in the thermodynamic limit with the infinite-boundary-
# condition (IBC) window method (Infinite_Many_Body_Chain.
# td_dynamical_correlator -> pyitensor/idmrg_window.py), and compared
# against the EXACT free-fermion answer at every (x,t).
#
# WHY A FERMIONIC PAIR IS DIFFERENT. Under Jordan-Wigner the physical
# operator at site s is (prod_{j<s} F_j) O_s -- a semi-infinite object,
# not a single-site matrix. Both the perturbation applied to the window
# (the ket, before the evolution) and the measurement (the bra, at every
# time step) have to carry that string, and the two truncations have to
# agree; get it wrong and S(x,t) is not "less converged" but a different
# number entirely, with wrong signs. Until 2026-08-29 this path applied a
# bare C/Cdag matrix with no string at all, on both backends -- which is
# why this example exists.
#
# THE ORACLE. For a quadratic H = sum_ij h_ij c^dag_i c_j the exact answer
# is one N x N diagonalization, at any N and any time:
#
#     <c^dag_x(t) c_0> = sum_l [e^{+i h t}]_{xl} P[l,0],
#     P[l,m] = <c^dag_l c_m> (all negative single-particle levels filled),
#
# so this is a quantitative check against an external number, not a
# self-consistency check. The same exact S(x,t) is pushed through the
# *same* Fourier transform as the dmrgpy data for the S(k,omega) panel,
# so the two maps are directly comparable with no band-folding
# bookkeeping.
#
# NOTE ON RESOLUTION: the two S(k,omega) maps below are broadened by the
# simulated time alone (~1/(nt*dt) ~ 0.7 here), which is why the occupied
# band shows up as two blobs rather than a sharp dispersion. Raise nt to
# sharpen it -- both panels equally, since the exact reference goes
# through the very same transform, which is the point of showing them side
# by side rather than overlaying an analytic band.
#
# WHAT YOU SHOULD SEE (measured at the parameters below, maxm=20,
# n_window=8, t up to 1.5): the energy density lands on the exact band
# integral to 1e-10, S(x,0) on the exact one-body density matrix to ~1e-8,
# and the worst |S - S_exact| over the whole (x,t) grid is ~8e-2, growing
# smoothly with t -- that growth is the window's own finite size and bond
# dimension, the honest error of the method. The two S(k,omega) maps
# should look like the same object at slightly different resolution.
#
# MODEL: the dimerized (SSH-like) spinless chain, two sites per cell, with
# an extra next-cell A-A hopping t3. That t3 term is what makes the string
# load-bearing already in the Hamiltonian (its two endpoints have a site
# strictly between them), and the dimerization gap is what makes iDMRG
# converge well -- the uniform half-filled chain is critical and is this
# codebase's documented iDMRG convergence trap.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt

from dmrgpy import infinitechain, timedependent

FERMION_SITE = 0            # spinless fermion site: C / Cdag / F / N
t1, t2, t3 = 1.0, 0.4, 0.1  # intra-cell, inter-cell, next-cell A-A


# ---------------------------------------------------------------- exact --
def single_particle_matrix(L=400):
    """h for a ring of L cells (site index 2*cell + sublattice)."""
    N = 2*L
    h = np.zeros((N, N))
    for n in range(L):
        a, b, a2 = 2*n, 2*n+1, (2*n+2) % N
        h[a, b] += t1 ; h[b, a] += t1
        h[b, a2] += t2 ; h[a2, b] += t2
        h[a, a2] += t3 ; h[a2, a] += t3
    return h


h = single_particle_matrix()
w, v = np.linalg.eigh(h)
P = (v[:, w < 0].conj() @ v[:, w < 0].T).real
i0 = 2*(h.shape[0]//4)      # a bulk A site


# ------------------------------------------------------------- the chain --
ic = infinitechain.Infinite_Many_Body_Chain([FERMION_SITE]*2,
                                             itensor_version="python")
ic.gs_method = "idmrg"      # td_dynamical_correlator needs the growing
ic.maxm = 20                # algorithm's own converged environment
ic.maxiter = 40

C = [ic.get_operator("C", i, "C") for i in range(2)]
Cd = [ic.get_operator("Cdag", i, "C") for i in range(2)]
CR = [ic.get_operator("C", i, "R") for i in range(2)]
CdR = [ic.get_operator("Cdag", i, "R") for i in range(2)]
ic.set_hamiltonian(t1*(Cd[0]*C[1] + Cd[1]*C[0])
                   + t2*(Cd[1]*CR[0] + CdR[0]*C[1])
                   + t3*(Cd[0]*CR[0] + CdR[0]*C[0]))

e0 = ic.gs_energy()
print("energy density: %.10f  (exact %.10f)"
      % (e0, w[w < 0].sum()/h.shape[0]))

# ------------------------------------------------------- S(x,t), both ways --
n_window, dt, nt = 8, 0.1, 16
xs = list(range(-5, 6))

from dmrgpy.pyitensor import idmrg_window
ts, xarr, S = idmrg_window.dynamical_correlator_td(
    ic._result, n_window, "Cdag", "C", dt=dt, nt=nt, cutoff=1e-10,
    maxdim=60, x_values=xs, connected=False)

S_exact = np.zeros_like(S)
for it, t in enumerate(ts):
    U = scipy.linalg.expm(1j*h*t)
    for ix, x in enumerate(xarr):
        S_exact[it, ix] = U[i0 + int(x), :] @ P[:, i0]

err = np.abs(S - S_exact)
print("worst |S - S_exact| over the whole (x,t) grid: %.3e" % err.max())

# ------------------------------------------------------------- S(k,omega) --
es = np.linspace(-3.5, 3.5, 400)
ks = np.linspace(-np.pi, np.pi, 201)
kw = dict(ks=ks, es=es, delta=0.15, window=[-1, 10], factor=1)
ks_o, es_o, Skw = timedependent.sxt_to_skomega(ts, xarr, S, dt, **kw)
_ks, _es, Skw_exact = timedependent.sxt_to_skomega(ts, xarr, S_exact, dt, **kw)

# ------------------------------------------------------------------ plot --
fig = plt.figure(figsize=(13, 8))

ax = fig.add_subplot(2, 2, 1)
for x, style in ((0, "-"), (1, "--"), (3, ":")):
    ix = list(xarr).index(x)
    ax.plot(ts, S[:, ix].real, style, color="C0", lw=2,
            label="dmrgpy, x=%d" % x)
    ax.plot(ts, S_exact[:, ix].real, style, color="C3", lw=1,
            label="exact, x=%d" % x)
ax.set_xlabel("t") ; ax.set_ylabel(r"Re $S(x,t)$")
ax.set_title(r"$S(x,t)=\langle C^\dagger_x(t)\, C_0\rangle$ (blue) vs exact (red)")
ax.legend(fontsize=7, ncol=2)

ax = fig.add_subplot(2, 2, 2)
for x in (0, 1, 3, 5):
    ix = list(xarr).index(x)
    ax.semilogy(ts[1:], np.maximum(err[1:, ix], 1e-16), label="x=%d" % x)
ax.set_xlabel("t") ; ax.set_ylabel(r"$|S-S_{\rm exact}|$")
ax.set_title("error: exact at t=0 by construction,\ngrowing with the window's own truncation")
ax.legend(fontsize=8)

vmax = np.abs(Skw_exact).max()
for i, (data, name) in enumerate(((Skw, "dmrgpy (IBC window)"),
                                   (Skw_exact, "exact free fermions"))):
    ax = fig.add_subplot(2, 2, 3+i)
    ax.imshow(np.abs(data).T, origin="lower", aspect="auto", vmin=0, vmax=vmax,
              extent=[ks_o[0], ks_o[-1], es_o[0], es_o[-1]], cmap="inferno")
    ax.set_xlabel("k") ; ax.set_ylabel(r"$\omega$")
    ax.set_title(r"$|S(k,\omega)|$ -- %s" % name)

fig.suptitle("Fermionic dynamical correlator of an infinite chain "
              "(Jordan-Wigner string threaded through the IBC window)")
fig.tight_layout()
plt.show()
