# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time

import numpy as np  # conventional numpy library
import matplotlib.pyplot as plt  # plotting
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object
from dmrgpy.pyitensor import idmrg_excitations  # the excitation ansatz itself

###########################################################################
### How the tangent-space excitation ansatz's eigenproblem is solved:   ###
### assembling H_eff(k) vs. Lanczos on its action, and where the two    ###
### cross over (pyitensor/idmrg_excitations.py)                         ###
###########################################################################
# excitation_energies solves H_eff(k)[X] = lambda*X, an eigenproblem of
# dimension dim = D*D*(d_g-1). Two solvers are available and are picked by
# size (idmrg_excitations._DENSE_EIG_MAX):
#
#   dense      -- build the matrix one basis vector at a time, then eigh.
#                 Costs `dim` applications of _h_eff_action.
#   iterative  -- Lanczos (scipy eigsh) straight on _h_eff_action, never
#                 assembling anything. Costs a few hundred applications
#                 regardless of `dim`.
#
# So the dense path wins at small `dim` and loses at large `dim`, and this
# script measures where the crossing actually is. It also checks the two
# agree, which is the point: the fast path must not cost accuracy.
#
# Both are exercised here by monkeypatching the threshold, so the whole
# sweep stays cheap enough to run by hand; in ordinary use one or the other
# is selected automatically and the caller never sees this choice.
#
# H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i] = -sigma^x_i sigma^x_{i+1} - g*sigma^z_i
# (the same transverse-field Ising convention as excitation_gap_tfim).

J, g = 1.0, 2.5  # g > 1: gapped paramagnetic phase
ks = np.linspace(0, np.pi, 7)
bond_dims = [2, 3, 4, 6, 8]


def ground_state(D):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "vumps"  # required by the excitation ansatz
    ic.maxm = D
    ic.etol = 1e-12
    ic.vumps_nrestarts = 4
    ic.set_hamiltonian(-4.0*J*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
    ic.gs_energy()
    return ic


def scan(env, threshold):
    """(wall time, dispersion) for a full momentum scan with the solver
    threshold forced one way or the other."""
    original = idmrg_excitations._DENSE_EIG_MAX
    idmrg_excitations._DENSE_EIG_MAX = threshold
    env.resolvent_cache.clear()  # time the resolvent builds too, both ways
    try:
        t0 = time.time()
        band = np.array([idmrg_excitations.excitation_energies(env, k, n=1)[0].real
                          for k in ks])
        return time.time() - t0, band
    finally:
        idmrg_excitations._DENSE_EIG_MAX = original


dims, t_dense, t_iter, disagreement = [], [], [], []
exact = 2*np.sqrt(J**2 + g**2 - 2*J*g*np.cos(ks))

for D in bond_dims:
    ic = ground_state(D)
    env = ic._get_excitation_environment()
    dim = env.D*env.D*(env.d_g - 1)
    td, band_d = scan(env, 10**9)  # force dense
    ti, band_i = scan(env, 0)      # force iterative
    dims.append(dim)
    t_dense.append(td)
    t_iter.append(ti)
    disagreement.append(np.max(np.abs(band_d - band_i)))
    print("D={:2d}  dim={:4d}  dense {:7.3f}s   iterative {:7.3f}s   "
          "speedup {:5.2f}x   max|dense-iter| {:.2e}   |E-exact| {:.2e}".format(
              D, dim, td, ti, td/max(ti, 1e-9), disagreement[-1],
              np.max(np.abs(band_d - exact))))

# The two solvers must agree to well below any tolerance the dispersion is
# ever asserted at -- that is what makes the fast path usable at all.
assert max(disagreement) < 1e-8, disagreement

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

ax.plot(dims, t_dense, "o-", color="tab:blue", label="dense (assemble + eigh)")
ax.plot(dims, t_iter, "s-", color="tab:red", label="iterative (Lanczos on the action)")
ax.axvline(idmrg_excitations._DENSE_EIG_MAX, color="black", ls="--", lw=1,
           label="_DENSE_EIG_MAX (the automatic choice)")
ax.set_xlabel(r"eigenproblem dimension $D^2(d_g-1)$")
ax.set_ylabel("wall time for a {}-point $k$ scan [s]".format(len(ks)))
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_title("Cost of the two solvers")
ax.legend(fontsize=8)

ax2.plot(dims, disagreement, "o-", color="tab:green")
ax2.set_xlabel(r"eigenproblem dimension $D^2(d_g-1)$")
ax2.set_ylabel(r"$\max_k |E_{\rm dense}(k) - E_{\rm iterative}(k)|$")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_title("The fast path costs no accuracy")

fig.tight_layout()
fig.savefig("excitation_solver_scaling.png", dpi=150)
print("\nsaved plot to excitation_solver_scaling.png")
