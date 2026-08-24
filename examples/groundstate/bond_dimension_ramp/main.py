# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground-state bond-dimension ramp (Many_Body_Chain.bond_ramp).
#
# Instead of running every DMRG sweep at the full target bond dimension
# self.maxm, the sweep schedule now spends the first
# self.bond_ramp_fraction of the sweeps growing the bond dimension
# geometrically from self.bond_ramp_start up to self.maxm, and holds it
# there for the rest (at the default fraction 0.5 the second half -- and
# hence the returned energy -- always runs at the full self.maxm). Since
# two-site DMRG costs ~O(m^3) per bond, those early sweeps -- which mostly
# just locate the right variational subspace -- become much cheaper, and
# the expensive full-maxm sweeps start from an already-good state. The
# noise term decays by self.bond_ramp_noise_decay per ramping sweep and is
# off entirely once the schedule reaches the full self.maxm.
#
# The model here is a 30-site *spatially inhomogeneous* Heisenberg-Hubbard
# chain: spinful fermions with a site-dependent hopping t_i, a
# site-dependent Hubbard U_i and a site-dependent Heisenberg exchange J_i
# between neighboring spins, all incommensurately modulated (periods 7, 11
# and 5 sites) so that no two bonds are equivalent and the ground state
# has genuine spatial structure. 30 spinful sites is 4^30 states, far
# beyond ED -- so correctness is checked the way it has to be at this
# size: the ramped and un-ramped runs must agree on the energy, and on the
# site-resolved density and magnetization profiles.
#
# The script sweeps the target bond dimension, times ramp vs. flat at each
# one, and plots both the speedup and the energy agreement.
import time

import numpy as np

from dmrgpy import fermionchain

n = 30  # number of spinful fermionic sites (4^30 states -- DMRG only)


def build(maxm, nsweeps, ramp):
    """A 30-site Heisenberg-Hubbard chain with spatial variations."""
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        t = 1.0 + 0.3 * np.cos(2 * np.pi * i / 7.0)   # modulated hopping
        J = 0.4 + 0.3 * np.sin(2 * np.pi * i / 5.0)   # modulated exchange
        h = h + t * fc.Cdagup[i] * fc.Cup[i + 1]
        h = h + t * fc.Cdagdn[i] * fc.Cdn[i + 1]
        h = h + J * (fc.Sx[i] * fc.Sx[i + 1] + fc.Sy[i] * fc.Sy[i + 1]
                     + fc.Sz[i] * fc.Sz[i + 1])
    for i in range(n):
        U = 2.0 + 1.5 * np.cos(2 * np.pi * i / 11.0)  # modulated Hubbard U
        h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
    h = h + h.get_dagger()
    h = 0.5 * h  # undo the double counting from the h.c. above
    fc.maxm = maxm
    fc.nsweeps = nsweeps
    fc.bond_ramp = ramp
    fc.set_hamiltonian(h)
    return fc


nsweeps = 16
maxms = [30, 60, 90]

t_flat, t_ramp, e_flat, e_ramp = [], [], [], []
chains = {}
for maxm in maxms:
    out = {}
    for ramp in (False, True):
        fc = build(maxm, nsweeps, ramp)
        t0 = time.time()
        e = fc.gs_energy()
        out[ramp] = (time.time() - t0, e)
        chains[ramp] = fc  # keep the last size's chains for the profiles
        print("maxm=%3d  %s  E = %.10f   t = %6.1f s"
              % (maxm, "ramp" if ramp else "flat", e, out[ramp][0]),
              flush=True)
    t_flat.append(out[False][0]) ; e_flat.append(out[False][1])
    t_ramp.append(out[True][0])  ; e_ramp.append(out[True][1])
    print("          speedup = %.2fx   dE = %.2e"
          % (t_flat[-1] / t_ramp[-1], abs(e_ramp[-1] - e_flat[-1])))

# The ramp is a scheduling change, not a different variational problem:
# both runs must land on the same ground state. Check that on the
# site-resolved observables too, not only on the total energy (the
# largest-maxm pair, still converged and in hand from the loop above).
fc_flat, fc_ramp = chains[False], chains[True]
d_flat = np.array(fc_flat.get_density()).real
d_ramp = np.array(fc_ramp.get_density()).real
# (mx, my, mz) of every orbital; column 2 is <S^z_i>.
m_flat = np.array(fc_flat.get_magnetization()).real
m_ramp = np.array(fc_ramp.get_magnetization()).real
print("max |dn_i| =", np.max(np.abs(d_ramp - d_flat)))
print("max |dm_i| =", np.max(np.abs(m_ramp - m_flat)))
# Regression guards. The energy is variational, so it agrees far more
# tightly than the local observables do -- a density profile at a finite
# maxm on a 30-site chain is the looser of the two checks by design.
assert np.max(np.abs(np.array(e_ramp) - np.array(e_flat))) < 1e-5
assert np.max(np.abs(d_ramp - d_flat)) < 1e-5

import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(14, 4))

ax[0].plot(maxms, t_flat, marker="o", label="flat schedule")
ax[0].plot(maxms, t_ramp, marker="s", label="bond-dimension ramp")
ax[0].set_xlabel("target bond dimension maxm")
ax[0].set_ylabel("ground-state time [s]")
ax[0].set_title("n=%d inhomogeneous Heisenberg-Hubbard" % n)
ax[0].legend()

ax[1].plot(maxms, np.array(t_flat) / np.array(t_ramp), marker="o",
           color="C2")
ax[1].axhline(1.0, color="k", ls="--", lw=0.8)
ax[1].set_xlabel("target bond dimension maxm")
ax[1].set_ylabel("speedup (flat / ramp)")
ax[1].set_title("speedup from the ramp")

ax[2].plot(maxms, e_flat, marker="o", label="flat schedule")
ax[2].plot(maxms, e_ramp, marker="s", ls="--", label="ramp")
ax[2].set_xlabel("target bond dimension maxm")
ax[2].set_ylabel("ground-state energy")
ax[2].set_title("same energy, less time")
ax[2].legend()

plt.tight_layout()
plt.show()

# ...and the spatial structure the ramp must reproduce exactly.
fig2, ax2 = plt.subplots(1, 2, figsize=(10, 4))
ax2[0].plot(range(n), d_flat, marker="o", label="flat")
ax2[0].plot(range(n), d_ramp, marker="s", ls="--", label="ramp")
ax2[0].set_xlabel("site") ; ax2[0].set_ylabel("density $n_i$")
ax2[0].legend()
ax2[1].plot(range(n), m_flat[:, 2], marker="o", label="flat")
ax2[1].plot(range(n), m_ramp[:, 2], marker="s", ls="--", label="ramp")
ax2[1].set_xlabel("site") ; ax2[1].set_ylabel("magnetization $\\langle S^z_i\\rangle$")
ax2[1].legend()
plt.tight_layout()
plt.show()
