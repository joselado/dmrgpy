# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time
import numpy as np  # conventional numpy library
import matplotlib.pyplot as plt
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object

#####################################################################
### Couplings reaching further than one unit cell                 ###
#####################################################################
# An infinite chain's operators come as SxC[i] (this cell), SxR[i] (the
# next one) and SxL[i] (the previous one), so those three lists reach
# exactly one cell either way. For anything longer,
# get_operator(name, i, group=c) takes an INTEGER cell offset -- c=-1,0,1
# are exactly L, C and R, and c=2 is the cell after the next one -- which
# is how the next-nearest-neighbour term below is written on a ONE-site
# unit cell.
#
# The alternative, and the only route available before, is to rewrite the
# same chain on a 2-site cell, where the same coupling becomes an
# ordinary C-to-R bond. Both are computed here and must agree: they are
# the same physical model, solved by two different code paths
# (gs_method="vumps" routes the reach-2 version to the sequential
# multi-site solver, pyitensor/vumps_ms.py, and the reach-1 version to the
# grouped one). The cost is what differs, and grows with the range: the
# sequential solver pays one extra automaton channel per site of reach --
# linear -- while making a range-R coupling reach-1 by hand costs the
# grouped solver a d**R supersite.

G = 2.5          # transverse field, deep in the paramagnetic phase
J1 = 1.0         # nearest-neighbour Ising coupling
D = 8            # bond dimension


def ising_j2(n_uc, J2):
    """H = -G sum sigma^z - J1 sum sigma^x_i sigma^x_{i+1}
           - J2 sum sigma^x_i sigma^x_{i+2}, on an n_uc-site cell."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc)
    ic.maxm, ic.maxiter, ic.etol = D, 300, 1e-10
    ic.vumps_nrestarts = 2
    h = 0
    for i in range(n_uc):
        h = h - 2.0 * G * ic.SzC[i]                     # sigma^z = 2*Sz
        for r, J in ((1, J1), (2, J2)):
            k = i + r
            # site i+r of this cell if it fits, otherwise site (i+r)%n_uc
            # of the cell (i+r)//n_uc steps to the right
            other = (ic.SxC[k] if k < n_uc
                     else ic.get_operator("Sx", k % n_uc, group=k // n_uc))
            h = h - 4.0 * J * ic.SxC[i] * other          # sigma^x = 2*Sx
    ic.set_hamiltonian(h)
    return ic


j2s = np.linspace(0.0, 1.0, 9)
energies, times = {1: [], 2: []}, {1: [], 2: []}
for J2 in j2s:
    for n_uc in (1, 2):
        t0 = time.time()
        ic = ising_j2(n_uc, J2)
        energies[n_uc].append(ic.gs_energy())
        times[n_uc].append(time.time() - t0)
    print("J2 = %5.3f   n_uc=1 (reach 2) e = %.10f   n_uc=2 (reach 1) e = %.10f"
          "   difference %.2e"
          % (J2, energies[1][-1], energies[2][-1],
             abs(energies[1][-1] - energies[2][-1])))

def tight_ylim(ax, values, margin=0.03, log=False):
    """Fit the y-axis to the data, so the curve reaches close to the top
    of the frame instead of floating in the middle of a padded range."""
    lo, hi = float(np.min(values)), float(np.max(values))
    if log:
        ax.set_ylim(lo / (1.0 + 10 * margin), hi * (1.0 + 10 * margin))
    else:
        pad = margin * max(hi - lo, abs(hi) * 1e-12, 1e-12)
        ax.set_ylim(lo - pad, hi + pad)


fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(14, 4))

ax0.plot(j2s, energies[1], "o-", label="$n_{uc}=1$, reach 2 (sequential)")
ax0.plot(j2s, energies[2], "s--", label="$n_{uc}=2$, reach 1 (grouped)")
ax0.set_xlabel("$J_2$")
ax0.set_ylabel("energy density")
ax0.set_title("Long-range Ising chain, $D=%d$" % D)
ax0.legend()
tight_ylim(ax0, energies[1] + energies[2])

diffs = np.abs(np.array(energies[1]) - np.array(energies[2])) + 1e-18
ax1.semilogy(j2s, diffs, "o-", color="C3")
ax1.set_xlabel("$J_2$")
ax1.set_ylabel("|difference| between the two cells")
ax1.set_title("Same model, two code paths")
tight_ylim(ax1, diffs, log=True)

ax2.plot(j2s, times[1], "o-", label="$n_{uc}=1$ (sequential)")
ax2.plot(j2s, times[2], "s--", label="$n_{uc}=2$ (grouped)")
ax2.set_xlabel("$J_2$")
ax2.set_ylabel("wall time (s)")
ax2.set_title("Cost of each route")
ax2.legend()
tight_ylim(ax2, times[1] + times[2])

plt.tight_layout()
plt.show()
