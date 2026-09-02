# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Conserved-sector targeting on the ED backend (mode="ED").
#
# Many_Body_Chain.set_conserved_sector(Nf=k) / (Sz=m) confines every
# calculation on a chain to one quantum-number sector. The DMRG backends
# do that structurally (itensor_version=3's block-sparse QN tensors) or
# variationally (itensor_version="python"'s charge penalty); ED does it
# the simplest way of the three -- a conserved charge is diagonal in the
# ED product basis, so a sector is a *set of basis states* and confining
# a calculation to it means restricting every assembled operator to the
# corresponding submatrix (see edtk/edchain.py).
#
# What that buys is the same thing it buys on DMRG: E(N) as a family of
# genuine ground-state energies, one per sector, instead of something
# extracted by tuning a chemical potential. Here that is the addition
# spectrum of a spinless t-V chain and the charge gap
#
#     Delta(N) = E(N+1) + E(N-1) - 2 E(N)
#
# read off from it. The second panel does the same for the total-Sz
# sectors of a Heisenberg chain, whose E(Sz) curve is the magnetization
# curve in disguise: its Legendre transform in a field B is what
# minimizes E(Sz) - B*Sz.
#
# Every sector energy is checked against an independent reference built
# by hand -- the full-Hilbert-space Hamiltonian restricted to the basis
# states carrying the requested charge -- so this doubles as a regression
# test in the spirit of tests/.

import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import fermionchain, spinchain
from dmrgpy.multioperator import MO2matrix


def sector_energy_by_hand(chain, h, charge, target):
    """Reference: assemble H on the *full* Hilbert space, then take the
    submatrix of the basis states with the requested charge. This never
    touches the sector machinery being checked -- the module-level
    MO2matrix reaches the ED object's single-site building blocks, not
    its sector-aware entry points."""
    ed = chain.get_ED_obj()
    H = np.array(MO2matrix(h, ed).todense())
    Q = np.array(MO2matrix(charge, ed).todense())
    sel = np.abs(np.diag(Q).real - target) < 1e-9
    return float(np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0])


############################################################
# 1) addition spectrum of a spinless t-V chain, sector by sector
############################################################

n, V = 8, 2.0
fc = fermionchain.Fermionic_Chain(n)
fc.mode = "ED"                      # this chain is solved by ED
h = 0
for i in range(n - 1):
    h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
    h = h + V * fc.N[i] * fc.N[i + 1]
fc.set_hamiltonian(h)
number = sum(fc.N)

e_global = fc.gs_energy(mode="ED")   # the unconstrained answer
nfs = np.arange(n + 1)
efs, refs = [], []
for nf in nfs:
    fc.set_conserved_sector(Nf=int(nf))
    efs.append(fc.gs_energy(mode="ED"))
    refs.append(sector_energy_by_hand(fc, h, number, nf))
    # the state really is in that sector, not merely close in energy
    assert abs(fc.vev(number, mode="ED").real - nf) < 1e-8
efs, refs = np.array(efs), np.array(refs)
print("t-V chain, max |E_sector - reference| =", np.max(np.abs(efs - refs)))
assert np.max(np.abs(efs - refs)) < 1e-8
# the lowest sector must reproduce the unconstrained ground state
assert abs(np.min(efs) - e_global) < 1e-8

# charge gap Delta(N) = E(N+1)+E(N-1)-2E(N)
gaps = efs[2:] + efs[:-2] - 2 * efs[1:-1]

############################################################
# 2) total-Sz sectors of a Heisenberg chain
############################################################

ns = 8
sc = spinchain.Spin_Chain([2] * ns)
sc.mode = "ED"
hs = 0
for i in range(ns - 1):
    hs = hs + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
            + sc.Sz[i] * sc.Sz[i + 1]
sc.set_hamiltonian(hs)
sz = sum(sc.Sz)

two_szs = np.arange(-ns, ns + 1, 2)   # ITensor's integer 2*Sz units
ess, refss = [], []
for two_sz in two_szs:
    sc.set_conserved_sector(Sz=int(two_sz))
    ess.append(sc.gs_energy(mode="ED"))
    refss.append(sector_energy_by_hand(sc, hs, sz, two_sz / 2.0))
    assert abs(sc.vev(sz, mode="ED").real - two_sz / 2.0) < 1e-8
ess, refss = np.array(ess), np.array(refss)
print("Heisenberg chain, max |E_sector - reference| =",
      np.max(np.abs(ess - refss)))
assert np.max(np.abs(ess - refss)) < 1e-8

# a charge-changing operator is refused, not silently zeroed: restricting
# it to the sector would give an exact, meaningless zero
sc.set_conserved_sector(Sz=0)
try:
    sc.vev(sc.Sx[0], mode="ED")
    raise AssertionError("Sx should not be measurable inside an Sz sector")
except ValueError as e:
    print("Sx in an Sz=0 sector ->", str(e)[:70], "...")

############################################################
# plots
############################################################

fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))

axs[0].plot(nfs, refs, "o-", label="restricted-ED reference")
axs[0].plot(nfs, efs, "x", ms=10, mew=2, label='set_conserved_sector, mode="ED"')
axs[0].axhline(e_global, color="gray", ls="--", lw=1,
               label="unconstrained ground state")
axs[0].set_xlabel("particle number $N_f$")
axs[0].set_ylabel("$E_0(N_f)$")
axs[0].set_title("$t$-$V$ chain, $n=%d$, $V=%g$" % (n, V))
axs[0].legend(fontsize=8)

axs[1].plot(nfs[1:-1], gaps, "s-", color="C3")
axs[1].axhline(0.0, color="gray", lw=0.8)
axs[1].set_xlabel("particle number $N_f$")
axs[1].set_ylabel(r"$E_0(N{+}1)+E_0(N{-}1)-2E_0(N)$")
axs[1].set_title("charge gap")

axs[2].plot(two_szs / 2.0, refss, "o-", label="restricted-ED reference")
axs[2].plot(two_szs / 2.0, ess, "x", ms=10, mew=2, label='mode="ED" sector')
axs[2].set_xlabel(r"total $S^z$")
axs[2].set_ylabel(r"$E_0(S^z)$")
axs[2].set_title("Heisenberg chain, $n=%d$" % ns)
axs[2].legend(fontsize=8)

fig.tight_layout()
plt.savefig("ed_sector_addition_spectrum.png", dpi=130)
plt.show()
