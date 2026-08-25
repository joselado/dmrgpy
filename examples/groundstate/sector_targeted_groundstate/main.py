# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground states of one *chosen* quantum-number sector
# (Many_Body_Chain.set_conserved_sector, itensor_version=3 only).
#
# By default DMRG here searches the whole Hilbert space and returns the
# global ground state, whatever particle number or total Sz that happens
# to have. set_conserved_sector(Nf=k) instead rebuilds the chain's sites
# with QN-carrying indices, so the entire calculation -- the starting
# state, every sweep, the returned state -- lives at exactly k particles.
#
# That gives the grand-canonical-free way to get E(N): the addition
# spectrum, the charge gap, the compressibility, all as literal
# ground-state energies of separate sectors rather than as artefacts of a
# chemical potential. Two consequences worth knowing: every operator built
# on the chain must itself conserve the requested quantity (a bare C, or
# an Sx on an Sz-conserving chain, raises a ValueError naming it), and the
# cost changes -- the quantum numbers restore block sparsity, which is
# worth more than its bookkeeping only once the solve is big enough (see
# the third panel: slower below ~24 sites here, several times faster
# above).
#
# The model is a spinless t-V chain: nearest-neighbour hopping plus a
# nearest-neighbour repulsion V, which at half filling drives a
# charge-density-wave. The script sweeps the particle number over every
# sector, checks each against ED, and plots
#   (1) E(N) -- the addition spectrum -- with the ED reference,
#   (2) the discrete second derivative E(N+1)+E(N-1)-2E(N), i.e. the
#       charge gap, which peaks at half filling for large V,
#   (3) sector vs. dense wall time against chain length, i.e. where the
#       block sparsity starts paying for its own bookkeeping.
import time

import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import fermionchain
from dmrgpy.multioperator import MO2matrix

V = 2.0  # nearest-neighbour repulsion


def tV_chain(n):
    """A spinless t-V chain, its Hamiltonian and its number operator."""
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + V * fc.N[i] * fc.N[i + 1]
    return fc, h, sum(fc.N)


def ed_sector_energy(chain, h, number, target):
    """Exact reference: the number operator is diagonal in the ED product
    basis, so one sector is a set of basis states and its ground-state
    energy is the lowest eigenvalue of the corresponding submatrix. (Do
    NOT instead diagonalize the full H and filter eigenvectors by <N>:
    eigenvalues degenerate across sectors come back as arbitrary
    superpositions belonging to no sector at all.)"""
    ed = chain.get_ED_obj()
    H = np.array(MO2matrix(h, ed).todense())
    N = np.array(MO2matrix(number, ed).todense())
    sel = np.abs(np.diag(N).real - target) < 1e-9
    return float(np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0])


# --- (1)+(2) the addition spectrum of a 10-site chain ------------------
n = 10
fc, h, number = tV_chain(n)
fc.set_hamiltonian(h)

fillings = list(range(n + 1))
e_dmrg, e_ed, occupations = [], [], []
for nf in fillings:
    fc.set_conserved_sector(Nf=nf)
    e_dmrg.append(fc.gs_energy())
    e_ed.append(ed_sector_energy(fc, h, number, nf))
    occupations.append([fc.vev(fc.N[i]).real for i in range(n)])
    print("Nf = %2d   E = %.10f   (ED %.10f)   <N> = %.8f"
          % (nf, e_dmrg[-1], e_ed[-1], fc.vev(number).real))

e_dmrg = np.array(e_dmrg)
e_ed = np.array(e_ed)
assert np.max(np.abs(e_dmrg - e_ed)) < 1e-6, "sector DMRG disagrees with ED"

# the sector-free ground state must be the best of the sectors
fc.set_conserved_sector()
e_global = fc.gs_energy()
print("global ground state E = %.10f, best sector E = %.10f (Nf = %d)"
      % (e_global, e_dmrg.min(), int(np.argmin(e_dmrg))))
assert abs(e_global - e_dmrg.min()) < 1e-6

charge_gap = e_dmrg[2:] + e_dmrg[:-2] - 2 * e_dmrg[1:-1]

# --- (3) what the block sparsity is worth ------------------------------
# Conserving a quantum number makes every tensor block-sparse: bookkeeping
# per block, less arithmetic inside each. Which side wins depends on how
# big the solve is, so sweep the chain length at a fixed bond dimension.
# The model here is the particle-hole symmetric point (an extra -V*N), so
# the global ground state *is* half filling -- which makes the dense and
# sector runs solve for the same state and lets the comparison double as
# a correctness check. (Pin BLAS to one thread when reading these numbers:
# MKL_NUM_THREADS=1 OMP_NUM_THREADS=1, see src/dmrgpy/blasthreads.py.)
def tV_symmetric(n):
    fc, h, number = tV_chain(n)
    return fc, h - sum(V * fc.N[i] for i in range(n)), number


timing_maxm = 100
lengths = [16, 24, 32, 40]
t_dense, t_sector = [], []
for length in lengths:
    fcl, hl, _ = tV_symmetric(length)
    fcl.maxm = timing_maxm
    fcl.set_hamiltonian(hl)
    t0 = time.time() ; e_d = fcl.gs_energy() ; t_dense.append(time.time() - t0)

    fcl2, hl2, _ = tV_symmetric(length)
    fcl2.maxm = timing_maxm
    fcl2.set_hamiltonian(hl2)
    fcl2.set_conserved_sector(Nf=length // 2)
    t0 = time.time() ; e_s = fcl2.gs_energy() ; t_sector.append(time.time() - t0)

    # same state, so the two runs must agree -- the timing comparison is
    # only meaningful because of this
    assert abs(e_d - e_s) < 1e-6, (length, e_d, e_s)
    print("n = %2d   dense %6.2fs   sector %6.2fs   %.2fx   (E = %.8f)"
          % (length, t_dense[-1], t_sector[-1], t_dense[-1] / t_sector[-1], e_d))

# --- plots -------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

ax = axes[0]
ax.plot(fillings, e_ed, "o-", label="ED (sector-restricted)")
ax.plot(fillings, e_dmrg, "x--", label="DMRG, set_conserved_sector(Nf=N)")
ax.axhline(e_global, color="grey", ls=":", label="unconstrained ground state")
ax.set_xlabel("particle number N")
ax.set_ylabel("ground-state energy of the sector")
ax.set_title("Addition spectrum, t-V chain (n=%d, V=%g)" % (n, V))
ax.legend()

ax = axes[1]
ax.plot(fillings[1:-1], charge_gap, "s-")
ax.axvline(n / 2, color="grey", ls=":", label="half filling")
ax.set_xlabel("particle number N")
ax.set_ylabel(r"$E(N{+}1)+E(N{-}1)-2E(N)$")
ax.set_title("Charge gap")
ax.legend()

ax = axes[2]
ax.plot(lengths, t_dense, "o-", label="dense (default)")
ax.plot(lengths, t_sector, "s-", label="sector, Nf = n/2")
ax.set_yscale("log")
ax.set_xlabel("chain length")
ax.set_ylabel("ground-state wall time [s]")
ax.set_title("Same state, both ways (maxm=%d)" % timing_maxm)
ax.legend()

fig.tight_layout()
plt.show()
