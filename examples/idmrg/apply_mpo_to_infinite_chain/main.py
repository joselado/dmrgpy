# Applying a periodic (bounded) MPO to the converged infinite MPS --
# pyitensor.idmrg.apply_mpo, the infinite-chain analogue of the finite
# backends' applyMPO: contract a periodic MPO onto the converged unit
# cell (grow_by_mpo), then re-canonicalize/truncate the grown bond
# dimension back down (idmrg._canonicalize_periodic) via the standard
# two-sided fixed-point infinite-MPS canonicalization procedure.
#
# SCOPE: apply_mpo only works for *bounded* (non-extensive) periodic
# operators -- the same tensor reused at every unit cell, with no
# unconditional "keep accumulating forever" self-loop -- see
# pyitensor/idmrg.py's own "Applying a (bounded) MPO to the converged
# iMPS" section docstring for why idmrg.py's own Hamiltonian automaton
# (built by _build_periodic_mpo) is explicitly *not* a valid input here.
#
# Demonstrated below: (1) a chi_W=1 unitary single-site operator (Pauli-X
# at every site of a uniform Heisenberg chain) -- flips <Sz> and leaves
# <Sz(0)Sz(r)> unchanged, both exact, closed-form checks; (2) a genuinely
# bond-dimension>1 local gate (an SVD-split 2-site unitary rotation,
# tiled once per unit cell) -- exercises the actual bond-growth path,
# and being unitary must preserve the norm (apply_mpo's own `eta`
# diagnostic).

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor

################################################################
### (1) chi_W=1 unitary single-site operator: Pauli-X at every site ###
################################################################
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
h = ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]
ic.set_hamiltonian(h)
ic.maxm, ic.maxiter, ic.etol = 30, 60, 1e-9
ic.gs_energy()

sites_uc = ic._result.sites_uc
d = sites_uc.dim(1)
pauli_x = 2 * sites_uc.site_type(1).matrix("Sx")  # unitary for spin-1/2
link_l, link_r = Index(1, tags="Link"), Index(1, tags="Link")
s = sites_uc.si(1)
W_pauli_x = [ITensor((link_l, s, s.prime(1), link_r), pauli_x.reshape(1, d, d, 1))]

flipped = idmrg.apply_mpo(ic._result, W_pauli_x, cutoff=1e-12, maxdim=None)

print("=== Pauli-X applied to every site of a uniform Heisenberg chain ===")
print("norm diagnostic eta (unitary -> should be ~1):", flipped.eta)
sz0 = idmrg.onsite_expectation(ic._result, "Sz", 0)
sz0_flipped = idmrg.onsite_expectation(flipped, "Sz", 0)
print("<Sz> before:", sz0, " after:", sz0_flipped, " (expect after = -before)")
rs = range(1, 6)
c0_list, c1_list = [], []
for r in rs:
    c0 = idmrg.two_point_correlator(ic._result, "Sz", 0, "Sz", r)
    c1 = idmrg.two_point_correlator(flipped, "Sz", 0, "Sz", r)
    c0_list.append(c0)
    c1_list.append(c1)
    print("<Sz(0)Sz({})> before: {}  after: {}  (expect unchanged)".format(r, c0, c1))
    assert abs(c1 - c0) < 1e-6

assert abs(sz0_flipped - (-sz0)) < 1e-6
assert abs(flipped.eta - 1.0) < 1e-6

########################################################################
### (2) chi_W>1 local gate: a 2-site unitary rotation, one per unit cell ###
########################################################################
ic2 = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
h2 = (ic2.SxC[0] * ic2.SxC[1] + ic2.SyC[0] * ic2.SyC[1] + ic2.SzC[0] * ic2.SzC[1]
      + 0.4 * (ic2.SxC[1] * ic2.SxR[0] + ic2.SyC[1] * ic2.SyR[0] + ic2.SzC[1] * ic2.SzR[0]))
ic2.set_hamiltonian(h2)
# Deliberately modest maxm: apply_mpo's own truncation/regauging step
# solves a numerically delicate two-sided-fixed-point problem whose
# conditioning degrades the further the *raw* grown bond dimension
# (chi_A*chi_W, before any truncation) exceeds the state's real
# entanglement at that cut -- see idmrg.py's _canonicalize_periodic's own
# comment. A larger maxm here still works, just needs a looser internal
# tolerance to accommodate the resulting ill-conditioning.
ic2.maxm, ic2.maxiter, ic2.etol = 4, 100, 1e-12
ic2.gs_energy()

sites_uc2 = ic2._result.sites_uc
Sx, Sy, Sz = (sites_uc2.site_type(1).matrix(n) for n in ("Sx", "Sy", "Sz"))
H2 = (np.kron(Sx, Sx) + np.kron(Sy, Sy) + np.kron(Sz, Sz)).real
gate = expm(-1j * 0.37 * H2)  # d^2 x d^2 unitary 2-site rotation
gate4 = np.transpose(gate.reshape(d, d, d, d), (2, 0, 3, 1))  # (s0,s0',s1,s1')
U, S, Vh = np.linalg.svd(gate4.reshape(d * d, d * d), full_matrices=False)
keep = int(np.sum(S > 1e-12))
U, S, Vh = U[:, :keep], S[:keep], Vh[:keep, :]
a_half = (U * S[None, :] ** 0.5).reshape(d, d, keep)
b_half = (S[:, None] ** 0.5 * Vh).reshape(keep, d, d)

s0, s1 = sites_uc2.si(1), sites_uc2.si(2)
left_dummy, mid, right_dummy = Index(1, tags="Link"), Index(keep, tags="Link"), Index(1, tags="Link")
W0 = ITensor((left_dummy, s0, s0.prime(1), mid), a_half.reshape(1, d, d, keep))
W1 = ITensor((mid, s1, s1.prime(1), right_dummy), b_half.reshape(keep, d, d, 1))

gated = idmrg.apply_mpo(ic2._result, [W0, W1], cutoff=1e-10, maxdim=None)

print()
print("=== 2-site unitary gate applied once per unit cell ===")
print("original bond dims:", [u.array.shape for u in ic2._result.U_list])
print("gated bond dims:   ", [u.array.shape for u in gated.U_list])
print("norm diagnostic eta (unitary -> should be ~1):", gated.eta)
assert abs(gated.eta - 1.0) < 1e-4

print()
print("apply_mpo example PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

ax1.plot(list(rs), c0_list, "o-", color="black", label="before Pauli-X")
ax1.plot(list(rs), c1_list, "x--", color="tab:red", label="after Pauli-X")
ax1.set_xlabel("distance $r$")
ax1.set_ylabel(r"$\langle S_z(0)S_z(r)\rangle$")
ax1.set_title("(1) chi_W=1: Pauli-X on every site")
ax1.legend()

orig_dims = [u.array.shape[-1] for u in ic2._result.U_list]
gated_dims = [u.array.shape[-1] for u in gated.U_list]
bond_idx = np.arange(len(orig_dims))
width = 0.35
ax2.bar(bond_idx - width/2, orig_dims, width, label="original")
ax2.bar(bond_idx + width/2, gated_dims, width, label="gated")
ax2.set_xlabel("bond index (per unit-cell site)")
ax2.set_ylabel("bond dimension")
ax2.set_xticks(bond_idx)
ax2.set_title("(2) chi_W>1: 2-site unitary gate")
ax2.legend()

fig.suptitle("apply_mpo: bounded-MPO application to a converged iMPS")
fig.tight_layout()
fig.savefig("apply_mpo_to_infinite_chain.png", dpi=150)
print("saved plot to apply_mpo_to_infinite_chain.png")
