# Applying an operator to every unit cell of a converged infinite MPS --
# Infinite_Many_Body_Chain.apply_operator, the infinite-chain analogue of
# the finite backends' applyMPO: the operator is a MultiOperator on the
# central cell (ic.SxC[0], ic.SzC[0]*ic.SzC[1], a sum of such terms), the
# same one in every cell, and the result is a NEW chain holding
# (prod over cells W)|psi>, whose vev/correlator measure that state and
# whose .eta is <psi|W^dag W|psi> per unit cell (1 for a unitary
# operator). Underneath, pyitensor.idmrg.apply_mpo contracts the periodic
# MPO onto the converged unit cell and re-canonicalizes/truncates the
# grown bond dimension via the standard two-sided fixed-point procedure
# (Orus & Vidal, PRB 78, 155117 (2008)).
#
# SCOPE: a per-cell operator is a bounded operator, which is what the
# construction needs; the Hamiltonian's own extensive sum over cells is
# not of this kind (see pyitensor/idmrg.py's own "Applying a (bounded)
# MPO to the converged iMPS" section docstring for why), and every factor
# must sit on the central cell (a factor on SxR/SxL raises ValueError).
#
# Demonstrated below: (1) a chi_W=1 unitary single-site operator (Pauli-X
# at every site of a uniform Heisenberg chain) -- flips <Sz> and leaves
# <Sz(0)Sz(r)> unchanged, both exact, closed-form checks; (2) a genuinely
# bond-dimension>1 local gate, exp(-i*theta*S0.S1) on the two sites of a
# two-site cell, written as alpha*Id + beta*S0.S1 (S.S has eigenvalues 1/4
# and -3/4) -- exercises the actual bond-growth path, and being unitary
# must preserve the norm (eta=1) and commute with S0.S1, so the intra-cell
# correlator is unchanged while the inter-cell one rotates with theta.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import infinitechain

################################################################
### (1) chi_W=1 unitary single-site operator: Pauli-X at every site ###
################################################################
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "idmrg"
h = ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]
ic.set_hamiltonian(h)
ic.maxm, ic.maxiter, ic.etol = 30, 60, 1e-9
ic.gs_energy()

flipped = ic.apply_operator(2 * ic.SxC[0])  # Pauli-X is unitary for spin-1/2

print("=== Pauli-X applied to every site of a uniform Heisenberg chain ===")
print("norm diagnostic eta (unitary -> should be ~1):", flipped.eta)
sz0 = ic.vev("Sz", 0)
sz0_flipped = flipped.vev("Sz", 0)
print("<Sz> before:", sz0, " after:", sz0_flipped, " (expect after = -before)")
rs = range(1, 6)
c0_list, c1_list = [], []
for r in rs:
    c0 = ic.correlator("Sz", 0, "Sz", r)
    c1 = flipped.correlator("Sz", 0, "Sz", r)
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
ic2.gs_method = "idmrg"
h2 = (ic2.SxC[0] * ic2.SxC[1] + ic2.SyC[0] * ic2.SyC[1] + ic2.SzC[0] * ic2.SzC[1]
      + 0.4 * (ic2.SxC[1] * ic2.SxR[0] + ic2.SyC[1] * ic2.SyR[0] + ic2.SzC[1] * ic2.SzR[0]))
ic2.set_hamiltonian(h2)
# Deliberately modest maxm: apply_operator's own truncation/regauging step
# solves a numerically delicate two-sided-fixed-point problem whose
# conditioning degrades the further the *raw* grown bond dimension
# (chi_A*chi_W, before any truncation) exceeds the state's real
# entanglement at that cut -- see idmrg.py's _canonicalize_periodic's own
# comment. The gate below has chi_W=4, so the grown bond is 16; maxdim is
# set to keep it, since truncating back to maxm=4 would discard about
# 1e-3 of weight, more than the gauge fixing can absorb.
ic2.maxm, ic2.maxiter, ic2.etol = 4, 100, 1e-12
ic2.gs_energy()

ss = ic2.SxC[0] * ic2.SxC[1] + ic2.SyC[0] * ic2.SyC[1] + ic2.SzC[0] * ic2.SzC[1]


def heisenberg_gate(theta):
    """exp(-i*theta*S0.S1) as a MultiOperator: S.S has eigenvalues 1/4
    (triplet) and -3/4 (singlet), so the exponential is alpha*Id +
    beta*S.S with alpha + beta/4 = exp(-i*theta/4) and
    alpha - 3*beta/4 = exp(3i*theta/4)."""
    beta = np.exp(-1j * theta / 4) - np.exp(3j * theta / 4)
    alpha = np.exp(-1j * theta / 4) - beta / 4
    return alpha * ic2.get_operator("Id", 0) + beta * ss


thetas = [0.0, 0.2, 0.37, 0.6, 0.9]
intra0 = sum(ic2.correlator(n, 0, n, 1) for n in ("Sx", "Sy", "Sz")).real
inter0 = ic2.correlator("Sz", 1, "Sz", 1).real  # site 1 and site 0 of the next cell
intra_list, inter_list, eta_list = [], [], []
print()
print("=== 2-site unitary gate exp(-i theta S0.S1) applied once per unit cell ===")
for theta in thetas:
    gated = ic2.apply_operator(heisenberg_gate(theta), cutoff=1e-10, maxdim=16)
    intra = sum(gated.correlator(n, 0, n, 1) for n in ("Sx", "Sy", "Sz")).real
    inter = gated.correlator("Sz", 1, "Sz", 1).real
    intra_list.append(intra)
    inter_list.append(inter)
    eta_list.append(gated.eta.real if isinstance(gated.eta, complex) else gated.eta)
    print("theta={:.2f}  eta={:.10f}  <S0.S1> intra={:.8f} (expect {:.8f})  "
          "<Sz Sz> inter={:.8f}".format(theta, eta_list[-1], intra, intra0, inter))
    assert abs(eta_list[-1] - 1.0) < 1e-4
    assert abs(intra - intra0) < 1e-6  # the gate commutes with S0.S1
assert abs(inter_list[-1] - inter0) > 1e-3  # the inter-cell bond does move

print()
print("apply_operator example PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

ax1.plot(list(rs), c0_list, "o-", color="black", label="before Pauli-X")
ax1.plot(list(rs), c1_list, "x--", color="tab:red", label="after Pauli-X")
ax1.set_xlabel("distance $r$")
ax1.set_ylabel(r"$\langle S_z(0)S_z(r)\rangle$")
ax1.set_title("(1) chi_W=1: Pauli-X on every site")
ax1.legend()

ax2.plot(thetas, intra_list, "o-", color="tab:blue", label=r"intra-cell $\langle S_0\cdot S_1\rangle$ (unchanged)")
ax2.plot(thetas, inter_list, "s--", color="tab:orange", label=r"inter-cell $\langle S_z(1) S_z(2)\rangle$")
ax2.set_xlabel(r"gate angle $\theta$")
ax2.set_ylabel("correlator")
ax2.set_title("(2) chi_W>1: 2-site unitary gate\n(eta = 1 at every angle)")
ax2.legend()

fig.suptitle("apply_operator: a per-cell operator applied to a converged iMPS")
fig.tight_layout()
fig.savefig("apply_mpo_to_infinite_chain.png", dpi=150)
print("saved plot to apply_mpo_to_infinite_chain.png")
