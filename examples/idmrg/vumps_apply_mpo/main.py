# Applying a periodic (bounded) MPO to a converged VUMPS iMPS --
# pyitensor.vumps.apply_mpo, the VUMPS-mixed-gauge analogue of
# pyitensor.idmrg.apply_mpo (see examples/idmrg/apply_mpo_to_infinite_chain
# for the idmrg-side counterpart this mirrors, and
# examples/idmrg/vumps_imps_sum for the sibling VUMPS-mixed-gauge
# operation, imps_sum, which shares most of the same construction).
#
# Construction: group W_bulk into a single grouped-supersite MPO tensor
# (vumps._group_automaton, the same routine that groups VUMPS's own
# Hamiltonian automaton), grow the converged AL by it
# (idmrg.grow_by_mpo, on the trivial n_uc=1 "periodic chain" VUMPS already
# works at), re-canonicalize/truncate (idmrg._canonicalize_periodic, the
# standard two-sided fixed-point infinite-MPS procedure, Orus & Vidal, PRB
# 78, 155117 (2008)), then complete the truncated left-canonical tensor
# back to the full mixed gauge {AL,AR,C,AC} (vumps._complete_mixed_gauge,
# Vanderstraeten/Haegeman/Verstraete, arXiv:1810.07006, Eq.(9)-(17)).
#
# SCOPE: identical to idmrg.apply_mpo's own -- W_bulk must be a genuinely
# *bounded* (non-extensive) periodic operator, the same tensor reused at
# every unit cell with no unconditional "keep accumulating forever"
# self-loop (the Hamiltonian's own automaton is explicitly out of scope).
# `W_bulk` uses EXACTLY idmrg.apply_mpo's own convention (a list of n_uc
# ITensors, one per unit-cell sublattice site), so the identical W_bulk
# list built once below is fed to both idmrg.apply_mpo and
# vumps.apply_mpo to cross-check the two backends' own results against
# each other on the same physical operator.
#
# Demonstrated below: (1) a chi_W=1 unitary single-site operator
# (Pauli-X), applied to BOTH a converged idmrg iMPS and a converged vumps
# iMPS at the exact, D=1, field-polarized point -- flips <Sz> on both,
# cross-checked directly against each other (no convergence-tolerance
# ambiguity to hide a real discrepancy behind); (2) the same operator
# applied to a genuinely entangled D>1 VUMPS ground state (TFIM, gapped)
# -- flips <Sz>, leaves <Sz(0)Sz(r)> unchanged, preserves the norm
# diagnostic `eta`; (3) a genuinely bond-dimension>1 local gate (an
# SVD-split 2-site unitary rotation, tiled once per an n_uc=2 unit cell)
# -- exercises the actual bond-growth path, and being unitary must also
# preserve `eta`.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor import vumps
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor


def pauli_x_mpo(sites_uc):
    d = sites_uc.dim(1)
    pauli_x = 2 * sites_uc.site_type(1).matrix("Sx")  # unitary for spin-1/2
    link_l, link_r = Index(1, tags="Link"), Index(1, tags="Link")
    s = sites_uc.si(1)
    return [ITensor((link_l, s, s.prime(1), link_r), pauli_x.reshape(1, d, d, 1))]


def field_chain(B, gs_method):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = gs_method
    ic.maxm = 1
    ic.set_hamiltonian(-B * ic.SzC[0])
    ic.gs_energy()
    return ic


def tfim_chain(g, D):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = D
    ic.etol = 1e-11
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


############################################################################
### (1) chi_W=1 Pauli-X flip at the exact D=1 point -- vumps VS idmrg ###
############################################################################
ic_v1 = field_chain(5.0, "vumps")
ic_i1 = field_chain(5.0, "idmrg")
W_pauli_x = pauli_x_mpo(ic_v1._vumps_result.sites_uc)

flipped_v1 = vumps.apply_mpo(ic_v1._vumps_result, W_pauli_x, cutoff=1e-12, maxdim=None)
flipped_i1 = idmrg.apply_mpo(ic_i1._result, W_pauli_x, cutoff=1e-12, maxdim=None)

sz_v_before = vumps.onsite_expectation(ic_v1._vumps_result, "Sz", 0)
sz_v_after = vumps.onsite_expectation(flipped_v1, "Sz", 0)
sz_i_before = idmrg.onsite_expectation(ic_i1._result, "Sz", 0)
sz_i_after = idmrg.onsite_expectation(flipped_i1, "Sz", 0)

print("=== (1) Pauli-X on the exact D=1 field-polarized point ===")
print("<Sz> before -- vumps: {}  idmrg: {}".format(sz_v_before, sz_i_before))
print("<Sz> after  -- vumps: {}  idmrg: {}  (expect after=-before, vumps==idmrg)".format(
    sz_v_after, sz_i_after))
assert abs(sz_v_after - (-sz_v_before)) < 1e-8
assert abs(sz_i_after - (-sz_i_before)) < 1e-8
assert abs(sz_v_after - sz_i_after) < 1e-8
assert abs(flipped_v1.eta - 1.0) < 1e-8
print("PASSED (vumps and idmrg agree exactly)")

##############################################################################
### (2) chi_W=1 Pauli-X flip on a genuinely entangled D>1 VUMPS state (TFIM) ###
##############################################################################
D = 2
ic_tfim = tfim_chain(1.5, D)
flipped_tfim = vumps.apply_mpo(ic_tfim._vumps_result, W_pauli_x, cutoff=1e-12, maxdim=None)

sz0 = vumps.onsite_expectation(ic_tfim._vumps_result, "Sz", 0)
sz0_flipped = vumps.onsite_expectation(flipped_tfim, "Sz", 0)
rs = list(range(1, 6))
c0_list = [vumps.two_point_correlator(ic_tfim._vumps_result, "Sz", 0, "Sz", r).real for r in rs]
c1_list = [vumps.two_point_correlator(flipped_tfim, "Sz", 0, "Sz", r).real for r in rs]

print()
print("=== (2) Pauli-X on a genuinely entangled D={} TFIM VUMPS ground state ===".format(D))
print("<Sz> before:", sz0, " after:", sz0_flipped, " (expect after=-before)")
print("norm diagnostic eta (unitary -> should be ~1):", flipped_tfim.eta)
assert abs(sz0_flipped - (-sz0)) < 1e-6
assert abs(flipped_tfim.eta - 1.0) < 1e-6
for c0, c1 in zip(c0_list, c1_list):
    assert abs(c1 - c0) < 1e-6
print("PASSED")

########################################################################
### (3) chi_W>1 local gate: a 2-site unitary rotation, one per unit cell ###
########################################################################
ic2 = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
ic2.gs_method = "vumps"
h2 = (ic2.SxC[0] * ic2.SxC[1] + ic2.SyC[0] * ic2.SyC[1] + ic2.SzC[0] * ic2.SzC[1]
      + 0.4 * (ic2.SxC[1] * ic2.SxR[0] + ic2.SyC[1] * ic2.SyR[0] + ic2.SzC[1] * ic2.SzR[0]))
ic2.set_hamiltonian(h2)
# Deliberately modest maxm -- apply_mpo's own truncation/regauging step
# solves a numerically delicate two-sided fixed-point problem whose
# conditioning degrades the further the *raw* grown bond dimension
# (chi_A*chi_W, before any truncation) exceeds the state's real
# entanglement at that cut, same as idmrg.apply_mpo's own example.
ic2.maxm, ic2.maxiter, ic2.etol = 3, 120, 1e-10
ic2.gs_energy()

sites_uc2 = ic2._vumps_result.sites_uc
Sx, Sy, Sz = (sites_uc2.site_type(1).matrix(n) for n in ("Sx", "Sy", "Sz"))
d2 = sites_uc2.dim(1)
H2 = (np.kron(Sx, Sx) + np.kron(Sy, Sy) + np.kron(Sz, Sz)).real
gate = expm(-1j * 0.37 * H2)  # d^2 x d^2 unitary 2-site rotation
gate4 = np.transpose(gate.reshape(d2, d2, d2, d2), (2, 0, 3, 1))  # (s0,s0',s1,s1')
U, S, Vh = np.linalg.svd(gate4.reshape(d2 * d2, d2 * d2), full_matrices=False)
keep = int(np.sum(S > 1e-12))
U, S, Vh = U[:, :keep], S[:keep], Vh[:keep, :]
a_half = (U * S[None, :] ** 0.5).reshape(d2, d2, keep)
b_half = (S[:, None] ** 0.5 * Vh).reshape(keep, d2, d2)

s0, s1 = sites_uc2.si(1), sites_uc2.si(2)
left_dummy, mid, right_dummy = Index(1, tags="Link"), Index(keep, tags="Link"), Index(1, tags="Link")
W0 = ITensor((left_dummy, s0, s0.prime(1), mid), a_half.reshape(1, d2, d2, keep))
W1 = ITensor((mid, s1, s1.prime(1), right_dummy), b_half.reshape(keep, d2, d2, 1))

gated = vumps.apply_mpo(ic2._vumps_result, [W0, W1], cutoff=1e-10, maxdim=None)

print()
print("=== (3) 2-site unitary gate applied once per unit cell (VUMPS, n_uc=2) ===")
print("original bond dim D:", ic2._vumps_result.D, " gated bond dim D:", gated.D)
print("norm diagnostic eta (unitary -> should be ~1):", gated.eta)
assert abs(gated.eta - 1.0) < 1e-4
print("PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

ax1.plot(rs, c0_list, "o-", color="black", label="before Pauli-X")
ax1.plot(rs, c1_list, "x--", color="tab:red", label="after Pauli-X")
ax1.set_xlabel("distance $r$")
ax1.set_ylabel(r"$\langle S_z(0)S_z(r)\rangle$")
ax1.set_title("(2) chi_W=1: Pauli-X on a D={} TFIM state".format(D))
ax1.legend()

labels = ["original", "gated"]
ax2.bar(labels, [ic2._vumps_result.D, gated.D], color=["tab:blue", "tab:orange"])
ax2.set_ylabel("VUMPS bond dimension D")
ax2.set_title("(3) chi_W>1: 2-site unitary gate\n(eta={:.6f})".format(gated.eta.real))

fig.suptitle("vumps.apply_mpo: bounded-MPO application to a converged VUMPS iMPS")
fig.tight_layout()
fig.savefig("vumps_apply_mpo.png", dpi=150)
print("saved plot to vumps_apply_mpo.png")
