# Cross-checks the ITensor v3 C++ port of VUMPS's apply_mpo
# (mpscpp3/chain_session.h's Chain::vumps_apply_mpo) against
# pyitensor.vumps.apply_mpo -- see examples/idmrg/vumps_apply_mpo for the
# itensor_version="python"-only version of this example, which this one
# extends with a direct v3 comparison at every step. There is no
# Infinite_Many_Body_Chain-level wrapper for apply_mpo on EITHER backend
# (docs/user_guide.md's own "Applying an operator/gate to a converged
# VUMPS iMPS" section), so v3's own apply_mpo is reached directly via
# ic._session3.vumps_apply_mpo(...), working against ic._session3's own
# converged VUMPS snapshot; the result is written back via
# vumps_load_uniform_state so ic.vev()/ic.correlator() see it (dmrgpy's
# own Python-level dispatch is otherwise untouched by this).
#
# Demonstrated: (1) a chi_W=1 unitary Pauli-X flip at the exact D=1
# field-polarized point -- v3 and pyitensor cross-checked directly against
# each other (no convergence-tolerance ambiguity to hide a real
# discrepancy behind); (2) the same flip swept over bond dimension D on a
# genuinely entangled TFIM ground state -- <Sz> and a multi-r correlator
# on both backends, plus the norm diagnostic eta (should stay ~1, unitary
# operator); (3) a genuinely bond-dimension>1 local gate (an SVD-split
# 2-site unitary rotation, tiled once per an n_uc=2 unit cell) on both
# backends -- exercises the actual bond-growth path.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor

_PAULI_X = np.array([[0, 1], [1, 0]], dtype=complex)
_SX = 0.5 * _PAULI_X
_SY = 0.5 * np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = 0.5 * np.array([[1, 0], [0, -1]], dtype=complex)


def field_chain(B, itensor_version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(-B * ic.SzC[0])
    ic.gs_energy()
    return ic


def tfim_chain(g, D, itensor_version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "vumps"
    ic.maxm = D
    ic.etol = 1e-11
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


def pauli_x_mpo_python(sites_uc):
    d = sites_uc.dim(1)
    link_l, link_r = Index(1, tags="Link"), Index(1, tags="Link")
    s = sites_uc.si(1)
    return [ITensor((link_l, s, s.prime(1), link_r), _PAULI_X.reshape(1, d, d, 1))]


def apply_pauli_x_v3(ic):
    """Applies the chi_W=1 Pauli-X to ic._session3's own converged VUMPS
    snapshot and loads the result back -- see this file's own module
    docstring."""
    W = [_PAULI_X.flatten().tolist()]
    D, d_g, AL, AR, C, AC, eta = ic._session3.vumps_apply_mpo(W, [1], [1], 1e-12, 0)
    ic._session3.vumps_load_uniform_state(
        D, d_g, AL.flatten().tolist(), AR.flatten().tolist(), C.flatten().tolist())
    return eta


############################################################################
### (1) chi_W=1 Pauli-X flip at the exact D=1 point -- v3 VS pyitensor ###
############################################################################
ic_py1 = field_chain(5.0, "python")
ic_v31 = field_chain(5.0, 3)
W_py = pauli_x_mpo_python(ic_py1._vumps_result.sites_uc)

sz_py_before = vumps.onsite_expectation(ic_py1._vumps_result, "Sz", 0)
sz_v3_before = ic_v31.vev("Sz", 0)

flipped_py1 = vumps.apply_mpo(ic_py1._vumps_result, W_py, cutoff=1e-12, maxdim=None)
eta_v31 = apply_pauli_x_v3(ic_v31)

sz_py_after = vumps.onsite_expectation(flipped_py1, "Sz", 0)
sz_v3_after = ic_v31.vev("Sz", 0)

print("=== (1) Pauli-X on the exact D=1 field-polarized point ===")
print("<Sz> before -- v3: {}  pyitensor: {}".format(sz_v3_before, sz_py_before))
print("<Sz> after  -- v3: {}  pyitensor: {}  (expect after=-before, v3==pyitensor)".format(
    sz_v3_after, sz_py_after))
assert abs(sz_v3_after - sz_py_after) < 1e-8
assert abs(eta_v31 - flipped_py1.eta) < 1e-8
print("PASSED (v3 and pyitensor agree exactly)")

##############################################################################
### (2) chi_W=1 Pauli-X flip on entangled TFIM VUMPS states, swept over D ###
##############################################################################
Ds = [1, 2, 3, 4]
rs = list(range(1, 6))
sz_before_v3, sz_after_v3, eta_v3_list = [], [], []
sz_before_py, sz_after_py, eta_py_list = [], [], []
corr_before_D2, corr_after_D2_v3, corr_after_D2_py = None, None, None

for D in Ds:
    ic_py = tfim_chain(1.5, D, "python")
    ic_v3 = tfim_chain(1.5, D, 3)

    sz0_py = vumps.onsite_expectation(ic_py._vumps_result, "Sz", 0)
    sz0_v3 = ic_v3.vev("Sz", 0)
    c0 = [vumps.two_point_correlator(ic_py._vumps_result, "Sz", 0, "Sz", r).real for r in rs]

    flipped_py = vumps.apply_mpo(ic_py._vumps_result, W_py, cutoff=1e-12, maxdim=None)
    eta_v3 = apply_pauli_x_v3(ic_v3)

    sz1_py = vumps.onsite_expectation(flipped_py, "Sz", 0)
    sz1_v3 = ic_v3.vev("Sz", 0)
    c1_v3 = [ic_v3.correlator("Sz", 0, "Sz", r).real for r in rs]
    c1_py = [vumps.two_point_correlator(flipped_py, "Sz", 0, "Sz", r).real for r in rs]

    assert abs(sz1_v3 - sz1_py) < 1e-6
    assert abs(eta_v3 - flipped_py.eta) < 1e-6
    for a, b in zip(c1_v3, c1_py):
        assert abs(a - b) < 1e-6

    sz_before_v3.append(sz0_v3); sz_after_v3.append(sz1_v3); eta_v3_list.append(eta_v3.real)
    sz_before_py.append(sz0_py); sz_after_py.append(sz1_py); eta_py_list.append(flipped_py.eta.real)
    if D == 2:
        corr_before_D2, corr_after_D2_v3, corr_after_D2_py = c0, c1_v3, c1_py

print()
print("=== (2) Pauli-X on entangled TFIM VUMPS ground states, D={} ===".format(Ds))
print("<Sz> after (v3): ", sz_after_v3)
print("<Sz> after (pyitensor): ", sz_after_py)
print("PASSED (v3 matches pyitensor at every D)")

########################################################################
### (3) chi_W>1 local gate: a 2-site unitary rotation, one per unit cell ###
########################################################################
def build_ic2(itensor_version):
    ic2 = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=itensor_version)
    ic2.gs_method = "vumps"
    h2 = (ic2.SxC[0] * ic2.SxC[1] + ic2.SyC[0] * ic2.SyC[1] + ic2.SzC[0] * ic2.SzC[1]
          + 0.4 * (ic2.SxC[1] * ic2.SxR[0] + ic2.SyC[1] * ic2.SyR[0] + ic2.SzC[1] * ic2.SzR[0]))
    ic2.set_hamiltonian(h2)
    ic2.maxm, ic2.maxiter, ic2.etol = 3, 120, 1e-10
    ic2.gs_energy()
    return ic2


d2 = 2
H2 = (np.kron(_SX, _SX) + np.kron(_SY, _SY) + np.kron(_SZ, _SZ)).real
gate = expm(-1j * 0.37 * H2)
gate4 = np.transpose(gate.reshape(d2, d2, d2, d2), (2, 0, 3, 1))
U, S, Vh = np.linalg.svd(gate4.reshape(d2 * d2, d2 * d2), full_matrices=False)
keep = int(np.sum(S > 1e-12))
U, S, Vh = U[:, :keep], S[:keep], Vh[:keep, :]
a_half = (U * S[None, :] ** 0.5).reshape(d2, d2, keep)
b_half = (S[:, None] ** 0.5 * Vh).reshape(keep, d2, d2)

ic2_py = build_ic2("python")
ic2_v3 = build_ic2(3)

s0, s1 = ic2_py._vumps_result.sites_uc.si(1), ic2_py._vumps_result.sites_uc.si(2)
left_dummy, mid, right_dummy = Index(1, tags="Link"), Index(keep, tags="Link"), Index(1, tags="Link")
W0_py = ITensor((left_dummy, s0, s0.prime(1), mid), a_half.reshape(1, d2, d2, keep))
W1_py = ITensor((mid, s1, s1.prime(1), right_dummy), b_half.reshape(keep, d2, d2, 1))
gated_py = vumps.apply_mpo(ic2_py._vumps_result, [W0_py, W1_py], cutoff=1e-10, maxdim=None)

W0_v3 = a_half.reshape(1, d2, d2, keep).flatten().tolist()
W1_v3 = b_half.reshape(keep, d2, d2, 1).flatten().tolist()
D_v3, d_g_v3, AL_v3, AR_v3, C_v3, AC_v3, eta_v3_gate = ic2_v3._session3.vumps_apply_mpo(
    [W0_v3, W1_v3], [1, keep], [keep, 1], 1e-10, 0)

print()
print("=== (3) 2-site unitary gate applied once per unit cell (n_uc=2) ===")
print("original bond dim D (pyitensor/v3):", ic2_py._vumps_result.D, ic2_v3._session3.vumps_get_snapshot()[0])
print("gated bond dim D    (pyitensor/v3):", gated_py.D, D_v3)
print("norm diagnostic eta (pyitensor/v3, unitary -> should be ~1):", gated_py.eta, eta_v3_gate)
assert abs(gated_py.eta - 1.0) < 1e-4
assert abs(eta_v3_gate - 1.0) < 1e-4
print("PASSED")

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

axes[0].plot(Ds, sz_after_v3, "o-", color="tab:blue", label="v3")
axes[0].plot(Ds, sz_after_py, "x--", color="tab:red", label="pyitensor")
axes[0].set_xlabel("VUMPS bond dimension D")
axes[0].set_ylabel(r"$\langle S_z\rangle$ after Pauli-X")
axes[0].set_title("(2) TFIM($g=1.5$): flipped $\\langle S_z\\rangle$ vs $D$")
axes[0].legend()

axes[1].plot(rs, corr_before_D2, "o-", color="black", label="before Pauli-X")
axes[1].plot(rs, corr_after_D2_v3, "s--", color="tab:blue", label="after (v3)")
axes[1].plot(rs, corr_after_D2_py, "x:", color="tab:red", label="after (pyitensor)")
axes[1].set_xlabel("distance $r$")
axes[1].set_ylabel(r"$\langle S_z(0)S_z(r)\rangle$")
axes[1].set_title("(2) chi_W=1 at D=2")
axes[1].legend()

labels = ["orig\n(pyitensor)", "orig\n(v3)", "gated\n(pyitensor)", "gated\n(v3)"]
values = [ic2_py._vumps_result.D, ic2_v3._session3.vumps_get_snapshot()[0], gated_py.D, D_v3]
axes[2].bar(labels, values, color=["tab:red", "tab:blue", "tab:red", "tab:blue"])
axes[2].set_ylabel("VUMPS bond dimension D")
axes[2].set_title("(3) chi_W>1 2-site gate, n_uc=2\neta: pyitensor={:.6f}, v3={:.6f}".format(
    gated_py.eta.real, eta_v3_gate.real))

fig.suptitle("Chain::vumps_apply_mpo (v3) VS pyitensor.vumps.apply_mpo: bounded-MPO application")
fig.tight_layout()
fig.savefig("vumps_apply_mpo_v3_VS_python.png", dpi=150)
print("saved plot to vumps_apply_mpo_v3_VS_python.png")
