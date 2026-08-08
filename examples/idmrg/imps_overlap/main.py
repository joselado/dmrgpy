# Per-site overlap/fidelity between two converged infinite MPS --
# pyitensor.idmrg.imps_overlap: the thermodynamic-limit generalization of
# a finite MPS inner product <phi|psi>. A literal <phi|psi> over an
# infinite chain is 0 or infinite (it scales as eta^N over N unit cells)
# unless the dominant mixed-transfer-matrix eigenvalue eta has magnitude
# exactly 1 -- so imps_overlap returns the *per-site fidelity*
# eta_ab/sqrt(eta_aa*eta_bb) by default (normalize=True), which has
# magnitude 1 iff the two iMPS represent the same physical state (any
# gauge/normalization convention), and magnitude < 1 otherwise.
#
# Demonstrated below: (1) self-overlap and overlap with an identity-MPO'd
# copy (idmrg.apply_mpo) are both exactly 1 -- a gauge-independent
# cross-check that apply_mpo's own canonicalization didn't change the
# physical state; (2) two independently-converged, oppositely-polarized
# product states (a strong onsite field pinning <Sz> to +-0.5) are
# orthogonal, |overlap| ~ 0.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor
import numpy as np
import matplotlib.pyplot as plt

################################################################
### (1) self-overlap and gauge-invariance under apply_mpo(Id) ###
################################################################
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "idmrg"
h = ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]
ic.set_hamiltonian(h)
ic.maxm, ic.maxiter, ic.etol = 30, 60, 1e-9
ic.gs_energy()

print("=== self-overlap of a converged Heisenberg ground state ===")
ov_self = idmrg.imps_overlap(ic._result, ic._result)
print("<psi|psi> per site:", ov_self, " (expect exactly 1)")
assert abs(ov_self - 1.0) < 1e-6

sites_uc = ic._result.sites_uc
d = sites_uc.dim(1)
link_l, link_r = Index(1, tags="Link"), Index(1, tags="Link")
s = sites_uc.si(1)
W_id = [ITensor((link_l, s, s.prime(1), link_r), np.eye(d, dtype=complex).reshape(1, d, d, 1))]
same_state = idmrg.apply_mpo(ic._result, W_id, cutoff=1e-12, maxdim=None)

print()
print("=== overlap with an identity-MPO'd (re-gauged) copy ===")
ov_gauge = idmrg.imps_overlap(ic._result, same_state)
print("|<psi|apply_mpo(Id, psi)>| per site:", abs(ov_gauge), " (expect exactly 1)")
assert abs(abs(ov_gauge) - 1.0) < 1e-6

########################################################################
### (2) two independently-converged, oppositely-polarized product states ###
########################################################################
ic_up = infinitechain.Infinite_Spin_Chain(["1/2"])
ic_up.gs_method = "idmrg"
ic_up.maxm, ic_up.maxiter, ic_up.etol = 4, 50, 1e-12
ic_up.set_hamiltonian(-10.0 * ic_up.SzC[0])
ic_up.gs_energy()

ic_down = infinitechain.Infinite_Spin_Chain(["1/2"])
ic_down.gs_method = "idmrg"
ic_down.maxm, ic_down.maxiter, ic_down.etol = 4, 50, 1e-12
ic_down.set_hamiltonian(10.0 * ic_down.SzC[0])
ic_down.gs_energy()

print()
print("=== overlap between opposite-field-polarized product states ===")
print("<Sz> up-chain:", ic_up.vev("Sz", 0), " down-chain:", ic_down.vev("Sz", 0))
ov_orth = idmrg.imps_overlap(ic_up._result, ic_down._result)
print("|<down|up>| per site:", abs(ov_orth), " (expect ~0, orthogonal states)")
assert abs(ov_orth) < 1e-6

print()
print("imps_overlap example PASSED")

labels = ["self-overlap\n(same state)", "gauge-invariance\n(apply_mpo(Id))",
          "orthogonality\n(opposite-polarized)"]
values = [abs(ov_self), abs(ov_gauge), abs(ov_orth)]
fig, ax = plt.subplots(figsize=(6.5, 4.5))
bars = ax.bar(labels, values, color=["tab:blue", "tab:green", "tab:red"])
ax.axhline(1.0, color="gray", lw=0.8, ls=":")
ax.set_ylabel("per-site fidelity |imps_overlap|")
ax.set_title("imps_overlap: self / gauge-invariance / orthogonality checks")
for bar, v in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2, v, "{:.3f}".format(v),
            ha="center", va="bottom")
fig.tight_layout()
fig.savefig("imps_overlap.png", dpi=150)
print("saved plot to imps_overlap.png")
