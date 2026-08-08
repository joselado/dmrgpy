# Direct sum of two converged VUMPS iMPS -- pyitensor.vumps.imps_sum, the
# VUMPS-mixed-gauge analogue of pyitensor.idmrg.imps_sum (see
# examples/idmrg/imps_sum/main.py). Block-diagonal in the bond space (VUMPS
# already works at a single grouped-supersite level, so there is only ever
# one cut to sum across, unlike idmrg's own per-sublattice periodic
# construction), re-canonicalized/truncated via idmrg._canonicalize_periodic
# and completed to the full mixed gauge {AL,AR,C,AC} via
# vumps._complete_mixed_gauge (the standard "bringing a uniform MPS to
# canonical form" construction, Vanderstraeten, Haegeman, Verstraete,
# arXiv:1810.07006, Eq.(9)-(17)).
#
# SCOPE: identical to idmrg.imps_sum's own -- tiled to the thermodynamic
# limit, this "+" is only well-posed (a single dominant branch surviving,
# the correct infinite-volume answer) when the two summands have a genuine
# per-site norm mismatch (their own self-overlap transfer eigenvalues eta
# differ). Two *ordinary* converged VUMPSResults are always individually
# normalized to eta=1 exactly on BOTH the left and right transfer
# eigenvalue (mixed-gauge construction), so summing two of them reliably
# ties -- a "cat state" superposition of two macroscopically distinct
# branches that this module's single-fixed-point machinery cannot
# represent -- and vumps.imps_sum raises RuntimeError there rather than
# silently returning one arbitrary branch.
#
# Demonstrated below: (1) the common, out-of-scope case -- summing two
# independently-converged, oppositely field-polarized VUMPS ground states
# raises RuntimeError, not a silently wrong answer; (2) the well-posed
# case, at a genuinely entangled bond dimension D -- summing a TFIM ground
# state with a deliberately norm-rescaled copy of a *different* TFIM
# ground state (different field g) reproduces the larger-norm branch's own
# onsite magnetization and two-point correlator exactly, with the
# smaller-norm branch fully (not just partially) discarded; (3) a cheap
# parameter sweep over the dominant branch's own field g, confirming
# imps_sum's reconstruction matches the untouched branch across the sweep,
# not just at one accidentally-favorable point.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps

D = 2
RS = list(range(1, 5))


def tfim_chain(g, maxm=D):
    """H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i] -- the standard Pauli-matrix
    TFIM convention, g>1: gapped paramagnetic phase (VUMPS converges
    reliably here)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = maxm
    ic.etol = 1e-11
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


def field_chain(B):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(-B * ic.SzC[0])
    ic.gs_energy()
    return ic


def rescaled_copy(result, factor):
    """A UniformMPS identical to `result` but with every mixed-gauge
    tensor rescaled by `factor` (self-overlap eta -> factor**2) -- ordinary
    VUMPSResults never carry eta!=1 on their own, so this is purely to
    exercise imps_sum's well-posed branch explicitly."""
    return vumps.UniformMPS(
        result.sites_uc, result.n_uc, result.D, result.d_g,
        result.AL * factor, result.AR * factor,
        result.C * factor, result.AC * factor ** 2,
        eta=factor ** 2)


################################################################
### (1) the common, degenerate case: two ordinary VUMPSResults ###
################################################################
ic_up = field_chain(10.0)
ic_down = field_chain(-10.0)

print("=== summing two ordinary (both eta=1) converged VUMPS states ===")
try:
    vumps.imps_sum(ic_up._vumps_result, ic_down._vumps_result)
    raise SystemExit("expected RuntimeError, but imps_sum returned normally")
except RuntimeError as e:
    print("raised RuntimeError as expected (degenerate dominant eigenvalue):")
    print(" ", str(e)[:120], "...")

##########################################################################
### (2) the well-posed case: sum against a deliberately norm-rescaled state ###
##########################################################################
ic_dom = tfim_chain(1.5)
ic_other = tfim_chain(2.0)
scaled_other = rescaled_copy(ic_other._vumps_result, 0.9)

result = vumps.imps_sum(ic_dom._vumps_result, scaled_other, cutoff=1e-12, maxdim=None)

print()
print("=== summing the (unscaled) g=1.5 TFIM state with a 0.9x-rescaled g=2.0 state ===")
print("result.eta (expect ~1, the unscaled/larger-norm branch's own value):", result.eta)
sz_dom = vumps.onsite_expectation(ic_dom._vumps_result, "Sz", 0)
sz_result = vumps.onsite_expectation(result, "Sz", 0)
print("<Sz> of the original g=1.5 state:", sz_dom)
print("<Sz> of the summed-and-truncated result:", sz_result, " (expect equal)")
assert abs(result.eta - 1.0) < 1e-6
assert abs(sz_result - sz_dom) < 1e-6
print("bond dim of original:", ic_dom._vumps_result.D, " of result:", result.D,
      " (expect equal -- the smaller-norm g=2.0 branch is fully discarded)")
assert result.D == ic_dom._vumps_result.D

corr_dom = [vumps.two_point_correlator(ic_dom._vumps_result, "Sz", 0, "Sz", r).real for r in RS]
corr_result = [vumps.two_point_correlator(result, "Sz", 0, "Sz", r).real for r in RS]
for r, cd, cr in zip(RS, corr_dom, corr_result):
    assert abs(cd - cr) < 1e-6, (r, cd, cr)

print()
print("imps_sum well-posed case PASSED")

#######################################################################
### (3) sweep the dominant branch's own field g -- cheap parameter axis ###
#######################################################################
G_VALUES = [1.2, 1.5, 2.0, 3.0]
sz_dom_sweep, sz_result_sweep = [], []
for g in G_VALUES:
    ic_g = tfim_chain(g)
    scaled = rescaled_copy(ic_other._vumps_result, 0.9)
    res_g = vumps.imps_sum(ic_g._vumps_result, scaled, cutoff=1e-12, maxdim=None)
    sz_dom_sweep.append(vumps.onsite_expectation(ic_g._vumps_result, "Sz", 0).real)
    sz_result_sweep.append(vumps.onsite_expectation(res_g, "Sz", 0).real)
    assert abs(sz_dom_sweep[-1] - sz_result_sweep[-1]) < 1e-6

print()
print("imps_sum field sweep PASSED (matches at every g in", G_VALUES, ")")

fig, axs = plt.subplots(1, 3, figsize=(13, 4.5))

labels = ["original\ng=1.5 TFIM", "summed +\ntruncated"]
axs[0].bar(labels, [sz_dom, sz_result], color=["tab:blue", "tab:orange"])
axs[0].set_ylabel(r"$\langle S_z\rangle$")
axs[0].set_title("<Sz> preserved by imps_sum\n(degenerate case: RuntimeError, nothing to plot)")

axs[1].plot(RS, corr_dom, "o-", color="tab:blue", label="original g=1.5")
axs[1].plot(RS, corr_result, "x--", color="tab:orange", label="summed + truncated")
axs[1].set_xlabel(r"separation $r$")
axs[1].set_ylabel(r"$\langle S_z(0) S_z(r)\rangle$")
axs[1].set_title("two-point correlator preserved")
axs[1].legend()

axs[2].plot(G_VALUES, sz_dom_sweep, "o-", color="tab:blue", label="original branch")
axs[2].plot(G_VALUES, sz_result_sweep, "x--", color="tab:orange", label="summed + truncated")
axs[2].set_xlabel(r"dominant branch field $g$")
axs[2].set_ylabel(r"$\langle S_z\rangle$")
axs[2].set_title("agreement across a field sweep")
axs[2].legend()

fig.suptitle("vumps.imps_sum: well-posed (norm-mismatch) case")
fig.tight_layout()
fig.savefig("vumps_imps_sum.png", dpi=150)
print("saved plot to vumps_imps_sum.png")
