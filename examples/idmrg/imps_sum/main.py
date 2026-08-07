# Direct sum of two converged infinite MPS -- pyitensor.idmrg.imps_sum,
# the periodic-chain analogue of mpsalgebra.sum: block-diagonal in the
# bond space at every cut (there is no open boundary to instead
# concatenate along, unlike a finite chain), re-canonicalized/truncated
# via the same two-sided fixed-point procedure idmrg.apply_mpo already
# uses (idmrg._canonicalize_periodic).
#
# SCOPE: tiled to the thermodynamic limit, this "+" is only well-posed
# (a single dominant branch surviving, the correct infinite-volume
# answer) when the two summands have a genuine per-site norm mismatch
# (their own self-overlap transfer eigenvalues eta differ). Two
# *ordinary* converged IDMRGResults are always individually normalized to
# eta=1 exactly, so summing two of them (the most natural reason to want
# this operation) reliably ties -- a "cat state" superposition of two
# macroscopically distinct branches that this module's single-fixed-point
# machinery cannot represent, and idmrg.imps_sum raises RuntimeError
# there rather than silently returning one arbitrary branch. See
# pyitensor/idmrg.py's own "Summing two converged iMPS" section docstring
# for the full derivation.
#
# Demonstrated below: (1) the common, out-of-scope case -- summing two
# independently-converged, oppositely Sz-polarized product states raises
# RuntimeError, not a silently wrong answer; (2) the well-posed case --
# summing a converged state with a deliberately norm-rescaled copy of a
# *different* converged state reproduces the larger-norm branch's own
# observables exactly, with the smaller-norm branch fully (not just
# partially) discarded.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor.tensor import ITensor
import matplotlib.pyplot as plt

################################################################
### (1) the common, degenerate case: two ordinary IDMRGResults ###
################################################################
ic_up = infinitechain.Infinite_Spin_Chain(["1/2"])
ic_up.maxm, ic_up.maxiter, ic_up.etol = 4, 50, 1e-12
ic_up.set_hamiltonian(-10.0 * ic_up.SzC[0])
ic_up.gs_energy()

ic_down = infinitechain.Infinite_Spin_Chain(["1/2"])
ic_down.maxm, ic_down.maxiter, ic_down.etol = 4, 50, 1e-12
ic_down.set_hamiltonian(10.0 * ic_down.SzC[0])
ic_down.gs_energy()

print("=== summing two ordinary (both eta=1) converged states ===")
try:
    idmrg.imps_sum(ic_up._result, ic_down._result)
    raise SystemExit("expected RuntimeError, but imps_sum returned normally")
except RuntimeError as e:
    print("raised RuntimeError as expected (degenerate dominant eigenvalue):")
    print(" ", str(e)[:120], "...")

##########################################################################
### (2) the well-posed case: sum against a deliberately norm-rescaled state ###
##########################################################################
ic_heis = infinitechain.Infinite_Spin_Chain(["1/2"])
h = ic_heis.SxC[0] * ic_heis.SxR[0] + ic_heis.SyC[0] * ic_heis.SyR[0] + ic_heis.SzC[0] * ic_heis.SzR[0]
ic_heis.set_hamiltonian(h)
ic_heis.maxm, ic_heis.maxiter, ic_heis.etol = 16, 30, 1e-9
ic_heis.gs_energy()

# A different (XXZ-anisotropic) converged ground state, then rescaled by
# 0.9 in amplitude (self-overlap eta=0.81) -- ordinary IDMRGResults never
# carry eta!=1 on their own, this is purely to exercise imps_sum's
# well-posed branch explicitly.
ic_xxz = infinitechain.Infinite_Spin_Chain(["1/2"])
h_xxz = ic_xxz.SxC[0] * ic_xxz.SxR[0] + ic_xxz.SyC[0] * ic_xxz.SyR[0] + 0.3 * ic_xxz.SzC[0] * ic_xxz.SzR[0]
ic_xxz.set_hamiltonian(h_xxz)
ic_xxz.maxm, ic_xxz.maxiter, ic_xxz.etol = 16, 30, 1e-9
ic_xxz.gs_energy()

scaled_xxz = idmrg.PeriodicMPS(
    ic_xxz._result.sites_uc, 1,
    [ITensor(T.inds, T.array * 0.9) for T in ic_xxz._result.U_list], eta=0.81)

result = idmrg.imps_sum(ic_heis._result, scaled_xxz, cutoff=1e-12, maxdim=None)

print()
print("=== summing the (unscaled) Heisenberg state with a 0.9x-rescaled XXZ state ===")
print("result.eta (expect ~1, the unscaled/larger-norm branch's own value):", result.eta)
sz_heis = idmrg.onsite_expectation(ic_heis._result, "Sz", 0)
sz_result = idmrg.onsite_expectation(result, "Sz", 0)
print("<Sz> of the original Heisenberg state:", sz_heis)
print("<Sz> of the summed-and-truncated result:", sz_result, " (expect equal)")
assert abs(result.eta - 1.0) < 1e-6
assert abs(sz_result - sz_heis) < 1e-6
chi_heis = ic_heis._result.U_list[0].array.shape[-1]
chi_result = result.U_list[0].array.shape[-1]
print("bond dim of original:", chi_heis, " of result:", chi_result,
      " (expect equal -- the smaller-norm XXZ branch is fully discarded)")
assert chi_result == chi_heis

print()
print("imps_sum example PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4.5))

labels = ["original\nHeisenberg", "summed +\ntruncated"]
ax1.bar(labels, [sz_heis, sz_result], color=["tab:blue", "tab:orange"])
ax1.set_ylabel(r"$\langle S_z\rangle$")
ax1.set_title("<Sz> preserved by imps_sum")

ax2.bar(labels, [chi_heis, chi_result], color=["tab:blue", "tab:orange"])
ax2.set_ylabel("bond dimension")
ax2.set_title("bond dim preserved\n(smaller-norm branch discarded)")

fig.suptitle("imps_sum: well-posed case (case (1), degenerate eta=1 sum,\n"
              "correctly raises RuntimeError instead -- nothing to plot there)")
fig.tight_layout()
fig.savefig("imps_sum.png", dpi=150)
print("saved plot to imps_sum.png")
