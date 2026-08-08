# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import matplotlib.pyplot as plt
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object
from dmrgpy import spinchain      # finite chain, for the ED cross-check

###########################################################################
### "local superblock gap": a cheap, cruder alternative to             ###
### excitation_gap's tangent-space/quasiparticle ansatz                 ###
###########################################################################
# Unlike excitation_gap (see examples/idmrg/excitation_gap_xx), this does
# NOT require a product-state-like (D=1) ground state -- it simply
# re-diagonalizes the growing algorithm's own final, converged 2-site
# effective Hamiltonian for its second-lowest eigenvalue, instead of only
# its ground state. See pyitensor.idmrg.local_excitation_gap's own
# docstring for the algorithm (and why this is the exact, infinite-chain
# analogue of finite-DMRG's Lagrange-multiplier excited-state trick) and
# its accuracy caveat: this has no momentum label and is not guaranteed to
# match the true minimum-momentum gap the way excitation_gap does.
#
# This example uses a dimerized (alternating strong/weak bond) Heisenberg
# chain -- a genuinely entangled (D>1) ground state, exactly the case
# excitation_gap cannot handle yet -- and cross-checks the result against a
# large finite open chain's own exact-diagonalization gap.
#
# The second half of this script demonstrates the window= refinement
# (local_excitation_gap_windowed): growing the local diagonalization block
# with extra real, free physical sites systematically tightens this
# estimate -- see that section below for a transverse-field Ising example
# (window= only supports n_uc=1 chains).

j_strong, j_weak = 1.0, 0.4
ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
ic.gs_method = "idmrg"
H = (j_strong*(ic.SxC[0]*ic.SxC[1] + ic.SyC[0]*ic.SyC[1] + ic.SzC[0]*ic.SzC[1])
     + j_weak*(ic.SxC[1]*ic.SxR[0] + ic.SyC[1]*ic.SyR[0] + ic.SzC[1]*ic.SzR[0]))
ic.set_hamiltonian(H)
ic.maxm = 40
ic.maxiter = 200
ic.etol = 1e-9

e0 = ic.gs_energy()
print("iDMRG converged:", ic.converged)
print("ground-state energy density:", e0)
print()

local_gap = ic.local_excitation_gap()
print("local_excitation_gap():", local_gap)
print()

print("cross-check: finite open chains of the same dimerized model, ED gap")
print("(get_gap(mode=\"ED\")), at growing size:")
dimer_n_sites = (12, 14, 16)
dimer_ed_gaps = []
for n_sites in dimer_n_sites:
    sc = spinchain.Spin_Chain(["1/2"]*n_sites)
    h = 0
    for i in range(n_sites - 1):
        j = j_strong if i % 2 == 0 else j_weak
        h = h + j*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    gap = sc.get_gap(mode="ED")
    dimer_ed_gaps.append(gap)
    print("  n_sites={}  ED gap={:.6f}".format(n_sites, gap))

###########################################################################
### Tightening it further: window= (idmrg.local_excitation_gap_windowed) ##
###########################################################################
# local_excitation_gap(window=w) grows the local diagonalization block by
# w extra *free* physical sites on each side of the original 2, re-solving
# both the ground state and the deflated first excited state fresh within
# this larger block, instead of only ever using the frozen 2-site block
# above. Only n_uc=1 is supported (see
# pyitensor.idmrg.local_excitation_gap_windowed's own docstring for why),
# so this uses a transverse-field Ising chain instead of the dimerized
# model above -- also genuinely entangled (D>1), also gapped, but with a
# single-site unit cell.
print()
print("=" * 70)
print("window=: transverse-field Ising chain (n_uc=1, D>1, gapped)")
print("H = -4*J*Sx_i*Sx_{i+1} - 2*h*Sz_i  (Pauli convention, sigma=2S)")
print("=" * 70)

J, h_field = 1.0, 2.0
ic2 = infinitechain.Infinite_Spin_Chain(["1/2"])
ic2.gs_method = "idmrg"
H2 = -4*J*ic2.SxC[0]*ic2.SxR[0] - 2*h_field*ic2.SzC[0]
ic2.set_hamiltonian(H2)
ic2.maxm = 12
ic2.maxiter = 200
ic2.etol = 1e-10
ic2.gs_energy()
print("iDMRG converged:", ic2.converged)
print()

print("local_excitation_gap(window=w) as w grows:")
windows = (0, 1, 2, 3)
windowed_gaps = []
for w in windows:
    g = ic2.local_excitation_gap(window=w)
    windowed_gaps.append(g)
    print("  window={}  gap={:.6f}".format(w, g))
print()

print("cross-check: finite open TFIM chains, ED gap (get_gap(mode=\"ED\")),")
print("at growing size -- the windowed estimate above converges toward")
print("this at least as fast as growing the finite chain itself does:")
tfim_n_sites = (12, 14, 16, 18)
tfim_ed_gaps = []
for n_sites in tfim_n_sites:
    sc = spinchain.Spin_Chain(["1/2"]*n_sites)
    h = 0
    for i in range(n_sites - 1):
        h = h + (-4*J)*sc.Sx[i]*sc.Sx[i+1]
    for i in range(n_sites):
        h = h + (-2*h_field)*sc.Sz[i]
    sc.set_hamiltonian(h)
    gap = sc.get_gap(mode="ED")
    tfim_ed_gaps.append(gap)
    print("  n_sites={}  ED gap={:.6f}".format(n_sites, gap))

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

axes[0].plot(dimer_n_sites, dimer_ed_gaps, "o-", color="tab:blue", label="ED gap (finite chain)")
axes[0].axhline(local_gap, color="tab:red", ls="--", label="local_excitation_gap (iDMRG)")
axes[0].set_xlabel("$n_{sites}$")
axes[0].set_ylabel("gap")
axes[0].set_title("dimerized Heisenberg chain")
axes[0].legend()

axes[1].plot(windows, windowed_gaps, "o-", color="tab:green")
axes[1].set_xlabel("window $w$")
axes[1].set_ylabel("local_excitation_gap(window=w)")
axes[1].set_title("TFIM: windowed refinement")

axes[2].plot(tfim_n_sites, tfim_ed_gaps, "o-", color="tab:blue", label="ED gap (finite chain)")
axes[2].axhline(windowed_gaps[-1], color="tab:red", ls="--",
                 label="local_excitation_gap(window={})".format(windows[-1]))
axes[2].set_xlabel("$n_{sites}$")
axes[2].set_ylabel("gap")
axes[2].set_title("TFIM: ED convergence")
axes[2].legend()

fig.suptitle("local_excitation_gap: cheap 2-site estimate vs finite-chain ED gap")
fig.tight_layout()
fig.savefig("local_excitation_gap.png", dpi=150)
print("\nsaved plot to local_excitation_gap.png")
