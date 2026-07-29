# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

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

j_strong, j_weak = 1.0, 0.4
ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
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
for n_sites in (12, 14, 16):
    sc = spinchain.Spin_Chain(["1/2"]*n_sites)
    h = 0
    for i in range(n_sites - 1):
        j = j_strong if i % 2 == 0 else j_weak
        h = h + j*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    gap = sc.get_gap(mode="ED")
    print("  n_sites={}  ED gap={:.6f}".format(n_sites, gap))
