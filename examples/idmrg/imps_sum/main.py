# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Sum of two infinite MPS with idmrg.imps_sum. In the thermodynamic limit
# the sum is dominated by the branch with the larger norm per site, and the
# other branch is discarded. We sum a Heisenberg ground state with an XXZ ground state
# rescaled by 0.9, and check that we get back the Heisenberg state.
# Both states must be well converged, or the sum cannot be gauge fixed.
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor.tensor import ITensor

def ground_state(delta): # XXZ chain with anisotropy delta
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "idmrg"
    ic.set_hamiltonian(ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0]
            + delta*ic.SzC[0]*ic.SzR[0])
    ic.maxm, ic.maxiter, ic.etol = 16, 100, 1e-12 # bond dimension, iterations
    ic.gs_energy()
    return ic._result

a = ground_state(1.0) # Heisenberg chain
b = ground_state(0.3) # XXZ chain
cell, _ = idmrg._correlator_cell(b) # tensors of the repeated cell
b = idmrg.PeriodicMPS(b.sites_uc, 1, [ITensor(T.inds, 0.9*T.array)
        for T in cell], eta=0.81) # norm per site 0.81
c = idmrg.imps_sum(a, b, cutoff=1e-12, maxdim=None)

chi = [s.U_list[0].array.shape[-1] for s in (a, b, c)] # bond dimensions
zz = [idmrg.two_point_correlator(s, "Sz", 0, "Sz", 1).real for s in (a, b, c)]
print("Bond dimension of the Heisenberg, XXZ and summed states:", chi)
print("<Sz_0 Sz_1> of the three states:", zz)
print("Max difference =", abs(zz[2]-zz[0]))
assert chi[2] == chi[0] and abs(zz[2]-zz[0]) < 1e-6

import matplotlib.pyplot as plt
plt.bar(["Heisenberg", "XXZ", "sum"], chi)
plt.ylabel("bond dimension")
plt.show()
