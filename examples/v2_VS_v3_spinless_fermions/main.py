# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Compare ITensor v2 vs v3 for a spinless fermion hopping chain: ground
# state energy plus the full <Cdag_i C_j> correlation matrix (exercises
# fermionic operators/Jordan-Wigner strings and the correlation_matrix
# path on both backends).
import numpy as np
from dmrgpy import fermionchain

n = 6

def get_result(itensor_version):
    fc = fermionchain.Fermionic_Chain(n,itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1]
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)
    e = fc.gs_energy()
    wf = fc.get_gs()
    cm = wf.get_correlation_matrix()
    return e,cm

e2,cm2 = get_result(2)
e3,cm3 = get_result(3)

print("Ground state energy (ITensor v2) =",e2)
print("Ground state energy (ITensor v3) =",e3)
print("Energy difference =",abs(e2-e3))
print("Correlation matrix max abs difference =",np.max(np.abs(cm2-cm3)))
