# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Compare the ITensor v2 (mpscpp2) and v3 (mpscpp3) C++ DMRG backends on
# the simplest possible case: the ground state energy of a spin-1/2
# Heisenberg chain. Both backends are used in the same script/process --
# this works because mpscpp3/bindings.cc registers its pybind11 types with
# py::module_local(), so it never collides with mpscpp2's identically-named
# (but ABI-incompatible) Chain/MPS/MPO types.
import numpy as np
from dmrgpy import spinchain

n = 8 # small chain, fast to converge
spins = ["S=1/2" for i in range(n)]

def get_energy(itensor_version):
    sc = spinchain.Spin_Chain(spins,itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1]
        h = h + sc.Sy[i]*sc.Sy[i+1]
        h = h + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc.gs_energy()

e2 = get_energy(2)
e3 = get_energy(3)

print("Ground state energy (ITensor v2) =",e2)
print("Ground state energy (ITensor v3) =",e3)
print("Difference =",abs(e2-e3))
