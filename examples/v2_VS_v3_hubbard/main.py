# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Compare ITensor v2 vs v3 for a (spinful) Hubbard chain, built the same
# way as examples/hubbard_gap/main.py: two interleaved spinless-fermion
# sites (up/down) per physical site, plus an on-site U(n_up-1/2)(n_dn-1/2)
# interaction.
import numpy as np
from dmrgpy import fermionchain

ns = 3 # physical sites (2*ns fermionic sites total)
U = 4.0

def get_energy(itensor_version):
    n = ns*2
    fc = fermionchain.Fermionic_Chain(n,itensor_version=itensor_version)
    C = fc.C ; Cdag = fc.Cdag ; N = fc.N
    h = 0
    for i in range(ns-1):
        for j in range(2):
            h = h + Cdag[2*i+j]*C[2*(i+1)+j]
    h = h + h.get_dagger()
    for i in range(ns):
        h = h + U*(N[2*i]-0.5)*(N[2*i+1]-0.5)
    fc.set_hamiltonian(h)
    return fc.gs_energy()

e2 = get_energy(2)
e3 = get_energy(3)

print("Ground state energy (ITensor v2) =",e2)
print("Ground state energy (ITensor v3) =",e3)
print("Difference =",abs(e2-e3))
