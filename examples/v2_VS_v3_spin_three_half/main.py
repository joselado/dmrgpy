# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Compare ITensor v2 vs v3 vs the pure-Python backend (itensor_version=
# "python") for a spin-3/2 chain: this is one of dmrgpy's own custom
# (non-stock-ITensor) site types (extra/spinthreehalf.h in both mpscpp2
# and mpscpp3, SpinThreeHalfSite in pyitensor/sites/spin.py), so this also
# checks that both the ITensor-v3 port and the from-scratch pyitensor port
# of that custom site type reproduce the same physics. Only Sz/Sx/Sy are
# defined for this representation (no Sp/Sm), unlike spin-1/2 or spin-1.
import numpy as np
from dmrgpy import spinchain

n = 4
spins = ["S=3/2" for i in range(n)]

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
epy = get_energy("python")

print("Ground state energy (ITensor v2)     =",e2)
print("Ground state energy (ITensor v3)     =",e3)
print("Ground state energy (pure Python)    =",epy)
print("Difference v2 vs v3                  =",abs(e2-e3))
print("Difference v3 vs pure Python         =",abs(e3-epy))
