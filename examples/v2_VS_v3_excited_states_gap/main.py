# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Compare ITensor v2 vs v3 on excited_states()/get_excited() -- the
# transverse-field Ising chain's gap (mirrors examples/1dIsing/main.py),
# exercising the Weight-based excited-state DMRG path on both backends.
import numpy as np
from dmrgpy import spinchain

n = 8
bx = 0.6 # away from the bx=0/bx=1 critical points, gap is well-defined

def get_gap(itensor_version):
    sc = spinchain.Spin_Chain([2 for i in range(n)],itensor_version=itensor_version)
    h = 0
    for i in range(n-1): h = h + sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): h = h + bx*sc.Sx[i]
    sc.set_hamiltonian(h)
    es = sc.get_excited(mode="DMRG",n=2)
    return es[1]-es[0]

gap2 = get_gap(2)
gap3 = get_gap(3)

print("Gap (ITensor v2) =",gap2)
print("Gap (ITensor v3) =",gap3)
print("Difference =",abs(gap2-gap3))
