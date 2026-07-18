# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Regression test: ITensor v2 and v3 must keep agreeing on the ground
# state energy of a Heisenberg spin-1/2 chain at larger system sizes
# (n>=20), not just the small n=8 case covered by
# examples/v2_VS_v3_heisenberg_spin_half. ED is not feasible at these
# sizes, so this checks DMRG-vs-DMRG cross-backend consistency instead.
# Added after directly comparing dmrgpy's current in-process v2/v3
# backends against the old (pre-rewrite) file-based v2 backend at
# n=20..32 and finding agreement to ~1e-6..1e-14 (normal DMRG variational
# noise from differing internal sweep paths, not a regression) -- this
# test locks in that same v2-vs-v3 agreement going forward.
import numpy as np
from dmrgpy import spinchain

sizes = [16, 24, 32]

def get_energy(itensor_version, n):
    spins = ["S=1/2" for i in range(n)]
    sc = spinchain.Spin_Chain(spins, itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc.gs_energy(mode="DMRG")

tol = 1e-3 # generous: only meant to catch a real divergence, not to pin
           # down DMRG's sweep-to-sweep numerical noise
for n in sizes:
    e2 = get_energy(2, n)
    e3 = get_energy(3, n)
    diff = abs(e2-e3)
    print("n=%d  v2=%.10f  v3=%.10f  diff=%.2e"%(n, e2, e3, diff))
    assert diff<tol, "v2 vs v3 disagree by %g at n=%d (tol=%g)"%(diff, n, tol)

print("TEST PASSED")
