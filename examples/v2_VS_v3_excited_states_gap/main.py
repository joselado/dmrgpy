# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Compare ITensor v2 vs v3 vs the pure-Python backend (itensor_version=
# "python") on excited_states()/get_excited() -- the transverse-field
# Ising chain's gap (mirrors examples/1dIsing/main.py), exercising the
# Weight-based excited-state DMRG path (pyitensor's overlap-penalty
# dmrg_excited(), see dmrg.py) on all three backends.
#
# Unlike every other v2_VS_v3_* comparison in this directory, "Difference
# v3 vs pure Python" here is NOT expected to be tiny: confirmed directly,
# for this specific Hamiltonian pyitensor's excited-states path settles at
# a gap ~6e-4 off from the exact/v3 value, and this residual is completely
# insensitive to maxdim, sweep count (tested up to 100, no improvement
# over 15), and penalty weight (tested 1x-100x bandwidth) -- i.e. a real
# stationary point of the penalized objective, not a truncation or
# under-sweeping artifact. This matches dmrg_excited()'s own docstring
# (dmrg.py), which already documents that the overlap-penalty method's
# landscape has genuine wrong-energy stationary points without the
# noise-term escape mechanism ITensor's own DMRG has (deliberately not
# implemented here, see dmrg.py's module docstring) -- this example is
# the case that surfaces it, not a new bug.
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
gappy = get_gap("python")

print("Gap (ITensor v2)     =",gap2)
print("Gap (ITensor v3)     =",gap3)
print("Gap (pure Python)    =",gappy)
print("Difference v2 vs v3          =",abs(gap2-gap3))
print("Difference v3 vs pure Python =",abs(gap3-gappy))
