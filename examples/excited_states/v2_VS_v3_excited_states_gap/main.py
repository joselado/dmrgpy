# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Compare ITensor v2 vs v3 vs the pure-Python backend (itensor_version=
# "python") on excited_states()/get_excited() -- the transverse-field
# Ising chain's gap (mirrors examples/1dIsing/main.py), exercising the
# Weight-based excited-state DMRG path (pyitensor's overlap-penalty
# dmrg_excited(), see dmrg.py) on all three backends.
#
# All three backends agree to ~1e-14 here.
#
# This comment used to say the opposite, and the correction is worth
# keeping. It read: pyitensor's gap "settles ~6e-4 off from the exact/v3
# value", insensitive to maxdim, sweep count and penalty weight, and
# therefore "a real stationary point of the penalized objective, not a
# truncation or under-sweeping artifact ... not a new bug". Every
# measurement in that claim was accurate. The conclusion drawn from it was
# wrong: the insensitivity to every convergence knob was not evidence of a
# hard stationary point, it was evidence of a plain bug, which no amount
# of sweeping can fix.
#
# The cause was a double conjugation in dmrg_excited()'s overlap-penalty
# projector (pyitensor/dmrg.py::_bond_projections -- see its docstring for
# the derivation): the penalty projected onto the wrong direction whenever
# the states were complex, so it barely orthogonalized at all. It went
# unnoticed because it is exactly invisible for real-valued tensors, and
# because nothing in tests/ covered pyitensor's excited states against ED.
# It is fixed, and tests/test_excited_states_pyitensor.py now pins it.
#
# The lesson this example is kept for: "the error does not respond to any
# convergence parameter" is a reason to suspect a bug, not a reason to
# conclude the algorithm has converged to something.
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

# Regression guard for the fix described at the top of this file: before
# it, this difference was ~6e-4.
assert abs(gap3-gappy) < 1e-8
assert abs(gap2-gap3) < 1e-8

import matplotlib.pyplot as plt
labels = ["ITensor v2","ITensor v3","pure Python"]
gaps = [gap2,gap3,gappy]
plt.bar(labels,gaps)
plt.ylabel("Excited-state gap")
plt.title("Transverse-field Ising gap: backend comparison")
plt.show()
