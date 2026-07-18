# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Compare ITensor v2 vs v3 vs the pure-Python backend (itensor_version=
# "python") on the Kernel Polynomial Method dynamical correlator path
# (Chain::kpm_dynamical_correlator in chain_session.h; pyfermion/algebra/
# kpm.py-style moment recursion in pyitensor -- see kpmdmrg.py).
import numpy as np
from dmrgpy import spinchain

n = 6

def get_correlator(itensor_version):
    sc = spinchain.Spin_Chain([2 for i in range(n)],itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.get_gs()
    name = (sc.Sz[0],sc.Sz[0]) # local dynamical correlator <Sz_0(t)Sz_0(0)>
    (x,y) = sc.get_dynamical_correlator(mode="DMRG",submode="KPM",name=name)
    return np.array(x),np.array(y)

x2,y2 = get_correlator(2)
x3,y3 = get_correlator(3)
xpy,ypy = get_correlator("python")

print("Sample of the correlator (ITensor v2)  :",y2[:5])
print("Sample of the correlator (ITensor v3)  :",y3[:5])
print("Sample of the correlator (pure Python) :",ypy[:5])
# v3 and the python backend always return the same-length grid (both are
# real DMRG runs with the same nominal settings); print this comparison
# first since it's the one this example is actually about.
print("Max abs difference v3 vs pure Python =",np.max(np.abs(y3-ypy)))
# v2's ED fallback (when the mpscpp2 extension isn't compiled -- see
# mode.py) can return a differently-sized frequency grid than v3's own
# KPM path; not related to the python backend, so just guard against it
# instead of crashing the whole script.
if len(y2) == len(y3):
    print("Max abs difference v2 vs v3          =",np.max(np.abs(y2-y3)))
else:
    print("Max abs difference v2 vs v3          = N/A (different grid sizes: "
          "{} vs {} points -- likely v2's ED fallback, see mode.py)".format(len(y2),len(y3)))
