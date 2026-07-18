# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Compare ITensor v2 vs v3 vs the pure-Python backend (itensor_version=
# "python") for a spinless fermion hopping chain: ground state energy plus
# the full <Cdag_i C_j> correlation matrix (exercises fermionic operators/
# Jordan-Wigner strings and the correlation_matrix path on all three
# backends).
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

e3,cm3 = get_result(3)
epy,cmpy = get_result("python")

print("Ground state energy (ITensor v3)     =",e3)
print("Ground state energy (pure Python)    =",epy)
print("Energy difference v3 vs pure Python  =",abs(e3-epy))
print("Correlation matrix max abs diff v3 vs pure Python =",np.max(np.abs(cm3-cmpy)))

# v2's ED fallback (when the mpscpp2 extension isn't compiled -- see
# mode.py) returns a plain ED-backend Fermionic_State that never grew a
# get_correlation_matrix() method (pyfermion/mbfermion.py); not related to
# the python backend, so just skip it gracefully instead of crashing the
# whole script.
try:
    e2,cm2 = get_result(2)
    print("Ground state energy (ITensor v2)     =",e2)
    print("Energy difference v2 vs v3           =",abs(e2-e3))
    print("Correlation matrix max abs diff v2 vs v3          =",np.max(np.abs(cm2-cm3)))
except AttributeError as e:
    print("ITensor v2 comparison skipped ({}) -- likely v2's ED fallback, "
          "see mode.py".format(e))
