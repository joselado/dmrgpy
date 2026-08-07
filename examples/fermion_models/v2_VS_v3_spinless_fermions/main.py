# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Compare ITensor v2 vs v3 vs the pure-Python backend (itensor_version=
# "python") for a spinless fermion hopping chain: ground state energy plus
# the full <Cdag_i C_j> correlation matrix (exercises fermionic operators/
# Jordan-Wigner strings and the correlation_matrix path on all three
# backends).
import numpy as np
from dmrgpy import fermionchain

n = 6

def get_result(itensor_version,n):
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

e3,cm3 = get_result(3,n)
epy,cmpy = get_result("python",n)

print("Ground state energy (ITensor v3)     =",e3)
print("Ground state energy (pure Python)    =",epy)
print("Energy difference v3 vs pure Python  =",abs(e3-epy))
print("Correlation matrix max abs diff v3 vs pure Python =",np.max(np.abs(cm3-cmpy)))

# v2's ED fallback (when the mpscpp2 extension isn't compiled -- see
# mode.py) returns a plain ED-backend Fermionic_State that never grew a
# get_correlation_matrix() method (pyfermion/mbfermion.py); not related to
# the python backend, so just skip it gracefully instead of crashing the
# whole script.
v2_available = True
try:
    e2,cm2 = get_result(2,n)
    print("Ground state energy (ITensor v2)     =",e2)
    print("Energy difference v2 vs v3           =",abs(e2-e3))
    print("Correlation matrix max abs diff v2 vs v3          =",np.max(np.abs(cm2-cm3)))
except AttributeError as e:
    v2_available = False
    print("ITensor v2 comparison skipped ({}) -- likely v2's ED fallback, "
          "see mode.py".format(e))

# sweep the chain length (cheap for this spinless hopping chain), comparing
# ground-state energies across backends
ns = [3,4,5,6,7]
e3s = [get_result(3,ni)[0] for ni in ns]
epys = [get_result("python",ni)[0] for ni in ns]
e2s = None
if v2_available:
    e2s = [get_result(2,ni)[0] for ni in ns]
for i,ni in enumerate(ns):
    line = "n = "+str(ni)+" v3 = "+str(e3s[i])+" python = "+str(epys[i])
    if e2s is not None:
        line += " v2 = "+str(e2s[i])
    print(line)

import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(ns,e3s,marker="x",linestyle="--",label="ITensor v3")
ax1.plot(ns,epys,marker="s",linestyle=":",label="pure Python")
if e2s is not None:
    ax1.plot(ns,e2s,marker="o",label="ITensor v2")
ax1.set_xlabel("Number of sites")
ax1.set_ylabel("Ground-state energy")
ax1.legend()
sites = list(range(n))
ax2.plot(sites,cm3[0,:].real,marker="x",linestyle="--",label="ITensor v3")
ax2.plot(sites,cmpy[0,:].real,marker="s",linestyle=":",label="pure Python")
ax2.set_xlabel("Site j")
ax2.set_ylabel(r"Re$\langle C^\dagger_0 C_j\rangle$")
ax2.legend()
plt.tight_layout()
plt.show()
