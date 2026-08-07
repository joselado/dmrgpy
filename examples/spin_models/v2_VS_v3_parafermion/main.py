# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test: ground state energy of a Z3 parafermion chain must
# agree across ITensor v2, v3, and ED. There was previously no
# v2_VS_v3_* comparison for the parafermion model at all (unlike spin,
# fermion, boson chains), and ED mode used to crash outright here with
# "KeyError: ('Id', 1)" -- multioperator.py's MO2matrix (the ED
# evaluation path) didn't filter the zero-coefficient "Id" placeholder
# term that "0 + MultiOperator" (the pervasive "h = 0; h = h + term"
# idiom) creates via __radd__, so it tried (and failed) to look up an
# operator ED's parafermion backend never registered. Fixed by filtering
# near-zero-coefficient terms in MO2matrix, matching what to_terms()
# already does for the DMRG path -- this test locks that fix in.
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import parafermionchain

n = 6

def get_energy(itensor_version, n, mode="DMRG"):
    pc = parafermionchain.Parafermionic_Chain(n)
    # Switch backend *before* seeding: constructing the pure-Python backend
    # (setup_python()) consumes draws from numpy's global RNG as a side
    # effect (its internal random-state initialization), unlike v2/v3 which
    # never touch np.random. Seeding before this point would desync the
    # "same" Hamiltonian across backends -- confirmed directly: an earlier
    # version of this test seeded first and saw the python backend land on
    # a wildly different (and wrongly suspected buggy) energy purely
    # because it solved a different random Hamiltonian, not because of any
    # real physics bug.
    if itensor_version!="python": pc.setup_cpp(itensor_version)
    else: pc.setup_python()
    np.random.seed(3) # same Hamiltonian for every backend
    h = 0
    for i in range(n):
        for j in range(n):
            h = h + pc.N[i]*pc.N[j]*np.random.random()
            h = h + pc.Sig[i]*pc.Sigd[j]*np.random.random()
            h = h + pc.Tau[i]*pc.Taud[j]*np.random.random()
    h = h + h.get_dagger() # Hermitize, this is what creates the stray
                            # zero "Id" term the fix above addresses
    pc.set_hamiltonian(h)
    return pc.gs_energy(mode=mode)

e2 = get_energy(2, n)
e3 = get_energy(3, n)
eed = get_energy(2, n, mode="ED") # was: KeyError('Id', 1) before the fix
epy = get_energy("python", n) # n=6 < 10, python backend is in scope

print("Ground state energy (ITensor v2)  =",e2)
print("Ground state energy (ITensor v3)  =",e3)
print("Ground state energy (ED)          =",eed)
print("Ground state energy (pure Python) =",epy)

tol = 1e-2
for name,e in [("v3",e3),("ED",eed),("python",epy)]:
    diff = abs(e2-e)
    print("Difference v2 vs %s = %.2e"%(name,diff))
    assert diff<tol, "v2 vs %s disagree by %g (tol=%g)"%(name,diff,tol)

print("TEST PASSED")

# sweep the chain length and overlay all four solvers -- they should all
# agree at every size
ns = [3,4,5,6]
e2s = [get_energy(2, ni) for ni in ns]
e3s = [get_energy(3, ni) for ni in ns]
eeds = [get_energy(2, ni, mode="ED") for ni in ns]
epys = [get_energy("python", ni) for ni in ns]

plt.plot(ns,e2s,marker="o",label="ITensor v2")
plt.plot(ns,e3s,marker="s",label="ITensor v3")
plt.plot(ns,eeds,marker="^",label="ED")
plt.plot(ns,epys,marker="x",label="pure Python")
plt.xlabel("Chain length")
plt.ylabel("Ground state energy")
plt.legend()
plt.show()
