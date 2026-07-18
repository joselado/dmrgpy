# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Regression test: static two-point correlators <Sz_i Sz_j> (spin chain)
# and <Cdag_i C_j> (fermionic chain, exercises Jordan-Wigner strings)
# must agree between DMRG (v2, v3, pure-Python) and ED, pairwise across
# every (i,j). None of the existing correlator examples
# (examples/correlation_matrix*, examples/*_static_correlator,
# examples/entanglement_entropy/*) check a correlator against ED ground
# truth -- they only print/plot a single backend's result, or compare
# DMRG backends against each other (examples/correlation_matrix_python_VS_CPP).
import numpy as np
from dmrgpy import spinchain, fermionchain

tol = 1e-3

def correlation_matrix(chain, pairop, n, mode, itensor_version=None):
    if itensor_version is not None:
        if itensor_version!="python": chain.setup_cpp(itensor_version)
        else: chain.setup_python()
    cm = np.zeros((n,n),dtype=complex)
    for i in range(n):
        for j in range(n):
            cm[i,j] = chain.vev(pairop(i,j), mode=mode)
    return cm

# --- spin chain: <Sz_i Sz_j> ---
n = 6
spins = ["S=1/2" for i in range(n)]
sc = spinchain.Spin_Chain(spins)
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
spin_pairop = lambda i,j: sc.Sz[i]*sc.Sz[j]

cm_ed = correlation_matrix(sc, spin_pairop, n, mode="ED")
for v in [2,3,"python"]:
    sc.set_hamiltonian(h) # re-set: switching backend invalidates the
                           # cached ground state (see restart(), called
                           # by set_hamiltonian) so each backend gets a
                           # fresh DMRG run
    cm_dmrg = correlation_matrix(sc, spin_pairop, n, mode="DMRG", itensor_version=v)
    diff = np.max(np.abs(cm_dmrg-cm_ed))
    print("Spin <Sz_i Sz_j>: max |DMRG(v%s) - ED| = %.2e"%(v,diff))
    assert diff<tol, "v%s spin correlator disagrees with ED by %g (tol=%g)"%(v,diff,tol)

# --- fermionic chain: <Cdag_i C_j> (Jordan-Wigner strings) ---
nf = 6
fc = fermionchain.Fermionic_Chain(nf)
hf = 0
for i in range(nf-1):
    hf = hf + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
fc.set_hamiltonian(hf)
fermion_pairop = lambda i,j: fc.Cdag[i]*fc.C[j]

cm_ed = correlation_matrix(fc, fermion_pairop, nf, mode="ED")
for v in [2,3,"python"]:
    fc.set_hamiltonian(hf)
    cm_dmrg = correlation_matrix(fc, fermion_pairop, nf, mode="DMRG", itensor_version=v)
    diff = np.max(np.abs(cm_dmrg-cm_ed))
    print("Fermion <Cdag_i C_j>: max |DMRG(v%s) - ED| = %.2e"%(v,diff))
    assert diff<tol, "v%s fermion correlator disagrees with ED by %g (tol=%g)"%(v,diff,tol)

print("TEST PASSED")
