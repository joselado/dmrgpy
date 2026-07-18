# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Regression test: ground state energy of a bosonic Bose-Hubbard-like
# chain must agree across ITensor v2, v3, ED, and the pure-Python
# backend. There was previously no v2_VS_v3_* comparison for the boson
# model at all (unlike spin, fermion, parafermion chains) -- this fills
# that gap, at a size (n=6) beyond the n=3 covered by
# examples/bosonic_hubbard/main.py. n is capped at 6 (not pushed further,
# unlike the other v2_VS_v3_* tests) because ED needs the full Hilbert
# space dimension (Bosonic_Chain's default per-site dimension^n) under
# get_ED_obj()'s own 10000 cutoff -- 4**6=4096, 4**7=16384.
import numpy as np
from dmrgpy import bosonchain

n = 6

def get_energy(itensor_version, mode="DMRG"):
    bc = bosonchain.Bosonic_Chain(n)
    # Switch backend *before* seeding random couplings: constructing the
    # pure-Python backend consumes draws from numpy's global RNG as a side
    # effect (unlike v2/v3, which never touch np.random), so seeding
    # first would silently desync "the same" Hamiltonian across backends
    # -- see examples/v2_VS_v3_parafermion/main.py for the concrete
    # failure this caused when gotten wrong.
    if itensor_version!="python": bc.setup_cpp(itensor_version)
    else: bc.setup_python()
    np.random.seed(11)
    h = 0
    for i in range(n-1):
        h = h + np.random.random()*(bc.Adag[i]*bc.A[i+1] + bc.Adag[i+1]*bc.A[i])
    for i in range(n):
        den = bc.Adag[i]*bc.A[i]
        h = h + 0.3*den*den
    bc.set_hamiltonian(h)
    return bc.gs_energy(mode=mode)

e2 = get_energy(2)
e3 = get_energy(3)
eed = get_energy(2, mode="ED")
epy = get_energy("python") # n=8 < 10, python backend is in scope

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
