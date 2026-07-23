# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Regression test: ITensor v3 and the pure-Python backend must both honor
# a non-default, per-site Bosonic_Chain(maxnb=...) instead of silently
# truncating every site to the fixed 4-level BosonFourSite regardless of
# what was requested (a real, previously-silent DMRG/ED mismatch -- see
# mpscpp3/get_sites.h's 100+dim boson type-code range and
# extra/bosonfour.h's MaxOcc-driven site, mirrored on the pyitensor side
# by pyitensor/sites/boson.py's get_boson_site() factory). Only v3,
# python, and ED are compared here: v2 still only understands the single
# fixed-dim-4 boson type code (104), so it's out of scope for this test
# -- see examples/v2_VS_v3_boson for the default-dim=4 cross-backend
# comparison that still covers all four backends including v2.
import numpy as np
from dmrgpy import bosonchain

n = 3
maxnb = [3,5,3] # heterogeneous, non-default per-site dimension

def get_energy(itensor_version, mode="DMRG"):
    bc = bosonchain.Bosonic_Chain(n,maxnb=maxnb)
    if itensor_version!="python": bc.setup_cpp(itensor_version)
    else: bc.setup_python()
    bc.maxm = 120
    bc.nsweeps = 50
    np.random.seed(11)
    h = 0
    for i in range(n-1):
        h = h + np.random.random()*(bc.Adag[i]*bc.A[i+1] + bc.Adag[i+1]*bc.A[i])
    for i in range(n):
        h = h + 0.7*bc.N[i]*(bc.N[i]-1.0)
    bc.set_hamiltonian(h)
    return bc.gs_energy(mode=mode)

e3 = get_energy(3)
epy = get_energy("python")
eed = get_energy(3, mode="ED")

print("Ground state energy (ITensor v3,   maxnb=%s) ="%maxnb,e3)
print("Ground state energy (pure Python,  maxnb=%s) ="%maxnb,epy)
print("Ground state energy (ED,           maxnb=%s) ="%maxnb,eed)

tol = 1e-3
for name,e in [("python",epy),("ED",eed)]:
    diff = abs(e3-e)
    print("Difference v3 vs %s = %.2e"%(name,diff))
    assert diff<tol, "v3 vs %s disagree by %g (tol=%g)"%(name,diff,tol)

print("TEST PASSED")
