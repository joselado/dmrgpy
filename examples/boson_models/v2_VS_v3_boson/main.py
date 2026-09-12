# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

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

def get_energy(itensor_version, U, mode="DMRG"):
    bc = bosonchain.Bosonic_Chain(n)
    # Switch backend *before* seeding random couplings: constructing the
    # pure-Python backend consumes draws from numpy's global RNG as a side
    # effect (unlike v2/v3, which never touch np.random), so seeding
    # first would silently desync "the same" Hamiltonian across backends
    # -- see examples/v2_VS_v3_parafermion/main.py for the concrete
    # failure this caused when gotten wrong.
    if itensor_version!="python": bc.setup_cpp(itensor_version)
    else: bc.setup_python()
    # Pin the sweep schedule instead of relying on the library defaults
    # (nsweeps=15, maxm=30). Both compiled backends start DMRG from an
    # unseeded random MPS, and at the defaults their convergence tail
    # occasionally left one U point ~1e-2 above the exact answer -- which
    # made this assert-carrying script fail on a clean tree roughly 40% of
    # the time, symmetrically for v2 and v3 (the 2026-09 audit, finding
    # #24). At n=6 the exact MPS bond dimension is only 4**3=64, so
    # maxm=100 truncates nothing at all and the whole sweep schedule costs
    # about 25s for the whole script. Measured over 16 full runs of this
    # script, the largest disagreement between any two backends at any U
    # is 1.7e-7 at these settings, against 1e-3 at nsweeps=40 and 2e-2 at
    # the defaults -- which is what lets the assert below run at tol=1e-5
    # (60x above the worst observed) instead of the 1e-2 that was still
    # not loose enough to be deterministic. Do not lower these.
    bc.nsweeps = 80 # enough sweeps to converge from a random start
    bc.maxm = 100 # above 4**3=64, i.e. no truncation at all at n=6
    np.random.seed(11)
    h = 0
    for i in range(n-1):
        h = h + np.random.random()*(bc.Adag[i]*bc.A[i+1] + bc.Adag[i+1]*bc.A[i])
    for i in range(n):
        den = bc.Adag[i]*bc.A[i]
        h = h + U*den*den
    bc.set_hamiltonian(h)
    return bc.gs_energy(mode=mode)

U0 = 0.3 # onsite interaction strength used in the original example
e2 = get_energy(2,U0)
e3 = get_energy(3,U0)
eed = get_energy(2,U0,mode="ED")
epy = get_energy("python",U0) # n=6 is small enough for the python backend

print("Ground state energy (ITensor v2)  =",e2)
print("Ground state energy (ITensor v3)  =",e3)
print("Ground state energy (ED)          =",eed)
print("Ground state energy (pure Python) =",epy)

tol = 1e-5 # supported by the pinned nsweeps/maxm in get_energy()
for name,e in [("v3",e3),("ED",eed),("python",epy)]:
    diff = abs(e2-e)
    print("Difference v2 vs %s = %.2e"%(name,diff))
    assert diff<tol, "v2 vs %s disagree by %g (tol=%g)"%(name,diff,tol)

print("TEST PASSED")

# sweep the onsite interaction strength U -- a cheap parameter already
# present in get_energy() -- and check all backends keep agreeing across
# it, not just at the single U=0.3 point above
Us = [0.1,0.2,0.3,0.4,0.5]
e2s,e3s,eeds,epys = [],[],[],[]
for U in Us:
    e2U = get_energy(2,U)
    e3U = get_energy(3,U)
    eedU = get_energy(2,U,mode="ED")
    epyU = get_energy("python",U)
    print("U =",U,"  v2 =",e2U,"  v3 =",e3U,"  ED =",eedU,"  python =",epyU)
    for name,e in [("v3",e3U),("ED",eedU),("python",epyU)]:
        diff = abs(e2U-e)
        assert diff<tol, "U=%g: v2 vs %s disagree by %g (tol=%g)"%(U,name,diff,tol)
    e2s.append(e2U); e3s.append(e3U); eeds.append(eedU); epys.append(epyU)

import matplotlib.pyplot as plt
plt.plot(Us,e2s,marker="o",label="ITensor v2")
plt.plot(Us,e3s,marker="s",linestyle="--",label="ITensor v3")
plt.plot(Us,eeds,marker="^",linestyle=":",label="ED")
plt.plot(Us,epys,marker="d",linestyle="-.",label="pure Python")
plt.xlabel("Onsite interaction U")
plt.ylabel("Ground state energy")
plt.legend()
plt.show()
