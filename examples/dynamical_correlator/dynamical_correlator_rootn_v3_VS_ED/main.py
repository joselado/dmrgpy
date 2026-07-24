# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

# Validate the root-N Krylov-space correction-vector method
# (Nocera & Alvarez, arXiv:2204.03165) on the ITensor v3 DMRG backend
# (mode="DMRG", itensor_version=3) against exact ED (mode="ED"), the same
# check examples/dynamical_correlator/dynamical_correlator_rootn_ED
# already does for the ED-only prototype (src/dmrgpy/algebra/rootn.py).
#
# The v3 (and v2/"python") implementation lives in rootndmrg.py: it
# builds the correction vector via the same N-step fractional-power
# Lanczos recursion as algebra/rootn.py, but with the Krylov subspace
# built out of *global* MPS vectors (repeated truncated MPO application)
# instead of numpy vectors -- see rootndmrg.py's module docstring for why
# a per-bond local sweep (the more "natural" DMRG-style approach) turned
# out to give wrong results and was abandoned in favor of this simpler,
# provably-correct design.

n = 6
spins = ["S=1/2" for i in range(n)]
sc = spinchain.Spin_Chain(spins)
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

name = (sc.Sz[0],sc.Sz[0])
es = np.linspace(0.5,3.0,10)
delta = 2e-1

# --- Check 1: self-consistency at a Krylov dimension that spans the
# whole (2^6=64-dimensional) Hilbert space -- root-N (any N) must then
# reproduce the exact ED Lehmann-sum answer, since the Krylov subspace at
# every step is exact. This isolates the recursion's own correctness
# from the Krylov-truncation approximation, on the actual v3 MPS/MPO
# machinery (not just the ED prototype).
sc.itensor_version = 3
sc.setup_cpp(version=3)
sc.maxm = 60 # comfortably above what a 6-site chain ever needs
sc.cutoff = 1e-10

x_ed,y_ed = sc.get_dynamical_correlator(mode="ED",submode="ED",
        name=name,es=es,delta=delta)
x_v3,y_v3 = sc.get_dynamical_correlator(mode="DMRG",submode="ROOTN",
        name=name,es=es,delta=delta,N=6,nkry=25)

err_full = np.max(np.abs(y_v3.real-y_ed.real))
print("Self-consistency check (v3, nkry=25 out of dim=64, N=6):")
print("  max|v3 ROOTN - exact ED| =",err_full)
assert err_full<1e-8, "v3 ROOTN must reproduce exact ED when the Krylov subspace is effectively exhaustive"

# --- Check 2: accuracy at a genuinely restricted Krylov dimension, near
# the top of the spectral bandwidth -- same qualitative check as the
# ED-only example: does root-N (N>1) improve over the conventional
# single-shot Krylov correction vector (N=1) at the same fixed budget?
es_hi = np.linspace(2.0,3.5,15)
x_ed_hi,y_ed_hi = sc.get_dynamical_correlator(mode="ED",submode="ED",
        name=name,es=es_hi,delta=1.5e-1)

nkry = 4
results = {}
for N in [1,4,8]:
    x,y = sc.get_dynamical_correlator(mode="DMRG",submode="ROOTN",
            name=name,es=es_hi,delta=1.5e-1,N=N,nkry=nkry)
    err = np.abs(y.real-y_ed_hi.real)
    results[N] = (x,y,err)
    print("  N=%d  nkry=%d  mean|err|=%.3e  max|err|=%.3e"%(
        N,nkry,np.mean(err),np.max(err)))

assert np.mean(results[4][2])<np.mean(results[1][2]), \
        "root-N (N=4) should reduce mean error over conventional CV (N=1) at fixed Krylov budget"

import matplotlib.pyplot as plt
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
fig.subplots_adjust(wspace=0.3,bottom=0.15)

ax1.plot(x_ed,y_ed.real,c="black",label="exact ED",lw=2)
ax1.plot(x_v3,y_v3.real,c="red",label="v3 ROOTN (N=6,nkry=25)",
        marker="o",ms=4,ls="none")
ax1.axvspan(es_hi.min(),es_hi.max(),color="gray",alpha=0.2,label="zoom (right panel)")
ax1.legend(fontsize=8)
ax1.set_xlabel("frequency [J]")
ax1.set_ylabel("Dynamical correlator")
ax1.set_title("v3 ROOTN vs exact ED (n=%d sites)"%n)

ax2.plot(x_ed_hi,y_ed_hi.real,c="black",label="exact ED",lw=2)
colors = {1:"red",4:"orange",8:"green"}
for N in [1,4,8]:
    x,y,_ = results[N]
    ax2.plot(x,y.real,c=colors[N],label="v3 ROOTN N=%d, nkry=%d"%(N,nkry),marker="o",ms=3)
ax2.legend(fontsize=8)
ax2.set_xlabel("frequency [J]")
ax2.set_title("Upper bandwidth, restricted Krylov dim")

plt.savefig("rootn_v3_vs_ed.png",dpi=150)
print("\nPlot saved to rootn_v3_vs_ed.png")
plt.show()
