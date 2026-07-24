# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

# Validate the root-N Krylov-space correction-vector method
# (Nocera & Alvarez, arXiv:2204.03165) against exact ED, on a small
# Heisenberg chain where the exact Lehmann-sum answer (submode="ED") is
# cheap to obtain and can serve as ground truth.
#
# Two checks:
#  1) Self-consistency: with a Krylov dimension nkry equal to (or larger
#     than) the Hilbert space size, the Lanczos subspace built at each
#     root-N step is exact, so root-N (any N) must reproduce submode="ED"
#     exactly, up to floating point. This isolates bugs in the root-N
#     recursion itself from the Krylov-truncation approximation.
#  2) The actual point of the method: at a small, fixed Krylov dimension
#     nkry (analogous to a bond-dimension cap in the MPS/DMRG setting),
#     compare the conventional Krylov-space correction vector (N=1, one
#     shot) against root-N with N=4,8,20 at the same nkry, and see
#     whether splitting the same "budget" into more, smaller steps
#     improves accuracy against the exact answer -- as reported in the
#     paper for MPS bond dimension.

n = 6
spins = ["S=1/2" for i in range(n)]
sc = spinchain.Spin_Chain(spins)
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

i,j = 0,0
name = (sc.Sz[i],sc.Sz[j])
es_full = np.linspace(-0.5,5.0,80)
delta = 2e-1

# Hilbert space size for this chain (spin-1/2, n sites) = 2**n
dim = 2**n

x_exact,y_exact = sc.get_dynamical_correlator(mode="ED",submode="ED",
        name=name,es=es_full,delta=delta)

# --- Check 1: self-consistency at full Krylov dimension ---
x_full,y_full = sc.get_dynamical_correlator(mode="ED",submode="ROOTN",
        name=name,es=es_full,delta=delta,N=6,nkry=dim)
err_full = np.max(np.abs(y_full-y_exact))
print("Self-consistency check (nkry=full Hilbert space, N=6):")
print("  max|rootN - exact| =",err_full)
assert err_full<1e-6, "root-N must reproduce the exact answer when nkry spans the full Hilbert space"

# --- Check 2: accuracy at fixed, restricted Krylov dimension, near the TOP
# of the spectral bandwidth -- this is exactly where the paper reports
# root-N's advantage over the conventional (N=1) Krylov correction vector,
# since that is where a bond-dimension- (here: Krylov-dimension-) limited
# single-shot correction vector struggles most.
emax = np.max(sc.get_ED_obj().get_diagonalized_hamiltonian()[0])
e0 = sc.gs_energy(mode="ED")
es_hi = np.linspace(0.55*(emax-e0),0.99*(emax-e0),30) # upper part of the bandwidth
x_exact_hi,y_exact_hi = sc.get_dynamical_correlator(mode="ED",submode="ED",
        name=name,es=es_hi,delta=delta)

nkry = 4 # deliberately small, well below dim, to expose truncation error
x_n1,y_n1 = sc.get_dynamical_correlator(mode="ED",submode="ROOTN",
        name=name,es=es_hi,delta=delta,N=1,nkry=nkry) # conventional Krylov CV
x_n4,y_n4 = sc.get_dynamical_correlator(mode="ED",submode="ROOTN",
        name=name,es=es_hi,delta=delta,N=4,nkry=nkry)
x_n8,y_n8 = sc.get_dynamical_correlator(mode="ED",submode="ROOTN",
        name=name,es=es_hi,delta=delta,N=8,nkry=nkry)

err_n1 = np.abs(y_n1-y_exact_hi)
err_n4 = np.abs(y_n4-y_exact_hi)
err_n8 = np.abs(y_n8-y_exact_hi)
print("\nAccuracy near the top of the bandwidth, fixed Krylov dimension nkry=%d (Hilbert space dim=%d):"%(nkry,dim))
print("  N=1  (conventional CV) mean|err| =",np.mean(err_n1),"max|err| =",np.max(err_n1))
print("  N=4  (root-N)          mean|err| =",np.mean(err_n4),"max|err| =",np.max(err_n4))
print("  N=8  (root-N)          mean|err| =",np.mean(err_n8),"max|err| =",np.max(err_n8))

import matplotlib.pyplot as plt
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
fig.subplots_adjust(wspace=0.3,bottom=0.15)

ax1.plot(x_exact,y_exact.real,c="black",label="exact (Lehmann sum)",lw=2)
ax1.axvspan(es_hi.min(),es_hi.max(),color="gray",alpha=0.2,label="zoom (right panel)")
ax1.legend(fontsize=8)
ax1.set_xlabel("frequency [J]")
ax1.set_ylabel("Dynamical correlator")
ax1.set_title("Full spectrum (n=%d sites)"%n)

ax2.plot(x_exact_hi,y_exact_hi.real,c="black",label="exact (Lehmann sum)",lw=2)
ax2.plot(x_n1,y_n1.real,c="red",label="conventional CV (N=1), nkry=%d"%nkry,marker="o",ms=3)
ax2.plot(x_n4,y_n4.real,c="orange",label="root-N, N=4, nkry=%d"%nkry,marker="o",ms=3)
ax2.plot(x_n8,y_n8.real,c="green",label="root-N, N=8, nkry=%d"%nkry,marker="o",ms=3)
ax2.legend(fontsize=8)
ax2.set_xlabel("frequency [J]")
ax2.set_title("Upper bandwidth, restricted Krylov dim")

plt.savefig("rootn_vs_ed.png",dpi=150)
print("\nPlot saved to rootn_vs_ed.png")
plt.show()
