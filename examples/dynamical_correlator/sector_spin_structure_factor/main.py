# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# The dynamical spin structure factor S(q,omega) of a Heisenberg chain,
# built from the eigenstates of the Sz=+1 sector (submode="SECTOR", see
# src/dmrgpy/sectordc.py).
#
# The selection rule is what makes this cheap. S+ raises Sz by one, so
# every state that carries transverse spectral weight lives in the Sz=+1
# sector and nowhere else -- and within it, the matrix elements for ALL
# sites come from the same set of states:
#
#     S(q,w) = sum_n |m_n(q)|^2 L(w - (E_n - E_0)),
#     m_n(q) = sqrt(2/(L+1)) sum_i sin(q(i+1)) <n|S^+_i|0>
#
# (the sine transform is the right one for an open chain, whose normal
# modes are standing waves with q = pi k/(L+1)). So one set of sector
# solves yields the whole (q,w) map: sectordc.py caches those solves
# across calls, and only the L matrix-element vectors are recomputed.
#
# The physics to look for: spectral weight concentrated just above the
# des Cloizeaux-Pearson lower boundary w = (pi/2) J |sin q|, which is the
# two-spinon continuum's lower edge, drawn on the map for comparison. The
# weight piles up at q=pi, the antiferromagnetic wavevector.
import time
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain, sectordc
from dmrgpy.operatornames import name2MO

L = 12
nex = 30
es = np.linspace(0.0,4.0,400)
delta = 0.12

sc = spinchain.Spin_Chain([2 for i in range(L)])
h = 0
for i in range(L-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 60 ; sc.nsweeps = 10
Sp,Sm = name2MO("Sp",sc),name2MO("Sm",sc)

# One matrix-element vector per site: m[n,i] = <n|S^+_i|0>. Only the first
# call solves the two sectors; the rest hit sectordc.py's cache.
t0 = time.time()
mel,energies = [],None
for i in range(L):
    e,w,info = sectordc.sector_poles(sc,name=(Sm[0],Sp[i]),nex=nex,quiet=True)
    if energies is None: energies = e
    mel.append(info["matrix_elements"])
    print("site %2d done (%5.1f s), captured %.4f"%(i,time.time()-t0,info["captured"]))
mel = np.array(mel).T                     # (state, site)
print("target sector",info["target_sector"],"with",len(energies),"states")

# sine transform to the open chain's normal modes
ks = np.arange(1,L+1)
qs = np.pi*ks/(L+1)
sine = np.sqrt(2.0/(L+1))*np.sin(np.outer(qs,np.arange(1,L+1)))   # (q, site)
mq = mel@sine.T                                                   # (state, q)

S = np.zeros((len(qs),len(es)))
for iq in range(len(qs)):
    for n in range(len(energies)):
        S[iq] += abs(mq[n,iq])**2*delta/np.pi/((es-energies[n])**2+delta**2)

# the weight must sit above the des Cloizeaux-Pearson lower edge
dcp = 0.5*np.pi*np.abs(np.sin(qs))
peak = es[np.argmax(S,axis=1)]
print("q/pi      ",np.round(qs/np.pi,3))
print("peak omega",np.round(peak,3))
print("dCP edge  ",np.round(dcp,3))
assert np.all(peak > 0.75*dcp-3*delta), \
    "spectral weight found well below the two-spinon lower edge"
assert np.argmax(S.max(axis=1))>=L-3, \
    "the strongest response should sit near q=pi for an antiferromagnet"
print("TEST PASSED")

# --- plots: the S(q,w) map, and cuts at a few wavevectors ---
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,4.5))
fig.subplots_adjust(wspace=0.3,bottom=0.15)

im = ax1.pcolormesh(qs/np.pi,es,S.T,shading="auto",cmap="magma")
qq = np.linspace(0.02,np.pi,200)
ax1.plot(qq/np.pi,0.5*np.pi*np.sin(qq),c="cyan",ls="--",lw=1.5,
         label=r"$\frac{\pi}{2}J|\sin q|$ (des Cloizeaux-Pearson)")
ax1.plot(qq/np.pi,np.pi*np.abs(np.sin(qq/2)),c="white",ls=":",lw=1.2,
         label=r"$\pi J|\sin(q/2)|$ (upper edge)")
ax1.set_xlabel(r"$q/\pi$") ; ax1.set_ylabel(r"$\omega$ [J]")
ax1.set_title(r"$S^{+-}(q,\omega)$, L=%d Heisenberg chain"%L)
ax1.legend(fontsize=7,loc="upper left") ; fig.colorbar(im,ax=ax1)

for iq in [len(qs)//4,len(qs)//2,len(qs)-1]:
    ax2.plot(es,S[iq],label=r"$q=%.2f\pi$"%(qs[iq]/np.pi))
ax2.set_xlabel(r"$\omega$ [J]") ; ax2.set_ylabel(r"$S^{+-}(q,\omega)$")
ax2.set_title("cuts at fixed wavevector") ; ax2.legend(fontsize=8)

plt.savefig("sector_spin_structure_factor.png",dpi=150)
print("Plot saved to sector_spin_structure_factor.png")
plt.show()
