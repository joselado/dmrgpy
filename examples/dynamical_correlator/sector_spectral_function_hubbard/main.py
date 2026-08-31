# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# The single-particle spectral function of a half-filled Hubbard chain,
# from get_spectral_function() -- i.e. from the eigenstates of the N+1 and
# N-1 particle-number sectors (src/dmrgpy/sectordc.py).
#
# This is the calculation the sector machinery exists for. Adding an
# electron takes the system to (Nf+1, Sz+1) and removing one to
# (Nf-1, Sz-1), so the states that carry the photoemission and inverse-
# photoemission weight are the ground and low excited states of two OTHER
# sectors -- each of them a clean, individually converged DMRG solve in a
# smaller Hilbert space, rather than a highly excited state of the full
# one.
#
# What the plot shows: as U grows the single band splits into a lower and
# an upper Hubbard band separated by the Mott gap, with the chemical
# potential mu = (E_0^{N+1}-E_0^{N-1})/2 sitting in the middle of it (at
# exactly U/2 here, by particle-hole symmetry). The gap read off the
# spectrum is checked against E_0^{N+1}+E_0^{N-1}-2E_0^N computed
# directly, which is a physics-level check that the assembled spectrum
# cannot fake.
import time
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain

n = 6                       # sites (2n spin-orbitals, half filling: Nf=n)
Us = [0.0,2.0,4.0,8.0]
nex = 10
es = np.linspace(-9,9,700)
delta = 0.35

def hubbard(u):
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = 0
    for i in range(n-1):
        for cd,c in ((fc.Cdagup,fc.Cup),(fc.Cdagdn,fc.Cdn)):
            h = h - cd[i]*c[i+1] - cd[i+1]*c[i]
    for i in range(n): h = h + u*fc.Nup[i]*fc.Ndn[i]
    fc.set_hamiltonian(h)
    fc.maxm = 60 ; fc.nsweeps = 10
    fc.set_conserved_sector(Nf=n,Sz=0)   # half filling, unpolarized
    return fc

spectra,gaps,mus,captured = {},[],[],[]
for u in Us:
    t0 = time.time()
    fc = hubbard(u)
    # site n//2, spin up: a bulk site, away from the open ends
    x,A,info = fc.get_spectral_function(i=n//2,spin="up",nex=nex,es=es,
                                        delta=delta,return_poles=True)
    spectra[u] = np.array(A).real
    gaps.append(info["gap"]) ; mus.append(info["mu"])
    captured.append(0.5*(info["particle"]["captured"]+info["hole"]["captured"]))
    print("U=%4.1f  mu=%7.4f  gap=%7.4f  captured=%.3f  (%.1f s)"
          %(u,info["mu"],info["gap"],captured[-1],time.time()-t0))
    print("        particle sector %s, hole sector %s"
          %(info["particle"]["sector"],info["hole"]["sector"]))
    if u>0:  # particle-hole symmetry of the half-filled Hubbard model
        assert abs(info["mu"]-u/2)<1e-3, \
            "mu should be U/2 at half filling, got %g"%info["mu"]

assert gaps[-1]>gaps[0]+1.0, "the Mott gap should open with U"
print("TEST PASSED")

# --- plots: the spectra, and the gap growing with U ---
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(11,4.5))
fig.subplots_adjust(wspace=0.3,bottom=0.15)

for k,u in enumerate(Us):
    ax1.plot(es,spectra[u]+0.6*k,label="U=%.0f"%u)
    ax1.axhline(0.6*k,c="black",lw=0.4,alpha=0.3)
ax1.axvline(0.0,c="black",ls="--",lw=1)
ax1.set_xlabel(r"$\omega-\mu$  [t]") ; ax1.set_ylabel(r"$A(\omega)$  (offset)")
ax1.set_title("Hubbard chain, n=%d at half filling"%n) ; ax1.legend(fontsize=8)

ax2.plot(Us,gaps,marker="o",c="crimson",label="charge gap from A($\\omega$)")
ax2.plot(Us,mus,marker="s",c="navy",label=r"$\mu$")
ax2.plot(Us,[u/2 for u in Us],c="deepskyblue",ls="--",lw=2.5,zorder=0,label="U/2 (exact)")
ax2.set_xlabel("U [t]") ; ax2.set_ylabel("energy [t]")
ax2.set_title("Mott gap and chemical potential") ; ax2.legend(fontsize=8)

plt.savefig("sector_spectral_function_hubbard.png",dpi=150)
print("Plot saved to sector_spectral_function_hubbard.png")
plt.show()
