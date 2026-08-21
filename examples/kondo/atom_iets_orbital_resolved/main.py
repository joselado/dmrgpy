# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

"""Orbital-resolved IETS spin-flip spectrum of a d^2 V(III) atom.

Two things are shown here, both of which used to be silently wrong:

1. `iorb` selects which orbital the STM tip couples to. The two entry
   points in atomtk/iets.py used to discard that argument and always
   compute orbital 0, so an in-plane (dxy) spectrum was returned when an
   out-of-plane (dz2) one was asked for. With a crystal field on, those
   are different spectra -- different weight *and* different peak
   position -- which is what the left panel plots.

2. `dex` in the ED dynamical correlator is a hard, equal-weight cutoff on
   what counts as the (degenerate) initial manifold. Sweeping a Zeeman
   field through a splitting comparable to `dex` makes the answer jump as
   levels cross the threshold. The right panel sweeps the field and
   overlays the equal-weight result against the Boltzmann-weighted
   finite-temperature one (pass a small `T` instead of tuning `dex`),
   which varies smoothly through the same sweep.
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import atom

orbs = ["dxy","dyz","dz2","dxz","dx2y2"] # d-shell, in this order
# crystal field: t2g (dxy,dyz,dxz) below eg (dz2,dx2y2), so the orbitals
# are genuinely inequivalent and `iorb` has something to select between
tij = np.diag([0.0,0.05,0.6,0.05,0.6])
es = np.linspace(0.,0.1,301) # bias window, in eV
delta = 1e-3 # broadening
dex = 1e-3 # equal-weight ground-state manifold cutoff


def spinflip(B=(0.,0.,0.),iorb=0,es=es,**kwargs):
    """Spin-flip IETS spectrum of the atom, for one orbital"""
    fc = atom.generate_atom(orbs=orbs,tij=tij,Ne=2,soc=0.03,U=4.,J=0.5,
                            B=list(B))
    with warnings.catch_warnings(): # the dex warning is the point of panel 2
        warnings.simplefilter("ignore")
        x,y = atom.get_spinflip(fc,es=es,iorb=iorb,delta=delta,**kwargs)
    return x,y


fig,(ax1,ax2) = plt.subplots(1,2,figsize=(11,4.2))

# --- panel 1: the spectrum really does depend on which orbital is probed
print("# orbital   spectral weight   peak (meV)")
for iorb in range(len(orbs)):
    x,y = spinflip(iorb=iorb,dex=dex)
    xp,yp = x[len(x)//2:],y[len(y)//2:] # positive-bias half
    weight = np.trapezoid(yp,xp) # integrated spectral weight
    peak = xp[np.argmax(yp)]*1e3 # peak position in meV
    print("%-9s %15.4f %13.2f"%(orbs[iorb],weight,peak))
    ax1.plot(x*1e3,y,label="%s (iorb=%d)"%(orbs[iorb],iorb))
ax1.set_xlabel("Bias [meV]") ; ax1.set_ylabel("IETS spin-flip signal")
ax1.set_title("Orbital-resolved spin flip\n(`iorb` is now respected)")
ax1.legend(fontsize=8)

# --- panel 2: dex is a step, a small T is not
Bs = np.linspace(0.,0.004,11) # Zeeman field sweep, splitting ~ through dex
esw = np.linspace(0.,0.1,81) # coarser grid: only the integral is used here
w_dex,w_T = [],[]
for B in Bs:
    x,y = spinflip(B=(0.,0.,B),iorb=0,es=esw,dex=dex) # equal-weight cutoff
    w_dex.append(np.trapezoid(y[len(y)//2:],x[len(x)//2:]))
    x,y = spinflip(B=(0.,0.,B),iorb=0,es=esw,T=1e-4) # Boltzmann-weighted instead
    w_T.append(np.trapezoid(y[len(y)//2:],x[len(x)//2:]))
ax2.plot(Bs*1e3,w_dex,"o-",label="equal weight, dex=%g"%dex)
ax2.plot(Bs*1e3,w_T,"s-",label="Boltzmann, T=1e-4")
ax2.set_xlabel("Zeeman field [meV]") ; ax2.set_ylabel("spectral weight")
ax2.set_title("Initial-manifold weighting\nacross a field sweep")
ax2.legend(fontsize=8)

plt.tight_layout()
plt.show()
