import sys; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import spinchain, bosonchain
from dmrgpy.lindblad import Lindblad_Chain
def compare(sc,jumps,label):
    lc=Lindblad_Chain(sc,jumps)
    lam,_,_=lc.get_spectrum(mode="ED")
    lv,_,_=lc.doubled.get_krylov_spectrum(G=lc.L,mode="ED")
    # match spectra as multisets
    d=max(np.min(np.abs(lv-x)) for x in lam); d2=max(np.min(np.abs(lam-x)) for x in lv)
    print(label,"dim",len(lam),len(lv)," max mismatch",d,d2)
n=3
sc=spinchain.Spin_Chain([2]*n,itensor_version="python")
H=sum(sc.Sx[i]*sc.Sy[i+1]-0.3*sc.Sy[i]*sc.Sz[i+1]+sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))+0.4*sc.Sy[1]; sc.set_hamiltonian(H)
Sp=lambda i: sc.Sx[i]+1j*sc.Sy[i]
compare(sc,[0.5*Sp(0)*sc.Sz[1], (0.3+0.2j)*sc.Sy[2], 0.7*(sc.Sx[1]-1j*sc.Sy[1])],"spin, Sy/S+-/products")
sc1=spinchain.Spin_Chain([3,3],itensor_version="python")
H1=sc1.Sx[0]*sc1.Sy[1]+0.3*sc1.Sz[0]*sc1.Sz[0]; sc1.set_hamiltonian(H1)
compare(sc1,[0.5*(sc1.Sx[0]+1j*sc1.Sy[0]), 0.4*sc1.Sy[1]*sc1.Sz[1]],"spin-1, same-site products")
bc=bosonchain.Bosonic_Chain(3,maxnb=[3,3,3],itensor_version="python")
Hb=sum(-1.0*(bc.Adag[i]*bc.A[i+1]+bc.Adag[i+1]*bc.A[i]) for i in range(2))+0.5*sum(bc.N[i]*(bc.N[i]-1) for i in range(3)); bc.set_hamiltonian(Hb)
compare(bc,[(bc.Adag[1]+bc.Adag[0])*(bc.A[1]-bc.A[0]), 0.3*bc.A[2], 0.2j*bc.N[0]],"bosons, dissipative hopping")
