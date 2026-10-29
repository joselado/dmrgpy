import sys, time; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import bosonchain
from dmrgpy.lindblad import Lindblad_Chain
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.mpsalgebra import randomMPS, _link_at
L,N=5,4
bc=bosonchain.Bosonic_Chain(L,maxnb=[N+1]*L,itensor_version="python"); bc.maxm=64
h=sum(-(bc.Adag[j+1]*bc.A[j]+bc.Adag[j]*bc.A[j+1]) for j in range(L-1))+0.5*sum(bc.N[j]*(bc.N[j]-1) for j in range(L)); bc.set_hamiltonian(h)
jumps=[(bc.Adag[j+1]+bc.Adag[j])*(bc.A[j+1]-bc.A[j]) for j in range(L-1)]
lc=Lindblad_Chain(bc,jumps); S=lc.doubled._session; G=S._mpo(lc.L.to_terms())
psi0=randomMPS(S.sites,64,rng=np.random.default_rng(0))
for one in (False,True):
    psi=psi0.copy(); psi.normalize()
    for m in range(3):
        t=time.time(); psi,_=clik._step(psi,G,0.5,1e-10,64,50,False,one); dt=time.time()-t
        print("one_site=%s step %d: %.1f s  bonds %s"%(one,m,dt,[_link_at(psi,i,i+1).dim for i in range(1,psi.length())]),flush=True)
