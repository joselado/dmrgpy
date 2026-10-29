import sys; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import bosonchain
from dmrgpy.lindblad import Lindblad_Chain
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.mpsalgebra import inner, _link_at
from dmrgpy.pyitensor.metts import product_state
L,N=3,2
bc=bosonchain.Bosonic_Chain(L,maxnb=[N+1]*L,itensor_version="python"); bc.maxm=64
h=sum(-(bc.Adag[j+1]*bc.A[j]+bc.Adag[j]*bc.A[j+1]) for j in range(L-1))+0.5*sum(bc.N[j]*(bc.N[j]-1) for j in range(L)); bc.set_hamiltonian(h)
jumps=[(bc.Adag[j+1]+bc.Adag[j])*(bc.A[j+1]-bc.A[j]) for j in range(L-1)]
lc=Lindblad_Chain(bc,jumps); S=lc.doubled._session; G=S._mpo(lc.L.to_terms())
p,u=(product_state(S.sites,v) for v in lc._seeds([1,1,0],[0,1,1]))
for name,psi0 in (("p",p),("u",u)):
    a=clik.evolve_branch(G,psi0,0.5,6,1e-10,64,one_site=False)
    b=clik.evolve_branch(G,psi0,0.5,6,1e-10,64,one_site=True)
    for m in range(6):
        ov=abs(inner(a[m],b[m]))
        dims=[_link_at(b[m],i,i+1).dim for i in range(1,b[m].length())]
        dimsa=[_link_at(a[m],i,i+1).dim for i in range(1,a[m].length())]
        print(name,m,"|<two|one>| = %.8f"%ov, "bonds one-site",dims,"two-site",dimsa)
