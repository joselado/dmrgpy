import sys; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import bosonchain
from dmrgpy.lindblad import Lindblad_Chain
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.tdvp import tdvp_step
from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor.mpsalgebra import inner, _link_at
from dmrgpy.pyitensor.metts import product_state
L,N=3,2
bc=bosonchain.Bosonic_Chain(L,maxnb=[N+1]*L,itensor_version="python"); bc.maxm=64
h=sum(-(bc.Adag[j+1]*bc.A[j]+bc.Adag[j]*bc.A[j+1]) for j in range(L-1))+0.5*sum(bc.N[j]*(bc.N[j]-1) for j in range(L)); bc.set_hamiltonian(h)
jumps=[(bc.Adag[j+1]+bc.Adag[j])*(bc.A[j+1]-bc.A[j]) for j in range(L-1)]
lc=Lindblad_Chain(bc,jumps); S=lc.doubled._session; G=S._mpo(lc.L.to_terms())
p,u=(product_state(S.sites,v) for v in lc._seeds([1,1,0],[0,1,1]))
a=clik.evolve_branch(G,u,0.5,4,1e-10,64,one_site=False)
x=a[2]; print("start bonds",[_link_at(x,i,i+1).dim for i in range(1,x.length())])
def nrm(s): s=s.copy(); s.normalize(); return s
for dt in (0.5,0.1,0.02):
    two=nrm(tdvp_step(x.copy(),G,1j*dt,cutoff=1e-10,maxdim=64,niter=50,hermitian=False))
    with bk.pad_bonds_suspended():
        one=nrm(tdvp_step(x.copy(),G,1j*dt,cutoff=1e-10,maxdim=64,niter=50,num_center=1,hermitian=False))
    print("dt=%g one-site vs two-site, no expansion: 1-|<.|.>| = %.2e  bonds one %s"%(dt,1-abs(inner(two,one)),[_link_at(one,i,i+1).dim for i in range(1,one.length())]))
