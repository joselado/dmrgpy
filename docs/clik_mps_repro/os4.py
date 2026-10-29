import sys; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import bosonchain
from dmrgpy.lindblad import Lindblad_Chain
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.tdvp import tdvp_step
from dmrgpy.pyitensor.gse import global_subspace_expand
from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor.mpsalgebra import inner, _link_at
from dmrgpy.pyitensor.metts import product_state
L,N=3,2
bc=bosonchain.Bosonic_Chain(L,maxnb=[N+1]*L,itensor_version="python"); bc.maxm=64
h=sum(-(bc.Adag[j+1]*bc.A[j]+bc.Adag[j]*bc.A[j+1]) for j in range(L-1))+0.5*sum(bc.N[j]*(bc.N[j]-1) for j in range(L)); bc.set_hamiltonian(h)
jumps=[(bc.Adag[j+1]+bc.Adag[j])*(bc.A[j+1]-bc.A[j]) for j in range(L-1)]
lc=Lindblad_Chain(bc,jumps); S=lc.doubled._session; G=S._mpo(lc.L.to_terms())
p,u=(product_state(S.sites,v) for v in lc._seeds([1,1,0],[0,1,1]))
ref=clik.evolve_branch(G,u,0.5,6,1e-10,64,one_site=False)
for gcut in (1e-10,1e-8,1e-6,1e-4):
    psi=u.copy(); psi.normalize(); out=[]
    for m in range(5):
        with bk.pad_bonds_suspended():
            psi=global_subspace_expand(G,psi,3,gcut,maxdim=64,bond_maxdim=64)
            psi=tdvp_step(psi.copy(),G,0.5j,cutoff=1e-10,maxdim=64,niter=50,num_center=1,hermitian=False)
        psi.normalize()
        out.append("%.1e/%d"%(1-abs(inner(ref[m+1],psi)),_link_at(psi,3,4).dim))
    print("gse cutoff %g: 1-overlap/middle bond per step:"%gcut," ".join(out),flush=True)
