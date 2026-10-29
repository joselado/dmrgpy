import sys, time; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import bosonchain
from dmrgpy.lindblad import Lindblad_Chain
L, N = 3, 2
bc = bosonchain.Bosonic_Chain(L, maxnb=[N+1]*L, itensor_version="python")
h = sum(-(bc.Adag[j+1]*bc.A[j] + bc.Adag[j]*bc.A[j+1]) for j in range(L-1)) + 0.5*sum(bc.N[j]*(bc.N[j]-1) for j in range(L))
bc.set_hamiltonian(h); bc.maxm = 64
jumps = [(bc.Adag[j+1]+bc.Adag[j])*(bc.A[j+1]-bc.A[j]) for j in range(L-1)]
lc = Lindblad_Chain(bc, jumps)
for start in (None, ([1,1,0],[0,1,1])):
    t=time.time(); lam,res,obs = lc.get_spectrum(nk=20, dt=0.5, start=start)
    print("start",start,"%.0fs"%(time.time()-t)); print("  ",np.round(lam[:5],5)); print("  res",np.round(res[:5],6)); print("  n",np.round([obs(bc.N[j],0).real for j in range(L)],6))
bc.set_conserved_sector(Nb=N)
lamE,_,obsE = lc.get_spectrum(mode="ED"); print("ED",np.round(lamE[:5],5)); print("  n",np.round([obsE(bc.N[j],0).real for j in range(L)],6))
