import sys, time; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import bosonchain
from dmrgpy.lindblad import Lindblad_Chain
L, N, J, U, kappa = 5, 4, 1.0, 1.0, 1.0
maxm = int(sys.argv[1]); nk=int(sys.argv[2]); dt=float(sys.argv[3])
bc = bosonchain.Bosonic_Chain(L, maxnb=[N+1]*L, itensor_version="python")
h = sum(-J*(bc.Adag[j+1]*bc.A[j] + bc.Adag[j]*bc.A[j+1]) for j in range(L-1))
h = h + 0.5*U*sum(bc.N[j]*(bc.N[j]-1) for j in range(L))
bc.set_hamiltonian(h)
jumps = [np.sqrt(kappa)*(bc.Adag[j+1]+bc.Adag[j])*(bc.A[j+1]-bc.A[j]) for j in range(L-1)]
bc.maxm = maxm
lc = Lindblad_Chain(bc, jumps)
t=time.time()
lam, res, obs = lc.get_spectrum(nk=nk, dt=dt, start=([1,1,1,1,0],[0,1,1,1,1]))
print("DMRG maxm=%d nk=%d dt=%g: %.0f s"%(maxm,nk,dt,time.time()-t), flush=True)
for k in range(6): print("  ", np.round(lam[k],5), "res %.1e"%res[k])
print("  <n_j>_ss", np.round([obs(bc.N[j],0).real for j in range(L)],5), flush=True)
if len(sys.argv)>4:
    bc.set_conserved_sector(Nb=N)
    t=time.time(); lamE,resE,obsE = lc.get_spectrum(mode="ED"); print("ED %.0f s"%(time.time()-t))
    for k in range(6): print("  ", np.round(lamE[k],5))
    print("  <n_j>_ss", np.round([obsE(bc.N[j],0).real for j in range(L)],5))
