import sys, time; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src"); sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/tests")
import numpy as np
from dmrgpy import spinchain
import test_nh_dmrg as T
n=6
sc=spinchain.Spin_Chain([2]*n,itensor_version="python")
h=sum(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1] for i in range(n-1)); sc.set_hamiltonian(h)
for label,(c,hh) in (("heis",(sc,h)),("asym",T.nh_asymmetric_hopping_chain(6,"python"))):
    t=time.time(); ev,ch,_=c.get_krylov_spectrum(); t1=time.time()-t
    t=time.time(); ev2,rx,_=c.get_krylov_spectrum(exact_residuals=True); t2=time.time()-t
    E=c.get_krylov_spectrum(mode="ED")[0]
    err=[np.min(np.abs(E-x)) for x in ev[:6]]
    print(label,"t=%.1fs (exact %.1fs)"%(t1,t2))
    for k in range(6): print("   k=%d true err %.1e  change %.1e  exact res %.1e"%(k,err[k],ch[k],rx[k]))
