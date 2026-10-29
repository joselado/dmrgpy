import sys; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import spinchain, timedependent as td
def run(mode,GAMMA):
    sc = spinchain.Spin_Chain([2]*4,itensor_version="python")
    h0 = 0.3*sc.Sz[0] - 0.2*sc.Sz[3] + 0.4*sc.Sx[1] + 0.5*sc.Sz[2] + 0.1*sc.Sx[2]
    sc.set_hamiltonian(h0); sc.maxm = 20
    wf = sc.get_gs(mode=mode)
    h = 0
    for i in range(3): h = h + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.5*sc.Sx[0] + GAMMA*sc.Sz[0]   # non-Hermitian onsite term
    out = td.evolve_and_measure(sc,operator=sc.Sz[2],nt=40,dt=0.05,h=h,wf=wf,mode=mode)
    return out
import sys as S; G=complex(S.argv[1]); a = run("DMRG",G); b = run("ED",G)
print(type(a), [np.shape(x) for x in a]); 
a=np.array(a[1]); b=np.array(b[1])
print(np.round(a[::8],4)); print(np.round(b[::8],4)); print("max diff", np.max(np.abs(a-b)))
