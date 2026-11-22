"""DMRG of the Anderson impurity in the site basis and in the Gaussian
basis, same Hamiltonian, same v3 backend, same sweeps, cold start."""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'src'))
import numpy as np, time, json
from dmrgpy import fermionchain, wilsonchain

def uniform(n, V=1.0, t=1.0):
    h = np.zeros((n,n)); h[0,1]=h[1,0]=-V
    for i in range(1,n-1): h[i,i+1]=h[i+1,i]=-t
    return h
def wilson(n, tD=0.2, Lambda=2.0):
    eps, t = wilsonchain.wilson_chain(n-1, Lambda=Lambda)
    h = np.diag(np.concatenate([[0.0], eps])); h[0,1]=h[1,0]=tD
    for k in range(n-2): h[k+1,k+2]=h[k+2,k+1]=t[k]
    return h
def hamiltonian(fc, h, U):
    n = h.shape[0]
    H = sum(h[i,j]*(fc.Cdagup[i]*fc.Cup[j]+fc.Cdagdn[i]*fc.Cdn[j]) for i in range(n) for j in range(n) if h[i,j]!=0)
    return H + U*(fc.Nup[0]-0.5)*(fc.Ndn[0]-0.5)
def dims(fc, wf):
    s = fc._session
    mpo = max(np.asarray(a).shape[-1] for a in s.mpo_site_arrays(s.hamiltonian_mpo()))
    mps = max(np.asarray(a).shape[-1] for a in s.mps_site_arrays(wf.cpp_handle))
    return int(mpo), int(mps)

case, n, U = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
maxms = [int(x) for x in sys.argv[4].split(",")]
h = uniform(n) if case=="uniform" else wilson(n)
out = dict(case=case, n=n, U=U, runs=[])
for basis in ["site", "gaussian"]:
    for m in maxms:
        fc = fermionchain.Spinful_Fermionic_Chain(n)
        fc.maxm = m; fc.nsweeps = 12
        H = hamiltonian(fc, h, U)
        t0 = time.time()
        if basis == "gaussian":
            b = fc.get_gaussian_basis(h=h, nf=n, fixed=[0])
            H = b.rotate(H)
        t1 = time.time()
        fc.set_hamiltonian(H)
        e = fc.gs_energy(); t2 = time.time()
        wf = fc.get_gs()
        S = [float(wf.get_bond_entropy(i)) for i in range(fc.ns-1)]
        mpo, mps = dims(fc, wf)
        d = float(fc.vev(fc.Nup[0]*fc.Ndn[0]).real)
        r = dict(basis=basis, maxm=m, E=float(e), docc=d, S=S, mpo=mpo, mps=mps, t_basis=t1-t0, t=t2-t1)
        out["runs"].append(r)
        print("%-8s n=%d maxm=%4d E=%.10f docc=%.6f Smax=%.3f mpo=%d mps=%d t_basis=%.1fs t=%.1fs"%(basis,n,m,e,d,max(S),mpo,mps,t1-t0,t2-t1), flush=True)
        json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'cmp_%s_n%d_U%g.json'%(case,n,U)),'w'))
