from common import *
from dmrgpy.pyitensor import kernels, dmrg
S, G, psi = lindblad_case()
tdvp.tdvp_step(psi, G, 1j*0.5, cutoff=1e-10, maxdim=64, niter=50, hermitian=False)
cap=[]
orig = dmrg.two_site_heff
def spy(L, Lbra, H, ket, i, R, Rbra):
    out = orig(L, Lbra, H, ket, i, R, Rbra)
    cap.append((L,Lbra,H,ket,i,R,Rbra,out)); return out
tdvp.two_site_heff = spy
orig_e = tdvp._arnoldi_expm_multiply
tdvp._arnoldi_expm_multiply = lambda *a, **k: a[1].copy()   # skip evolve, only capture
tdvp.tdvp_step(psi, G, 1j*0.5, cutoff=1e-10, maxdim=64, niter=50, hermitian=False)
big=[c for c in cap if c[-1][3].size>100000][0]
L,Lbra,H,ket,i,R,Rbra,(matvec,order_in,shape,x0)=big
print("bond", i, "shape", shape, "L", L.array.shape, "R", R.array.shape, "Hi", H.A(i).array.shape, "Hj", H.A(i+1).array.shape)
import pickle
x=x0.copy()
def tm(f, n=5):
    f(); t=time.time()
    for _ in range(n): f()
    return (time.time()-t)/n
print("matvec ms", 1e3*tm(lambda: matvec(x)))
# plan steps
steps, fp = kernels._plan_contraction_chain([p for p in (L,H.A(i),H.A(i+1),R)], order_in, [Lbra, ket.A(i).inds[0], ket.A(i).inds[0], Rbra] if False else order_in)
