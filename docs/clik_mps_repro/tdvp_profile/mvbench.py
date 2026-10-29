from common import *
from mpocomp import mpo_arrays
S, G, psi = lindblad_case()
Ws=mpo_arrays(G)
W1,W2=Ws[3],Ws[4]    # (14,5,5,22),(22,5,5,14)
for W in (W1,W2):
    nzp=np.abs(W).max(axis=(1,2))>1e-14
    print(W.shape,"nonzero (a,b) pairs",nzp.sum(),"of",nzp.size)
    # operator rank
    M=W.transpose(0,3,1,2).reshape(-1,25)
    s=np.linalg.svd(M,compute_uv=False); print(" operator-space rank", int((s>1e-12*s[0]).sum()), "sv", np.round(s[:12],3))
chi=64; d=5
rng=np.random.default_rng(1)
cr=lambda *s: rng.standard_normal(s)+1j*rng.standard_normal(s)
L=cr(chi,14,chi); R=cr(chi,14,chi); v=cr(chi,d,d,chi)
def tm(f,n=5):
    f(); t=time.time()
    for _ in range(n): f()
    return 1e3*(time.time()-t)/n
# V0 stepwise
def v0():
    x=L.reshape(chi*14,chi)@v.reshape(chi,d*d*chi)         # (bra,a,si,sj,r)
    x=x.reshape(chi,14,d,d,chi)
    x=np.tensordot(x,W1,axes=([1,2],[0,1]))                  # (bra,sj,r,s',b)
    x=np.tensordot(x,W2,axes=([4,1],[0,1]))                  # (bra,r,s_i',s_j',c)
    return np.tensordot(x,R,axes=([4,1],[1,2]))             # (bra,si',sj',rbra)
def steps():
    t0=time.time(); x=L.reshape(chi*14,chi)@v.reshape(chi,d*d*chi); x=x.reshape(chi,14,d,d,chi); t1=time.time()
    x=np.tensordot(x,W1,axes=([1,2],[0,1])); t2=time.time()
    x=np.tensordot(x,W2,axes=([4,1],[0,1])); t3=time.time()
    y=np.tensordot(x,R,axes=([4,1],[1,2])); t4=time.time()
    return [1e3*(t1-t0),1e3*(t2-t1),1e3*(t3-t2),1e3*(t4-t3)]
v0(); print("V0 total ms %.1f"%tm(v0)); 
st=np.mean([steps() for _ in range(5)],axis=0); print("steps ms L,W1,W2,R:",np.round(st,1))
