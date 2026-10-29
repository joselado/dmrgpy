import numpy as np, time
chi=64; d=5; rng=np.random.default_rng(1)
cr=lambda *s: rng.standard_normal(s)+1j*rng.standard_normal(s)
L=cr(chi,14,chi); R=cr(chi,14,chi); v=cr(chi,d,d,chi); W1=cr(14,5,5,22); W2=cr(22,5,5,14)
def mv():
    x=L.reshape(chi*14,chi)@v.reshape(chi,d*d*chi); x=x.reshape(chi,14,d,d,chi)
    x=np.tensordot(x,W1,axes=([1,2],[0,1])); x=np.tensordot(x,W2,axes=([4,1],[0,1]))
    return np.tensordot(x,R,axes=([4,1],[1,2]))
Q=cr(102400,35); w=cr(102400)
def gs():
    h=Q.conj().T@w; u=w-Q@h; h2=Q.conj().T@u; u=u-Q@h2
def tm(f,n=5):
    f(); t=time.time()
    for _ in range(n): f()
    return 1e3*(time.time()-t)/n
print("matvec ms %.1f  gram-schmidt(k=35) ms %.1f"%(tm(mv),tm(gs)))
