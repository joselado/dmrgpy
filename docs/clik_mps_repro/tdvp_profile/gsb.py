import numpy as np, time
rng=np.random.default_rng(1)
cr=lambda *s: rng.standard_normal(s)+1j*rng.standard_normal(s)
n=102400; k=35
Qc=np.empty((n,50),complex); Qc[:]=cr(n,50)       # current layout (n,m), view [:, :k]
Qr=np.ascontiguousarray(cr(50,n))                 # row-major (m,n)
w=cr(n)
def cur():
    Q=Qc[:,:k]; h=Q.conj().T@w; u=w-Q@h; h2=Q.conj().T@u; u=u-Q@h2
def rowmajor_2pass():
    Q=Qr[:k]; h=(Q@w.conj()).conj(); u=w-h@Q; h2=(Q@u.conj()).conj(); u=u-h2@Q
def rowmajor_1pass():
    Q=Qr[:k]; h=(Q@w.conj()).conj(); u=w-h@Q
def tm(f,n=10):
    f(); t=time.time()
    for _ in range(n): f()
    return 1e3*(time.time()-t)/n
for f in (cur,rowmajor_2pass,rowmajor_1pass): print(f.__name__, "%.1f ms"%tm(f))
