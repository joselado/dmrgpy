from common import *
from dmrgpy.pyitensor import kernels, dmrg
S, G, psi = lindblad_case()
for _ in range(3): tdvp.tdvp_step(psi, G, 1j*0.5, cutoff=1e-10, maxdim=64, niter=50, hermitian=False)
cap=[]
orig = dmrg.two_site_heff
def spy(*a):
    out = orig(*a); cap.append(out); return out
tdvp.two_site_heff = spy
tdvp._arnoldi_expm_multiply = lambda *a, **k: a[1].copy()
tdvp.tdvp_step(psi, G, 1j*0.5, cutoff=1e-10, maxdim=64, niter=50, hermitian=False)
big=[c for c in cap if c[3].size>100000]
matvec,order_in,shape,x0=big[1]
tdvp._arnoldi_expm_multiply = tdvp.__dict__['_arnoldi_expm_multiply']
import importlib; importlib.reload(tdvp)
ref=tdvp._arnoldi_expm_multiply(matvec,x0,0.25,niter=120,errgoal=1e-14,tol=1e-14)
print("|x0|",np.linalg.norm(x0),"|ref|",np.linalg.norm(ref))
cnt=[0]
def mv2(x): cnt[0]+=1; return matvec(x)
for eg in (1e-12,1e-10,1e-8,1e-6,1e-4):
    cnt[0]=0; t=time.time(); y=tdvp._arnoldi_expm_multiply(mv2,x0,0.25,niter=50,errgoal=eg); dt=time.time()-t
    print("errgoal %.0e matvecs %d  time %.2fs  rel err %.2e"%(eg,cnt[0],dt,np.linalg.norm(y-ref)/np.linalg.norm(ref)))
# ritz-ish spectral extent via Hessenberg of 30 steps
