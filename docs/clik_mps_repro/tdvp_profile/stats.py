from common import *
from dmrgpy.pyitensor import kernels
S, G, psi = lindblad_case()
# sparsity of MPO tensors
for i in range(1, G.length()+1):
    a = G.A(i).array
    print(i, a.shape, a.dtype, "nnz frac %.3f"%(np.count_nonzero(np.abs(a)>1e-14)/a.size), "imag?", np.abs(a.imag).max()>1e-14)
runs=[]
orig = tdvp._arnoldi_run
def spy(matvec, basis, coeff, m, n, tol, errgoal):
    t=time.time(); r = orig(matvec, basis, coeff, m, n, tol, errgoal)
    runs.append((n, r[3], r[4], coeff, time.time()-t)); return r
tdvp._arnoldi_run = spy
tdvp.tdvp_step(psi, G, 1j*0.5, cutoff=1e-10, maxdim=64, niter=50, hermitian=False)
tdvp.tdvp_step(psi, G, 1j*0.5, cutoff=1e-10, maxdim=64, niter=50, hermitian=False)
for r in runs[-34:]: print("n=%d k=%d conv=%s coeff=%s t=%.2f ms/mv=%.1f"%(r[0],r[1],r[2],r[3],r[4], 1e3*r[4]/r[1]))
print("total runs", len(runs), "unconverged-first-pass", sum(1 for r in runs if not r[2]))
