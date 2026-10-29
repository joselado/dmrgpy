import sys; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np, scipy.linalg as sl
from dmrgpy.pyitensor import tdvp
def dense(matvec, v0, coeff, niter=30, **kw):
    n=v0.size; A=np.zeros((n,n),dtype=complex)
    for j in range(n):
        e=np.zeros(n,dtype=complex); e[j]=1; A[:,j]=matvec(e)
    return sl.expm(coeff*A)@v0
if sys.argv[1]=="dense": tdvp._lanczos_expm_multiply = dense
exec(open("nh_tdvp2.py").read().split("import sys as S")[0])
for G in (0, 0.7j):
    a=np.array(run("DMRG",G)[1]); print(G, np.round(a[::8].real,6))
