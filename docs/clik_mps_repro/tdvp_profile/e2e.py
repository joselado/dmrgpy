from common import *
import fastarn
from dmrgpy.pyitensor.mpsalgebra import inner
S, G, psi0 = lindblad_case()
def run(label, expm_fn, nsteps=3):
    psi = psi0.copy(); ts=[]
    if expm_fn is not None: tdvp._arnoldi_expm_multiply = expm_fn
    for _ in range(nsteps):
        t=time.time(); tdvp.tdvp_step(psi, G, 1j*0.5, cutoff=1e-10, maxdim=64, niter=50, hermitian=False)
        psi.normalize(); ts.append(time.time()-t)
    print(label, "step times", np.round(ts,1), "total %.1f"%sum(ts), flush=True)
    return psi
base_fn = tdvp._arnoldi_expm_multiply
a = run("baseline", None)
b = run("fast(selective)", fastarn.fast_arnoldi_expm_multiply)
print("second-pass fraction", fastarn.STATS["second_pass"]/fastarn.STATS["iters"], flush=True)
fastarn.SELECTIVE=False
c = run("fast(2-pass)", fastarn.fast_arnoldi_expm_multiply)
def cmp(x,y,name):
    ov = inner(x,y); nx=inner(x,x); ny=inner(y,y)
    Ex=inner(x,G,x)/nx; Ey=inner(y,G,y)/ny
    print(name, "1-|<x|y>|/norms = %.2e"%(1-abs(ov)/np.sqrt(abs(nx)*abs(ny))), "<G> base %.10f fast %.10f diff %.1e"%(Ex.real, Ey.real, abs(Ex-Ey)))
cmp(a,b,"selective"); cmp(a,c,"2-pass")
