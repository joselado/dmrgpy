import sys, time; sys.path.insert(0, sys.path[0])
from cases import *
from dmrgpy import bosonchain
from dmrgpy.lindblad import Lindblad_Chain
from dmrgpy.pyitensor import clik
kappa, dt, meter = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3] == "on"
L, N = 3, 2
bc = bosonchain.Bosonic_Chain(L, maxnb=[N+1]*L, itensor_version="python")
bc.set_hamiltonian(sum(-(bc.Adag[j+1]*bc.A[j] + bc.Adag[j]*bc.A[j+1]) for j in range(L-1)) + 0.5*sum(bc.N[j]*(bc.N[j]-1) for j in range(L)))
bc.maxm = 64
lc = Lindblad_Chain(bc, [np.sqrt(kappa)*(bc.Adag[j+1]+bc.Adag[j])*(bc.A[j+1]-bc.A[j]) for j in range(L-1)])
if not meter:
    clik._step_accepted = lambda *a: True
_eb = clik.evolve_branch
def eb(*a, **k):
    st, n = _eb(*a, **k); print("  branch done, substeps", n, flush=True); return st, n
clik.evolve_branch = eb
t = time.time()
lam, res, obs = lc.get_spectrum(nk=20, dt=dt, alpha=0.2, start=bh_start(L, N))
t1 = time.time() - t
bc.set_conserved_sector(Nb=N)
lamE, _, obsE = lc.get_spectrum(mode="ED")
print("kappa", kappa, "dt", dt, "meter", meter, "time %.1f" % t1)
print(" err", np.array2string(errors(lam, lamE, 6), precision=1), " n err %.1e" % max(abs(obs(bc.N[j], 0) - obsE(bc.N[j], 0)) for j in range(L)))
