from common import *
import cProfile, pstats
S, G, psi = lindblad_case()
print("MPO bond dims", [ind.dim for i in range(1, G.length()) for ind in G.A(i).inds if ind.hastags("Link") and ind in G.A(i+1).inds] if hasattr(G,'length') else None)
print("site dims", [ind.dim for ind in psi.A(1).inds])
print("herm", mpsalgebra.mpo_is_hermitian(G))
t=time.time(); tdvp.tdvp_step(psi, G, 1j*0.5, cutoff=1e-10, maxdim=64, niter=50, hermitian=False); print("step1", time.time()-t)
print("bonds", [ [ind.dim for ind in psi.A(i).inds] for i in range(1,psi.length()+1)])
pr = cProfile.Profile(); pr.enable()
for _ in range(2):
    tdvp.tdvp_step(psi, G, 1j*0.5, cutoff=1e-10, maxdim=64, niter=50, hermitian=False)
pr.disable()
print("bonds", [ [ind.dim for ind in psi.A(i).inds] for i in range(1,psi.length()+1)])
pstats.Stats(pr).sort_stats("tottime").print_stats(18)
