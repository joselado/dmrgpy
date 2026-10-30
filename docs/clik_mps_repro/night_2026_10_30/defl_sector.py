import sys, time, itertools; sys.path.insert(0, sys.path[0])
from cases import *
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.mpsalgebra import sum_many, inner
from dmrgpy.pyitensor.metts import product_state
from dmrgpy.pyitensor.vectorize import identity_mps
mode = sys.argv[1]
L, N = int(sys.argv[2]), int(sys.argv[3]); dt = 0.5; nk = 20
bc, lc = bh(L, N)
dc = lc.doubled; S = dc._session
G = S._mpo(lc.L.to_terms())
def basis(k, j):
    v = np.zeros(S.sites.dim(k), complex); v[j] = 1; return v
from sector_id import sector_identity_mps
oneN = sector_identity_mps(S.sites, N)
p, u = lc._seeds(*bh_start(L, N))
starts = (product_state(S.sites, p), product_state(S.sites, u))
d = {"none": None, "full": identity_mps(S.sites), "sector": oneN}[mode]
from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor.mpsalgebra import sum as mps_sum
def eb(G, psi0, z, nk, cutoff, maxdim, niter=50, hermitian=None, one_site=False, deflate_=None):
    dfl = d if z != dt else None
    psi = psi0.copy(); psi.normalize(); states = [psi]
    for m in range(nk - 1):
        if m < clik._GSE_STEPS:
            with bk.pad_bonds_suspended():
                psi = clik.global_subspace_expand(G, psi, clik._GSE_ORDER, cutoff, maxdim=maxdim, bond_maxdim=maxdim)
        psi, _ = clik._step(psi, G, z, cutoff, maxdim, niter, hermitian, False)
        if dfl is not None:
            c = inner(dfl, psi) / inner(dfl, dfl)
            psi = mps_sum(psi, (-c) * dfl, cutoff=cutoff, maxdim=maxdim); psi.normalize()
        states.append(psi)
    return states
clik.evolve_branch = eb
t = time.time()
ks = clik.krylov_spectrum(G, S.sites, nk=nk, dt=dt, alpha=None, maxdim=64, cutoff=1e-10, starts=starts)
t1 = time.time() - t
tr = [abs(inner(oneN, s)) for s in ks.states[nk:]]
bc.set_conserved_sector(Nb=N)
lamE = lc.get_spectrum(mode="ED")[0]
print(mode, L, N, "time %.1f" % t1, "alpha %.3f" % ks.alpha, "traces branch 2", np.array2string(np.array(tr[::5]), precision=1))
print(" err", np.array2string(errors(ks.eigenvalues, lamE, 8), precision=1))
