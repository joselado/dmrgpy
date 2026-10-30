"""Two-site TDVP until the bonds stop changing, then one-site TDVP (no expansion),
with a two-site step every `refresh` steps. Accuracy and time against pure two-site."""
import sys, time; sys.path.insert(0, sys.path[0])
from cases import *
from dmrgpy.pyitensor import clik, backend as bk
from dmrgpy.pyitensor.tdvp import tdvp_step
from dmrgpy.pyitensor.mpsalgebra import mpo_is_hermitian
from dmrgpy.pyitensor.mpscontainer import _link_at
from dmrgpy.pyitensor.gse import global_subspace_expand
STATS = {"one": 0, "two": 0}
def bonds(x): return [_link_at(x, i, i+1).dim for i in range(1, x.length())]
def evolve_branch(G, psi0, z, nk, cutoff, maxdim, niter=50, hermitian=None, one_site=False, refresh=REFRESH if False else None):
    refresh = REFRESH
    if hermitian is None: hermitian = mpo_is_hermitian(G)
    psi = psi0.copy(); psi.normalize(); states = [psi]; saturated = False; since = 0
    for m in range(nk - 1):
        if m < clik._GSE_STEPS:
            with bk.pad_bonds_suspended():
                psi = global_subspace_expand(G, psi, clik._GSE_ORDER, cutoff, maxdim=maxdim, bond_maxdim=maxdim)
        use_one = saturated and refresh > 0 and since % refresh != refresh - 1
        b0 = bonds(psi)
        if use_one:
            out = tdvp_step(psi.copy(), G, 1j*z, cutoff=cutoff, maxdim=maxdim, niter=niter, num_center=1, hermitian=hermitian); STATS["one"] += 1; since += 1
        else:
            out = tdvp_step(psi.copy(), G, 1j*z, cutoff=cutoff, maxdim=maxdim, niter=niter, hermitian=hermitian); STATS["two"] += 1
            if m >= clik._GSE_STEPS:
                saturated = bonds(out) == b0
                since = 0
        out.normalize(); psi = out; states.append(psi)
    return states
case = sys.argv[1]; REFRESH = int(sys.argv[2])
if REFRESH > 0: clik.evolve_branch = evolve_branch
if case == "xxz":
    sc, lc = xxz(); kw = {}; ops = [sc.Sz[i] for i in range(4)]
elif case == "bh":
    sc, lc = bh(); kw = dict(nk=20, dt=0.5, start=bh_start(3, 2)); ops = [sc.N[i] for i in range(3)]
t = time.time(); lam, res, obs = lc.get_spectrum(**kw); t1 = time.time() - t
if case == "bh": sc.set_conserved_sector(Nb=2)
lamE, _, obsE = lc.get_spectrum(mode="ED")
print(case, "refresh", REFRESH, "time %.1f" % t1, STATS)
print(" err", np.array2string(errors(lam, lamE, 8), precision=1))
print(" obs err %.1e" % max(abs(obs(o, 0) - obsE(o, 0)) for o in ops))
