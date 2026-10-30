"""Does one-site TDVP + GSE leave the (N_ket, N_bra) sector? Seeded 3-site BH."""
import sys; sys.path.insert(0, sys.path[0])
from cases import *
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.mpsalgebra import inner
from dmrgpy.pyitensor.metts import product_state
from dmrgpy.pyitensor.mpscontainer import _link_at
import importlib; from dmrgpy.pyitensor import tdvp; svdmod = importlib.import_module("dmrgpy.pyitensor.svd")
if len(sys.argv) > 4 and sys.argv[4] == "rr":
    def rr_split(T, left_inds, tags="Link", orthonormal="left"):
        U, S_, V, _ = svdmod.svd(T, left_inds, cutoff=1e-26)
        if orthonormal == "left":
            u = next(i for i in U.inds if i not in left_inds)
            return U, S_ * V, u
        v = next(i for i in V.inds if i in (S_.inds))
        return U * S_, V, v
    tdvp.qr_split = rr_split
    print("rank-revealing split")
dt = float(sys.argv[1]); nsteps = int(sys.argv[2]); cutoff = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10
L, N = 3, 2
bc, lc = bh(L, N)
dc = lc.doubled; S = dc._session
G = S._mpo(lc.L.to_terms())
nk = sum(dc.N[2*i] for i in range(L)); nb = sum(dc.N[2*i+1] for i in range(L))
Q = S._mpo((((nk - N)*(nk - N)) + ((nb - N)*(nb - N))).to_terms())
p, u = lc._seeds(*bh_start(L, N))
psi0 = product_state(S.sites, u)
def bonds(psi): return [_link_at(psi, i, i+1).dim for i in range(1, psi.length())]
def rank(psi):
    x = psi.copy(); x.position(1); x.position(x.length(), cutoff=1e-12); return bonds(x)
runs = {}
for name, one in (("two", False), ("one", True)):
    psi = psi0.copy(); psi.normalize(); out = []
    for m in range(nsteps):
        if m < clik._GSE_STEPS and not one:
            psi = clik.global_subspace_expand(G, psi, clik._GSE_ORDER, cutoff, maxdim=64, bond_maxdim=64)
        psi, _ = clik._step(psi, G, dt, cutoff, 64, 50, False, one)
        out.append(psi.copy())
        dev = abs(inner(psi, Q, psi))
        print(name, m, "sector dev %.1e" % dev, "bonds", bonds(psi), "rank", rank(psi), flush=True)
    runs[name] = out
for m in range(nsteps):
    a, b = runs["two"][m], runs["one"][m]
    print("step", m, "1-|<two|one>| = %.1e" % (1 - abs(inner(a, b))))
from dmrgpy.pyitensor.vectorize import identity_mps
one_mps = identity_mps(S.sites)
ref = psi0.copy(); ref.normalize(); sub = 20
for m in range(nsteps):
    for k in range(sub):
        if m == 0 and k < 3:
            ref = clik.global_subspace_expand(G, ref, 3, 1e-12, maxdim=64, bond_maxdim=64)
        ref, _ = clik._step(ref, G, dt / sub, 1e-12, 64, 50, False, False)
    print("step", m, "ref-two %.1e" % (1 - abs(inner(ref, runs["two"][m]))), "ref-one %.1e" % (1 - abs(inner(ref, runs["one"][m]))),
          "trace two %.1e one %.1e ref %.1e" % tuple(abs(inner(one_mps, x)) for x in (runs["two"][m], runs["one"][m], ref)))
