"""Does one-site TDVP + GSE leave the (N_ket, N_bra) sector? Seeded 3-site BH."""
import sys; sys.path.insert(0, sys.path[0])
from cases import *
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.mpsalgebra import inner
from dmrgpy.pyitensor.metts import product_state
from dmrgpy.pyitensor.mpscontainer import _link_at
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
