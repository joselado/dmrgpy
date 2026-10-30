import sys, time; sys.path.insert(0, sys.path[0])
from cases import *
from dmrgpy.pyitensor import clik, tdvp, backend as bk
from dmrgpy.pyitensor.mpsalgebra import inner
from dmrgpy.pyitensor.metts import product_state
L, N = int(sys.argv[1]), int(sys.argv[2])
bc, lc = bh(L, N)
S = lc.doubled._session; G = S._mpo(lc.L.to_terms())
p, u = lc._seeds(*bh_start(L, N))
log = []
orig = tdvp._arnoldi_run
def run(matvec, basis, coeff, m, n, tol, errgoal):
    out = orig(matvec, basis, coeff, m, n, tol, errgoal)
    c = out[5]; log.append((out[3], np.linalg.norm(c)))
    return out
tdvp._arnoldi_run = run
res = {}
from sector_id import sector_identity_mps
oneN = sector_identity_mps(S.sites, N)
dc = lc.doubled
nk_ = sum(dc.N[2*i] for i in range(L)); nb_ = sum(dc.N[2*i+1] for i in range(L))
Q = S._mpo((((nk_ - N)*(nk_ - N)) + ((nb_ - N)*(nb_ - N))).to_terms())
EGS = [float(x) for x in sys.argv[3].split(",")]
for eg in EGS:
    d = list(tdvp._arnoldi_expm_multiply.__defaults__); d[-1] = eg
    tdvp._arnoldi_expm_multiply.__defaults__ = tuple(d)
    psi = product_state(S.sites, p); psi.normalize(); out = []; log.clear()
    for m in range(3):
        psi = clik.global_subspace_expand(G, psi, 3, 1e-10, maxdim=64, bond_maxdim=64)
        psi, ln = clik._step(psi, G, 0.5, 1e-10, 64, 50, False)
        out.append((psi, ln))
        print("   eg %.0e step %d: lognorm %.3f trace %.3e sector dev %.1e" % (eg, m, ln, abs(inner(oneN, psi)), abs(inner(psi, Q, psi))), flush=True)
    ks = np.array([k for k, _ in log]); cn = np.array([c for _, c in log])
    print("errgoal %.0e: k mean %.1f min %d; ||c|| min %.1e; lognorms %s" % (eg, ks.mean(), ks.min(), cn.min(), [round(x[1], 3) for x in out]), flush=True)
    res[eg] = out
for eg in EGS[:-1]:
    print(eg, "1-|<ref|x>| per step:", ["%.1e" % (1 - abs(inner(res[EGS[-1]][m][0], res[eg][m][0]))) for m in range(3)])
