import sys, time; sys.path.insert(0, sys.path[0])
from cases import *
from dmrgpy.pyitensor import tdvp
eg = float(sys.argv[2]) if len(sys.argv) > 2 else None
if eg:
    d = list(tdvp._arnoldi_expm_multiply.__defaults__); d[-1] = eg
    tdvp._arnoldi_expm_multiply.__defaults__ = tuple(d)
ks = []; orig = tdvp._arnoldi_run
def run(*a, **k):
    out = orig(*a, **k); ks.append(out[3]); return out
tdvp._arnoldi_run = run
case = sys.argv[1]
if case == "xxz":
    sc, lc = xxz(); kw = {}
else:
    sc, lc = bh(); kw = dict(nk=20, dt=0.5, start=bh_start(3, 2))
t = time.time(); lam, res, obs = lc.get_spectrum(**kw); t1 = time.time() - t
if case == "bh": sc.set_conserved_sector(Nb=2)
lamE, _, obsE = lc.get_spectrum(mode="ED")
print(case, "errgoal", eg, "time %.1f" % t1, "exps", len(ks), "mean k %.1f max %d" % (np.mean(ks), max(ks)))
print(" err", np.array2string(errors(lam, lamE), precision=1))
print(" obs err %.1e" % max(abs(obs(sc.Sz[i] if case=="xxz" else sc.N[i], 0) - obsE(sc.Sz[i] if case=="xxz" else sc.N[i], 0)) for i in range(len(sc.sites))))
