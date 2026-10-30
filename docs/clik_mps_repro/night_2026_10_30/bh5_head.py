import sys, time; sys.path.insert(0, sys.path[0])
from cases import *
maxm, nk, dt, alpha = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
L, N = 5, 4
from dmrgpy.pyitensor import clik, tdvp
if len(sys.argv) > 5:
    d = list(tdvp._arnoldi_expm_multiply.__defaults__); d[-1] = float(sys.argv[5])
    tdvp._arnoldi_expm_multiply.__defaults__ = tuple(d)
    print("errgoal", d[-1])
_add = clik.KrylovMatrices.add; T0 = time.time()
def add(self, psi):
    _add(self, psi); print("snapshot %d at %.0f s" % (len(self.states), time.time() - T0), flush=True)
clik.KrylovMatrices.add = add
bc, lc = bh(L, N, maxm)
t = time.time()
lam, res, obs = lc.get_spectrum(nk=nk, dt=dt, alpha=alpha, start=([1,1,1,1,0], [0,1,1,1,1]))
print("DMRG maxm=%d nk=%d dt=%g alpha=%g: %.0f s" % (maxm, nk, dt, alpha, time.time()-t), flush=True)
for k in range(min(8, len(lam))): print("  ", np.round(lam[k], 5), "res %.1e" % res[k], flush=True)
print("  <n_j>_ss", np.round([obs(bc.N[j], 0).real for j in range(L)], 6), flush=True)
bc.set_conserved_sector(Nb=N)
t = time.time(); lamE, _, obsE = lc.get_spectrum(mode="ED"); print("ED %.0f s" % (time.time()-t))
for k in range(8): print("  ", np.round(lamE[k], 5))
print("  <n_j>_ss", np.round([obsE(bc.N[j], 0).real for j in range(L)], 6))
print("  kept", len(lam)); print("  err", ["%.1e" % x for x in errors(lam, lamE, 8)])
