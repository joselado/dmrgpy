"""Exact Lorentzian dI/dV of an 8-site dimer on two metals from sparse ED
(Lanczos ground state, Lanczos continued fraction per vertex), against the
two TASK routes at truncating bond dimensions."""
import sys, time, warnings; import os; sys.path.insert(0, os.environ.get("DMRGPY_SRC", "../../../../src"))
import numpy as np, scipy.sparse.linalg as sla
from dmrgpy import ysr
warnings.simplefilter("ignore")
Js, n = float(sys.argv[1]), 3
JK = 0.3*ysr.lambda_factor(2.0)
es = np.linspace(-0.6, 0.6, 121); eta = 0.02
def chain(order, site, version=3):
    return ysr.spin_chain_on_metal(2, n, site=site, Js=Js, JK=JK, Lambda=2.0, order=order, itensor_version=version)
t0 = time.time()
ref = chain("linear", [0, 1])
ed = ref.get_ED_obj()
H = ed.get_operator(ref.hamiltonian).tocsr()
e0, v = sla.eigsh(H, k=3, which="SA", v0=np.ones(H.shape[0]), tol=1e-12); print("lowest three", e0)
e0, gs = e0[0], v[:, 0]
print("dim", H.shape[0], "E0", e0, "t=%.0f" % (time.time() - t0), flush=True)
def cf(b, z, m=400):
    nb = np.linalg.norm(b)
    if nb == 0: return np.zeros(len(z))
    q, qold, a, bb = b/nb, np.zeros_like(b), [], []
    for k in range(m):
        w = H @ q - (bb[-1]*qold if bb else 0)
        al = np.vdot(q, w).real; w = w - al*q
        a.append(al); be = np.linalg.norm(w)
        if be < 1e-12: break
        bb.append(be); qold, q = q, w/be
    g = np.zeros(len(z), complex)
    for k in reversed(range(len(a))):
        g = 1.0/(z - a[k] - (bb[k]**2*g if k < len(bb) else 0))
    return nb**2*g
A = np.zeros(len(es))
for X in ysr.tunneling_operators(ref, imp=ref.impurities[0]):
    Xm = ed.get_operator(X).tocsr(); Xd = ed.get_operator(X.get_dagger()).tocsr()
    ge = cf(Xm @ gs, es + e0 + 1j*eta)
    gh = cf(Xd @ gs, -es + e0 + 1j*eta)
    A += -ge.imag/np.pi - gh.imag/np.pi
print("exact peak %.4f t=%.0f" % (A.max(), time.time() - t0), flush=True)
for order, kw in (("linear", {}), ("sides", dict(expand_bonds=8, expand_rounds=3))):
    for m in (16, 32, 64):
        fc = chain(order, [0, 1]); fc.maxm, fc.nsweeps = m, 20
        d = ysr.get_stm_spectrum(fc, es, delta=eta, imp=fc.impurities[0], nkry=200, **kw)[1]
        print("Js %.2f %s maxm %d: max|TASK-exact|/peak %.3e" % (Js, order, m, np.max(np.abs(d - A))/A.max()), flush=True)
