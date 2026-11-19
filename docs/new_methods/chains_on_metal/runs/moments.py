"""Spectral moments of the tip vertices on the 44-site dimer: exact from the
MPS (b = X|GS>, m_k = <b|(H-E0)^k|b>) against the TASK poles of each route."""
import sys, os, warnings, time
sys.path.insert(0, os.environ.get("DMRGPY_SRC", "../../../../src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from dmrgpy import ysr, taskdmrg
from dimer import LAMBDA, N, JK
warnings.simplefilter("ignore")
order = sys.argv[1]
kw = {} if order == "linear" else dict(expand_bonds=8, expand_rounds=3)
fc = ysr.spin_chain_on_metal(2, N, site=[0, 1], Js=0.3, JK=JK, Lambda=LAMBDA, order=order)
fc.maxm, fc.nsweeps, fc.cutoff = 128, 30, 1e-10
t0 = time.time(); e0 = float(np.real(fc.gs_energy())); gs = fc.get_gs()
H = fc.hamiltonian
for X in ysr.tunneling_operators(fc, imp=fc.impurities[0]):
    for lab, (A, B) in (("e", (X.get_dagger(), X)), ("h", (X, X.get_dagger()))):
        b = fc.applyoperator(B, gs)
        Hb = fc.applyoperator(H, b)
        m0 = np.real(b.dot(b)); hb = np.real(b.dot(Hb)); hh = np.real(Hb.dot(Hb))
        m1 = hb - e0*m0; m2 = hh - 2*e0*hb + e0**2*m0
        r = taskdmrg.get_poles(fc, name=(A, B), nkry=100, **kw)
        w, wt = r["omega"], np.real(r["weight"]); el = np.real(r["elastic"])
        t0m, t1, t2 = np.sum(wt) + el, np.sum(wt*w), np.sum(wt*w**2)
        print(order, lab, "m0 %.6f/%.6f m1 %.6f/%.6f m2 %.6f/%.6f (exact/TASK)" % (m0, t0m, m1, t1, m2, t2), flush=True)
print("E0 %.10f t=%.0f" % (e0, time.time() - t0))
