"""Figure 15, correlator-methods_comparison: the local spin correlator
S^zz_00(omega) of the open spin-1/2 Heisenberg chain from every DMRG
submode (KPM, CVM, TD, TDZ, EX, SECTOR), each against the exact Lehmann sum
(mode="ED", submode="ED") with the same broadening.

Run with  bash docs/figures/run.sh correlator-methods_comparison

Arrays in data.npz:
    methods       (m,)        submode names, in plot order
    es            (ne,)       frequency grid, in J
    es_cvm        (nc,)       coarser grid of CVM, in J
    y_<method>    (ne,)/(nc,) real part of S_00(omega) from that submode
    ed            (ne,)       exact Lehmann sum on es (Lorentzian lines)
    ed_cvm        (nc,)       the same on es_cvm
    ed_kpm        (ne,)       mode="ED", submode="KPM" on es (same kernel as KPM)
    time_<method> ()          wall time of that submode, in s
    e0            ()          DMRG ground-state energy, total, in J
    e0_ed         ()          ED ground-state energy, total, in J
    params                    JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain

params = common.start(dict(
    model="spin-1/2 Heisenberg chain, open, J=1",
    n=int(os.environ.get("FIG_N", 10)),
    pair="(Sz[0], Sz[0])",
    maxm=60,
    nsweeps=14,
    cutoff=1e-12,
    kpmmaxm=60,
    cvm_maxm=60,
    delta=0.1,
    es=[-0.5, 3.5, 300],
    es_cvm=[-0.5, 3.5, 100],
    dt=0.05,
    nex_EX=20,
    nex_SECTOR=30,
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
n = params["n"]
delta = params["delta"]
es = np.linspace(*params["es"][:2], params["es"][2])
es_cvm = np.linspace(*params["es_cvm"][:2], params["es_cvm"][2])

sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=params["itensor_version"])
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
sc.set_hamiltonian(h)
sc.maxm = params["maxm"]
sc.nsweeps = params["nsweeps"]
sc.cutoff = params["cutoff"]
sc.kpmmaxm = params["kpmmaxm"]
sc.cvm_maxm = params["cvm_maxm"]
pair = (sc.Sz[0], sc.Sz[0])

e0 = sc.gs_energy()
e0_ed = sc.gs_energy(mode="ED")
print("E0 DMRG = %.10f  ED = %.10f" % (e0, e0_ed), flush=True)

_, ed = sc.get_dynamical_correlator(mode="ED", submode="ED", name=pair,
                                    es=es, delta=delta)
_, ed_cvm = sc.get_dynamical_correlator(mode="ED", submode="ED", name=pair,
                                        es=es_cvm, delta=delta)
_, ed_kpm = sc.get_dynamical_correlator(mode="ED", submode="KPM", name=pair,
                                        es=es, delta=delta)
ed, ed_cvm, ed_kpm = ed.real, ed_cvm.real, ed_kpm.real

runs = [
    ("KPM", dict(submode="KPM", es=es)),
    ("CVM", dict(submode="CVM", es=es_cvm)),
    ("TD", dict(submode="TD", es=es, dt=params["dt"])),
    ("TDZ", dict(submode="TDZ", es=es, dt=params["dt"])),
    ("EX", dict(submode="EX", es=es, nex=params["nex_EX"])),
    ("SECTOR", dict(submode="SECTOR", es=es, nex=params["nex_SECTOR"])),
]
out = {}
for name, kw in runs:
    t0 = time.time()
    _, y = sc.get_dynamical_correlator(name=pair, delta=delta, **kw)
    y = np.asarray(y).real
    t = time.time() - t0
    ref = ed_cvm if name == "CVM" else ed
    x = es_cvm if name == "CVM" else es
    res = np.max(np.abs(y - ref)) / np.max(np.abs(ref))
    msg = "%-7s max|y-ED|/max ED = %.3e  sum on window = %.4f (ED %.4f)" % (
        name, res, np.trapezoid(y, x), np.trapezoid(ref, x))
    if name == "KPM":
        msg += "  against ED-KPM %.3e" % (np.max(np.abs(y - ed_kpm))
                                          / np.max(np.abs(ed_kpm)))
    print(msg + "  (%.1f s)" % t, flush=True)
    out["y_" + name] = y
    out["time_" + name] = np.array(t)

common.save(__file__, params, methods=np.array([r[0] for r in runs]),
            es=es, es_cvm=es_cvm, ed=ed, ed_cvm=ed_cvm, ed_kpm=ed_kpm,
            e0=np.array(e0), e0_ed=np.array(e0_ed), **out)
