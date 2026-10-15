# Local spin correlator S^zz_00(omega) of the open spin-1/2 Heisenberg chain
# from every DMRG submode, each against the exact Lehmann sum with the same
# broadening.
import numpy as np
from dmrgpy import spinchain

n = 10
delta = 0.1
es = np.linspace(-0.5, 3.5, 300)
es_cvm = np.linspace(-0.5, 3.5, 100)  # CVM solves one linear system per point

sc = spinchain.Spin_Chain(["S=1/2"]*n)
sc.set_hamiltonian(sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                       + sc.Sz[i]*sc.Sz[i+1] for i in range(n-1)))
sc.maxm, sc.nsweeps, sc.cutoff = 60, 14, 1e-12
sc.kpmmaxm, sc.cvm_maxm = 60, 60
pair = (sc.Sz[0], sc.Sz[0])


def spectrum(**kw):
    return np.asarray(sc.get_dynamical_correlator(name=pair, delta=delta,
                                                  **kw)[1]).real


ed = spectrum(mode="ED", submode="ED", es=es)
ed_cvm = spectrum(mode="ED", submode="ED", es=es_cvm)
ed_kpm = spectrum(mode="ED", submode="KPM", es=es)  # the same kernel as KPM

runs = {"KPM": dict(es=es), "CVM": dict(es=es_cvm),
        "TD": dict(es=es, dt=0.05), "TDZ": dict(es=es, dt=0.05),
        "EX": dict(es=es, nex=20), "SECTOR": dict(es=es, nex=30)}
y = {"y_" + m: spectrum(submode=m, **kw) for m, kw in runs.items()}

np.savez("data.npz", methods=list(runs), es=es, es_cvm=es_cvm, ed=ed,
         ed_cvm=ed_cvm, ed_kpm=ed_kpm, **y)
