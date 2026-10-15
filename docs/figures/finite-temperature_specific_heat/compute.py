# Thermal energy and specific heat of the spin-1/2 Heisenberg chain from
# purification and from METTS, against the exact full-spectrum result.
import numpy as np
from dmrgpy import spinchain, thermal

n = 10
T_pur = np.geomspace(0.1, 5, 20)
T_metts = [0.25, 0.5, 1.0, 2.0, 4.0]
T_exact = np.geomspace(0.05, 10, 300)


def heisenberg(c):
    return sum(c.Sx[i]*c.Sx[i+1] + c.Sy[i]*c.Sy[i+1] + c.Sz[i]*c.Sz[i+1]
               for i in range(n-1))


# exact: <H>_T and C = (<H^2>_T - <H>_T^2)/T^2 from all 1024 eigenvalues
sc = spinchain.Spin_Chain(["S=1/2"]*n)
h = heisenberg(sc)
sc.set_hamiltonian(h)
hm = sc.get_ED_obj().get_hamiltonian()
ev = np.linalg.eigvalsh(hm.toarray() if hasattr(hm, "toarray") else hm)
E_exact, C_exact = [], []
for T in T_exact:
    w = np.exp(-(ev - ev.min())/T)
    e, e2 = (ev*w).sum()/w.sum(), (ev**2*w).sum()/w.sum()
    E_exact.append(e)
    C_exact.append((e2 - e**2)/T**2)

# purification, a fresh thermal chain per temperature
E_pur, C_pur = [], []
for T in T_pur:
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*n, T=T)
    tc.MBChain.maxm = 100
    ht = heisenberg(tc)
    tc.set_hamiltonian(ht)
    tc.get_gs()
    E_pur.append(tc.MBChain.vev(ht).real)
    C_pur.append(tc.MBChain.gs_energy_fluctuation()**2/T**2)

# METTS on the physical chain
E_metts, err_metts = [], []
for T in T_metts:
    e, err = sc.metts_vev(h, T, nsamples=300, nwarmup=30, seed=1)
    E_metts.append(e.real)
    err_metts.append(err)

np.savez("data.npz", n=n, T=T_pur, E_pur=E_pur, C_pur=C_pur,
         T_metts=T_metts, E_metts=E_metts, err_metts=err_metts,
         T_exact=T_exact, E_exact=E_exact, C_exact=C_exact)
