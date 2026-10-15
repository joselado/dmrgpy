# Local density of states A_i(omega) of a half-filled open chain of spinless
# fermions, at a bulk site and at the edge, for V=0 and V=3t (charge density
# wave), with the exact free-fermion spectrum as the reference at V=0.
import numpy as np
from dmrgpy import fermionchain

n = 20
N = n//2  # half filling
Vs = [0.0, 3.0]
sites = [n//2, 0]
nex, delta = 50, 0.1
es = np.linspace(-5.0, 5.0, 1000)

# free fermions on the same open chain, Lorentzian broadening delta
eps, phi = np.linalg.eigh(-(np.eye(n, k=1) + np.eye(n, k=-1)))
mu_free = (eps[N-1] + eps[N])/2
A_free = [sum(phi[i, k]**2*delta/np.pi/((es - eps[k] + mu_free)**2 + delta**2)
              for k in range(n)) for i in sites]

A, captured, gap = [], [], []
for V in Vs:
    fc = fermionchain.Fermionic_Chain(n)
    h = sum(-fc.Cdag[i]*fc.C[i+1] - fc.Cdag[i+1]*fc.C[i]
            + V*(fc.N[i] - 0.5)*(fc.N[i+1] - 0.5) for i in range(n-1))
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps, fc.cutoff = 80, 12, 1e-12
    fc.set_conserved_sector(Nf=N)
    for i in sites:  # the same chain, so the sector states are reused
        _, a, info = fc.get_spectral_function(i=i, nex=nex, es=es, delta=delta,
                                              return_poles=True)
        A.append(np.real(a))
        captured.append([info["particle"]["captured"], info["hole"]["captured"]])
    gap.append(info["gap"])

np.savez("data.npz", Vs=Vs, sites=sites, es=es,
         A=np.reshape(A, (len(Vs), len(sites), len(es))),
         captured=np.reshape(captured, (len(Vs), len(sites), 2)),
         gap=gap, A_free=A_free)
