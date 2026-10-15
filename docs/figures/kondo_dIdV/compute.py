# STM dI/dV of magnetic atoms and chains, in the third-order perturbation
# theory of Ternes (NJP 17, 063016). Energies in eV, temperatures in kelvin.
import numpy as np
from dmrgpy import spinchain

muB = 5.7883818066e-5  # eV/T

# (a) inelastic steps of an S=2 atom, D Sz^2 + E (Sx^2 - Sy^2) (Fe on Cu2N)
D, E = -1.55e-3, 0.31e-3
sc = spinchain.Spin_Chain(["S=2"])
sc.set_hamiltonian(D*sc.Sz[0]*sc.Sz[0] + E*(sc.Sx[0]*sc.Sx[0] - sc.Sy[0]*sc.Sy[0]))
eV_a, dIdV_a = sc.get_kondo_spectrum(np.linspace(-10e-3, 10e-3, 401), site=0,
                                     order=2, T=0.5)
# reference: excitations of the exact 5x5 Hamiltonian that S+, S-, Sz reach
m = np.arange(2.0, -3.0, -1.0)
sz, sp = np.diag(m), np.diag(np.sqrt(6 - m[1:]*(m[1:] + 1)), 1)
sx, sy = (sp + sp.T)/2, (sp - sp.T)/2j
w, v = np.linalg.eigh((D*sz @ sz + E*(sx @ sx - sy @ sy)).real)
weight = sum(f*np.abs(v.T @ op @ v[:, 0])**2
             for f, op in ((0.5, sp), (0.5, sp.T), (1.0, sz)))
eps_a = (w - w[0])[1:][weight[1:] > 1e-8]

# (b) Kondo peak of an S=1/2 in a field B (in tesla), in units of e^2 T0^2/h
B_b = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
eV_b = np.linspace(-4e-3, 4e-3, 161)
dIdV_b = []
for B in B_b:
    sc = spinchain.Spin_Chain(["S=1/2"])
    sc.set_hamiltonian(2.0*muB*B*sc.Sz[0])
    _, y = sc.get_kondo_spectrum(eV_b, site=0, Jrho_s=-0.05, U=0.0, T=1.0,
                                 order=3, omega0=20e-3)
    dIdV_b.append(y/(2*np.pi))

# (c) 3-site Heisenberg chain with a field on the first site, at T=0:
# DMRG (KPM and a two-time TDVP construction) against ED
sc = spinchain.Spin_Chain(["S=1/2"]*3)
sc.maxm, sc.nsweeps = 30, 15
sc.set_hamiltonian(2.0*muB*10.0*sc.Sz[0]
                   + sum(0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                               + sc.Sz[i]*sc.Sz[i+1]) for i in range(2)))
sc.get_gs()
eV_c = np.linspace(-2e-3, 2e-3, 41)
kw = dict(site=0, Jrho_s=0.05, U=0.0, T=0.0, order=3, omega0=2e-3, Gamma0=2e-3)
_, dIdV_ed = sc.get_kondo_spectrum(eV_c, mode="ED", **kw)
_, dIdV_dmrg = sc.get_kondo_spectrum(
    eV_c, mode="DMRG", submode="KPM", delta=2e-5,
    es=np.linspace(-3e-3, 3e-3, 800), dt2=150.0, n_t2_half=20, dtau=150.0,
    n_tau_half=20, **kw)

np.savez("data.npz", eV_a=eV_a, dIdV_a=dIdV_a, eps_a=eps_a, B_b=B_b, eV_b=eV_b,
         dIdV_b=dIdV_b, eV_c=eV_c, dIdV_ed=dIdV_ed, dIdV_dmrg=dIdV_dmrg)
