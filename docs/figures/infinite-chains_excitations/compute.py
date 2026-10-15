# Excitation bands of infinite chains from the quasiparticle ansatz on top of
# VUMPS, with the dynamical structure factor S(k,omega) of the same ansatz:
# (a) the transverse-field Ising chain, (b) the spin-1 Heisenberg chain.
import numpy as np
from dmrgpy import infinitechain


def vumps_chain(spin, h, maxm):
    ic = infinitechain.Infinite_Spin_Chain([spin], itensor_version="python")
    ic.set_hamiltonian(h(ic))
    ic.gs_method, ic.maxm = "vumps", maxm
    ic.maxiter, ic.vumps_nrestarts = 800, 6
    ic.gs_energy()
    return ic


# (a) H = -sum sx sx - g sum sz in Pauli matrices, paramagnetic phase
g = 2.5
ic = vumps_chain("1/2", lambda c: -4.0*c.SxC[0]*c.SxR[0] - 2.0*g*c.SzC[0], 4)
ks_a = np.linspace(-np.pi, np.pi, 61)
band_a = np.array([ic.excitation_energies(k, n=1)[0] for k in ks_a]).real
exact_a = 2*np.sqrt(1 + g**2 - 2*g*np.cos(ks_a))  # free-fermion band
_, es_a, S_a = ic.dynamical_structure_factor(
    "Sx", ks=ks_a, energies=np.linspace(1.5, 8.5, 400), delta=0.08, n=1)

# (b) spin-1 Heisenberg chain, 12 branches of the ansatz at each k
ich = vumps_chain("1", lambda c: c.SxC[0]*c.SxR[0] + c.SyC[0]*c.SyR[0]
                  + c.SzC[0]*c.SzR[0], 16)
ks_b = np.linspace(0, 2*np.pi, 97)
e_b, w_b = [], []
for k in ks_b:
    # energies E_a(k) and weights |<k,a|S^z_k|GS>|^2
    e, w = ich.spectral_weights("Sz", k, n=12)
    e_b.append(np.real(e))
    w_b.append(np.real(w))
e_b, w_b = np.array(e_b), np.array(w_b)
es_b, delta = np.linspace(0, 6, 400), 0.08
lorentz = (delta/np.pi)/((es_b[None, :, None] - e_b[:, None, :])**2 + delta**2)
S_b = np.einsum("kea,ka->ke", lorentz, w_b)
gap_ref = 0.4104789  # Haldane gap

np.savez("data.npz", ks_a=ks_a, band_a=band_a, exact_a=exact_a, es_a=es_a,
         S_a=np.real(S_a), ks_b=ks_b, e_b=e_b, es_b=es_b, S_b=S_b,
         gap_ref=gap_ref)
