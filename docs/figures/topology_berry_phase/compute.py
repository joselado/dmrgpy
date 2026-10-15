# (a) Many-body Berry phase of a dimerized spin-1/2 Heisenberg ring with a
# twist exp(i phi) on the closing bond, against the Zak phase of the SSH chain.
# (b) Fermion parity and total S^z of a superconducting chain with a magnetic
# impurity (a Zeeman field B on the first site). ED is the reference for both.
import numpy as np
from dmrgpy import spinchain, fermionchain, topology

n, nphi = 12, 16
deltas = np.linspace(-0.6, 0.6, 12).round(8)  # J_i = 1 + delta (-1)^i


def ring_phase(sc, delta, mode):
    """Berry phase/pi, the phase of prod_phi <gs(phi)|gs(phi+dphi)>."""
    J = 1 + delta*(-1.0)**np.arange(n)
    hbulk = sum(J[i]*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                      + sc.Sz[i]*sc.Sz[i+1]) for i in range(n-1))
    Spn, Sp0 = sc.Sx[n-1] + 1j*sc.Sy[n-1], sc.Sx[0] + 1j*sc.Sy[0]
    wfs = []
    for phi in np.linspace(0, 2*np.pi, nphi, endpoint=False):
        z = np.exp(1j*phi)
        sc.set_hamiltonian(hbulk + J[n-1]*(z*Spn*Sp0.get_dagger()/2
                                           + np.conj(z)*Sp0*Spn.get_dagger()/2
                                           + sc.Sz[0]*sc.Sz[n-1]))
        wfs.append(sc.get_gs(mode=mode))
    prod = np.prod([wfs[i].dot(wfs[(i+1) % nphi]) for i in range(nphi)])
    return (np.angle(prod)/np.pi) % 2


gamma, gamma_ed = [], []
for delta in deltas:
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    sc.maxm, sc.nsweeps = 64, 12
    gamma.append(ring_phase(sc, delta, "DMRG"))
    gamma_ed.append(ring_phase(sc, delta, "ED"))


def ssh(delta):  # Bloch Hamiltonian, k in units of 2 pi
    def hk(k):
        dx = (1+delta) + (1-delta)*np.cos(2*np.pi*k)
        dy = (1-delta)*np.sin(2*np.pi*k)
        return np.array([[0, dx-1j*dy], [dx+1j*dy, 0]])
    return hk


delta_ssh = np.linspace(-0.6, 0.6, 121).round(8)
zak = [(topology.berry_phase_matrix(ssh(d), nk=60)/np.pi) % 2
       for d in delta_ssh]


def impurity(B, L=6, mu=-0.5, D=0.5, U=0.5):
    fc = fermionchain.Spinful_Fermionic_Chain(L)
    h = 0
    for i in range(L-1):  # hopping t=1
        h = h - (fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagup[i+1]*fc.Cup[i])
        h = h - (fc.Cdagdn[i]*fc.Cdn[i+1] + fc.Cdagdn[i+1]*fc.Cdn[i])
    for i in range(L):
        h = h - mu*(fc.Nup[i] + fc.Ndn[i])
        h = h + D*(fc.Cup[i]*fc.Cdn[i] + fc.Cdagdn[i]*fc.Cdagup[i])  # pairing
        h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    fc.set_hamiltonian(h + B*fc.Sz[0])
    fc.maxm, fc.nsweeps = 60, 30
    return fc


bz = np.linspace(0, 4, 41).round(8)
parity, parity_ed, sz, sz_ed = [], [], [], []
for B in bz:
    fc = impurity(B)
    parity.append(fc.get_gs().get_fermionic_parity().real)
    sz.append(fc.vev(sum(fc.Sz)).real)
    parity_ed.append(fc.get_gs(mode="ED").get_fermionic_parity().real)
    sz_ed.append(fc.vev(sum(fc.Sz), mode="ED").real)

# field where the ED parity flips, by bisection
lo, hi = 0.0, 4.0
while hi - lo > 1e-5:
    mid = (lo + hi)/2
    if impurity(mid).get_gs(mode="ED").get_fermionic_parity().real > 0:
        lo = mid
    else:
        hi = mid

np.savez("data.npz", delta=deltas, gamma=gamma, gamma_ed=gamma_ed,
         delta_ssh=delta_ssh, zak=zak, bz=bz, parity=parity,
         parity_ed=parity_ed, sz=sz, sz_ed=sz_ed, b_c=(lo + hi)/2)
