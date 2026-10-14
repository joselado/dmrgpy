"""Figure 10, topology_berry_phase: three invariants that jump.

(a) The many-body Berry phase of a dimerized spin-1/2 Heisenberg ring,
    J_i = 1 + delta (-1)^i, with a twist exp(i phi) on the closing bond.
    Reference: the Zak phase of the SSH chain with the same dimerization,
    from topology.berry_phase_matrix, and ED on the same ring.
(b) The fermion parity and the total S^z of a short spinful chain with
    onsite pairing and a Zeeman field B on the first site (a magnetic
    impurity in a superconductor). Reference: ED on the same chain.

Run with  bash docs/figures/run.sh topology_berry_phase

Arrays in data.npz:
    delta      (k,)  dimerizations
    gamma      (k,)  Berry phase of the spin ring, DMRG, in units of pi, in [0,2)
    gamma_ed   (k,)  the same from ED
    delta_ssh  (q,)  a fine dimerization grid for the SSH reference
    zak        (q,)  SSH Zak phase from berry_phase_matrix, units of pi, in [0,2)
    bz         (m,)  Zeeman field on the impurity site
    parity     (m,)  fermion parity of the ground state, DMRG
    parity_ed  (m,)  the same from ED
    sz         (m,)  total <S^z> of the ground state, DMRG
    sz_ed      (m,)  the same from ED
    ntot       (m,)  total <N>, DMRG
    b_c        ()    field where the ED parity flips, by bisection
    params           JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain, fermionchain, topology

params = common.start(dict(
    ring=dict(model="spin-1/2 Heisenberg ring, J_i = 1 + delta (-1)^i, "
                    "twist exp(i phi) on the bond (n-1, 0)",
              n=12,
              nphi=16,
              delta=list(np.round(np.linspace(-0.6, 0.6, 12), 8)),
              maxm=64, nsweeps=12),
    ssh=dict(model="h(k) = [(1+delta)+(1-delta)cos k] sx + (1-delta) sin k sy",
             nk=60, delta=list(np.round(np.linspace(-0.6, 0.6, 121), 8))),
    impurity=dict(model="H = -sum_s (c+_is c_i+1s + h.c.) - mu N "
                        "+ D sum (c_iup c_idn + h.c.) "
                        "+ U sum (Nup_i-1/2)(Ndn_i-1/2) + B Sz_0",
                  n=6, mu=-0.5, D=0.5, U=0.5,
                  B=list(np.round(np.linspace(0, 4, 41), 8)),
                  maxm=60, nsweeps=30),
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
pr, ps, pi_ = params["ring"], params["ssh"], params["impurity"]


# (a) spin ring
def ring_phase(sc, delta, mode):
    n = pr["n"]
    hbulk = 0
    for i in range(n-1):
        J = 1. + delta*(-1)**i
        hbulk = hbulk + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                           + sc.Sz[i]*sc.Sz[i+1])
    J = 1. + delta*(-1)**(n-1)
    Spn = sc.Sx[n-1] + 1j*sc.Sy[n-1]
    Sp0 = sc.Sx[0] + 1j*sc.Sy[0]
    wfs = []
    for phi in np.linspace(0, 2*np.pi, pr["nphi"], endpoint=False):
        z = np.exp(1j*phi)
        h = hbulk + J*(z*Spn*Sp0.get_dagger()/2 + np.conj(z)*Sp0*Spn.get_dagger()/2
                       + sc.Sz[0]*sc.Sz[n-1])
        sc.set_hamiltonian(h)
        wfs.append(sc.get_gs(mode=mode))
    prod = 1.0
    for i in range(len(wfs)):
        prod = prod*wfs[i].dot(wfs[(i+1) % len(wfs)])
    return (np.angle(prod)/np.pi) % 2, abs(prod)


gamma, gamma_ed = [], []
for d in pr["delta"]:
    t0 = time.time()
    sc = spinchain.Spin_Chain(["S=1/2"]*pr["n"],
                              itensor_version=params["itensor_version"])
    sc.maxm, sc.nsweeps = pr["maxm"], pr["nsweeps"]
    g, a = ring_phase(sc, d, "DMRG")
    ge, ae = ring_phase(sc, d, "ED")
    gamma.append(g)
    gamma_ed.append(ge)
    print("delta=%+.4f  gamma/pi DMRG=%.6f |prod|=%.4f  ED=%.6f |prod|=%.4f"
          "  (%.1f s)" % (d, g, a, ge, ae, time.time()-t0), flush=True)


# SSH Zak phase, k in units of 2 pi as berry_phase_matrix expects
def ssh(delta):
    def hk(k):
        k = 2*np.pi*k
        dx = (1+delta) + (1-delta)*np.cos(k)
        dy = (1-delta)*np.sin(k)
        return np.array([[0, dx-1j*dy], [dx+1j*dy, 0]])
    return hk


zak = np.array([(topology.berry_phase_matrix(ssh(d), nk=ps["nk"])/np.pi) % 2
                for d in ps["delta"]])
print("SSH Zak phase/pi at delta=-0.3, +0.3:",
      zak[np.argmin(np.abs(np.array(ps["delta"])+0.3))],
      zak[np.argmin(np.abs(np.array(ps["delta"])-0.3))], flush=True)


# (b) magnetic impurity in a superconducting chain
def impurity(B):
    L = pi_["n"]
    fc = fermionchain.Spinful_Fermionic_Chain(
        L, itensor_version=params["itensor_version"])
    h = 0
    for i in range(L-1):
        h = h - (fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagup[i+1]*fc.Cup[i])
        h = h - (fc.Cdagdn[i]*fc.Cdn[i+1] + fc.Cdagdn[i+1]*fc.Cdn[i])
    for i in range(L):
        h = h - pi_["mu"]*(fc.Nup[i] + fc.Ndn[i])
        h = h + pi_["D"]*(fc.Cup[i]*fc.Cdn[i] + fc.Cdagdn[i]*fc.Cdagup[i])
        h = h + pi_["U"]*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    h = h + B*fc.Sz[0]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = pi_["maxm"], pi_["nsweeps"]
    return fc


parity, parity_ed, sz, sz_ed, ntot = [], [], [], [], []
for B in pi_["B"]:
    fc = impurity(B)
    wf = fc.get_gs()
    Sz = sum(fc.Sz)
    parity.append(wf.get_fermionic_parity().real)
    sz.append(fc.vev(Sz).real)
    ntot.append(fc.vev(sum(fc.N)).real)
    wfe = fc.get_gs(mode="ED")
    parity_ed.append(wfe.get_fermionic_parity().real)
    sz_ed.append(fc.vev(Sz, mode="ED").real)
    print("B=%.2f  P=%+.6f (ED %+.6f)  Sz=%+.6f (ED %+.6f)  N=%.5f"
          % (B, parity[-1], parity_ed[-1], sz[-1], sz_ed[-1], ntot[-1]),
          flush=True)

lo, hi = 0.0, 4.0
while hi - lo > 1e-5:
    mid = 0.5*(lo + hi)
    if impurity(mid).get_gs(mode="ED").get_fermionic_parity().real > 0:
        lo = mid
    else:
        hi = mid
b_c = 0.5*(lo + hi)
print("ED parity flips at B_c = %.5f" % b_c, flush=True)

common.save(__file__, params, delta=np.array(pr["delta"]),
            gamma=np.array(gamma), gamma_ed=np.array(gamma_ed),
            delta_ssh=np.array(ps["delta"]), zak=zak,
            bz=np.array(pi_["B"]), parity=np.array(parity),
            parity_ed=np.array(parity_ed), sz=np.array(sz),
            sz_ed=np.array(sz_ed), ntot=np.array(ntot),
            b_c=np.array(b_c))
