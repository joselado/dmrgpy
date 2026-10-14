"""Figure 18, kondo_dIdV: STM dI/dV spectra of magnetic atoms and chains,
in the third-order perturbation theory of Ternes, NJP 17, 063016 (2015).

(a) second-order inelastic steps of an S=2 atom with anisotropy
    D Sz^2 + E (Sx^2 - Sy^2), D=-1.55 meV, E=0.31 meV (Fe on Cu2N), B=0;
    reference: the excitation energies of the exact 5x5 Hamiltonian.
(b) third-order Kondo peak of an S=1/2 in a field B (paper's Fig. 3b);
    reference: the paper's B=0 values, 1.13 at zero bias and 0.886 at
    +-4 mV, in units of e^2 T0^2/h.
(c) third-order T=0 spectrum of a 3-site chain, mode="DMRG" (TDVP two-time
    construction) against mode="ED", as in
    examples/kondo/kondo_third_order_dmrg_VS_ED.

Run with  bash docs/figures/run.sh kondo_dIdV

Arrays in data.npz (energies in eV, dI/dV in the units of
get_kondo_spectrum unless noted):
    eV_a       (401,)    bias of panel a
    dIdV_a     (401,)    second-order dI/dV of the S=2 atom
    levels_a   (5,)      exact eigenvalues of the 5x5 Hamiltonian, minus E0
    eps_a      (k,)      excitation energies from the ground state that carry
                         second-order weight (the step positions)
    weight_a   (k,)      their weights |<f|S-|0>|^2/2+|<f|S+|0>|^2/2+|<f|Sz|0>|^2
    B_b        (6,)      fields of panel b, in tesla
    eV_b       (161,)    bias of panel b
    dIdV_b     (6,161)   third-order dI/dV divided by 2 pi, in e^2 T0^2/h
    eV_c       (41,)     bias of panel c
    dIdV_ed    (41,)     mode="ED"
    dIdV_dmrg  (41,)     mode="DMRG"
    params               JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain

MUB = 5.7883818066e-5   # eV/T
params = common.start(dict(
    units="H in eV, T in kelvin",
    a=dict(spin="S=2", D=-1.55e-3, E=0.31e-3, B=0.0, order=2, T=0.5,
           eV="np.linspace(-10e-3, 10e-3, 401)"),
    b=dict(spin="S=1/2", g=2.0, B=[0, 2, 4, 6, 8, 10], Jrho_s=-0.05, U=0.0,
           T=1.0, order=3, omega0=20e-3, eV="np.linspace(-4e-3, 4e-3, 161)"),
    c=dict(chain="3-site S=1/2, g muB B Sz_0 with B=10 T, "
                 "plus 0.01 eV Heisenberg couplings",
           Jrho_s=0.05, U=0.0, T=0.0, order=3, omega0=2e-3, Gamma0=2e-3,
           maxm=30, nsweeps=15, submode="KPM", delta=2e-5,
           es="np.linspace(-3e-3, 3e-3, 800)", dt2=150.0, n_t2_half=20,
           dtau=150.0, n_tau_half=20, itensor_version=3,
           eV="np.linspace(-2e-3, 2e-3, 41)"),
))

# (a) S=2 atom, second order
pa = params["a"]
sc = spinchain.Spin_Chain([pa["spin"]])
s = sc
sc.set_hamiltonian(pa["D"]*s.Sz[0]*s.Sz[0]
                   + pa["E"]*(s.Sx[0]*s.Sx[0] - s.Sy[0]*s.Sy[0]))
eV_a, dIdV_a = sc.get_kondo_spectrum(np.linspace(-10e-3, 10e-3, 401), site=0,
                                     order=pa["order"], T=pa["T"])
# exact 5x5 multiplet
m = np.arange(2, -3, -1).astype(float)
sz = np.diag(m)
sp = np.zeros((5, 5))
for k in range(1, 5):
    sp[k-1, k] = np.sqrt(6 - m[k]*(m[k] + 1))
sx, sy = (sp + sp.T)/2, (sp - sp.T)/2j
hm = (pa["D"]*sz @ sz + pa["E"]*(sx @ sx - sy @ sy)).real
w, v = np.linalg.eigh(hm)
eps_a, weight_a = [], []
for f in range(1, 5):
    wt = (0.5*abs(v[:, f] @ sp @ v[:, 0])**2 + 0.5*abs(v[:, f] @ sp.T @ v[:, 0])**2
          + abs(v[:, f] @ sz @ v[:, 0])**2)
    if wt > 1e-8:
        eps_a.append(w[f] - w[0])
        weight_a.append(wt)
print("(a) levels [meV]:", np.round((w - w[0])*1e3, 4))
print("(a) steps from the ground state [meV]:", np.round(np.array(eps_a)*1e3, 4),
      "weights", np.round(weight_a, 3), flush=True)

# (b) S=1/2 in a field, third order
pb = params["b"]
eV_b = np.linspace(-4e-3, 4e-3, 161)
dIdV_b = []
for B in pb["B"]:
    sc = spinchain.Spin_Chain([pb["spin"]])
    sc.set_hamiltonian(pb["g"]*MUB*B*sc.Sz[0])
    _, y = sc.get_kondo_spectrum(eV_b, site=0, Jrho_s=pb["Jrho_s"], U=pb["U"],
                                 T=pb["T"], order=pb["order"],
                                 omega0=pb["omega0"])
    y = y/(2*np.pi)
    dIdV_b.append(y)
    print("(b) B=%4.1f T  zero bias %.4f  tails %.4f %.4f  maximum at %+.3f meV"
          % (B, y[80], y[0], y[-1], eV_b[np.argmax(y)]*1e3), flush=True)
print("(b) paper at B=0: zero bias 1.13, tails 0.886", flush=True)

# (c) 3-site chain, DMRG against ED
pc = params["c"]
sc = spinchain.Spin_Chain(["S=1/2"]*3, itensor_version=pc["itensor_version"])
sc.maxm, sc.nsweeps = pc["maxm"], pc["nsweeps"]
h = 2.0*MUB*10.0*sc.Sz[0]
for i in range(2):
    h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)
sc.get_gs()
eV_c = np.linspace(-2e-3, 2e-3, 41)
t0 = time.time()
_, dIdV_ed = sc.get_kondo_spectrum(eV_c, site=0, Jrho_s=pc["Jrho_s"], U=pc["U"],
                                   T=pc["T"], order=3, mode="ED",
                                   omega0=pc["omega0"], Gamma0=pc["Gamma0"])
print("(c) ED done (%.1f s)" % (time.time() - t0), flush=True)
t0 = time.time()
_, dIdV_dmrg = sc.get_kondo_spectrum(
    eV_c, site=0, Jrho_s=pc["Jrho_s"], U=pc["U"], T=pc["T"], order=3,
    mode="DMRG", omega0=pc["omega0"], Gamma0=pc["Gamma0"],
    submode=pc["submode"], delta=pc["delta"], es=np.linspace(-3e-3, 3e-3, 800),
    dt2=pc["dt2"], n_t2_half=pc["n_t2_half"], dtau=pc["dtau"],
    n_tau_half=pc["n_tau_half"])
diff = np.max(np.abs(dIdV_dmrg - dIdV_ed))
print("(c) DMRG done (%.1f s); max |DMRG-ED| = %.4f on max |ED| = %.4f (%.2f%%)"
      % (time.time() - t0, diff, np.max(np.abs(dIdV_ed)),
         100*diff/np.max(np.abs(dIdV_ed))), flush=True)

common.save(__file__, params, eV_a=eV_a, dIdV_a=dIdV_a, levels_a=w - w[0],
            eps_a=np.array(eps_a), weight_a=np.array(weight_a),
            B_b=np.array(pb["B"], dtype=float), eV_b=eV_b,
            dIdV_b=np.array(dIdV_b), eV_c=eV_c, dIdV_ed=dIdV_ed,
            dIdV_dmrg=dIdV_dmrg)
