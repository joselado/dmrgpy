"""Figure 7, entanglement_central_charge: bond entropies S(l) of open
critical chains (XX, Heisenberg, critical Ising) against the gapped
spin-1 Heisenberg chain, and the central charge from Calabrese-Cardy.

Run with  bash docs/figures/run.sh entanglement_central_charge

Arrays in data.npz:
    names     (4,)    chain labels: xx, heisenberg, ising, spin1
    ell       (L-1,)  size of the left block, l = 1 .. L-1
    entropy   (4,L-1) DMRG bond entropy S(l) of each chain
    s_exact   (L-1,)  exact S(l) of the XX chain from free fermions
    c_lib     (4,)    wf.get_CFT_central_charge() of each chain
    c_even    (4,)    central charge fitted on even l only
    energy    (4,)    DMRG ground-state energy of each chain
    params            JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain

L = int(os.environ.get("FIG_L", 100))
params = common.start(dict(
    models={"xx": "sum SxSx+SySy, c=1",
            "heisenberg": "sum S.S, c=1",
            "ising": "sum SzSz + 1/2 sum Sx (critical), c=1/2",
            "spin1": "spin-1 Heisenberg, gapped"},
    boundary="open",
    L=L,
    maxm={"xx": 200, "heisenberg": 200, "ising": 120, "spin1": 100},
    nsweeps=20,
    cutoff=1e-12,
    itensor_version=3,
    fit="S = c/6 x + const on even l, x = ln[2L/pi sin(pi l/L)], "
        "l from 4 to L-4",
    seeds="none (ITensor v3 random initial MPS)",
))
names = ["xx", "heisenberg", "ising", "spin1"]
ell = np.arange(1, L)
x = np.log(2*L/np.pi*np.sin(np.pi*ell/L))


def chain(name):
    spin = "S=1" if name == "spin1" else "S=1/2"
    sc = spinchain.Spin_Chain([spin]*L,
                              itensor_version=params["itensor_version"])
    if name == "xx":
        h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] for i in range(L-1))
    elif name == "ising":
        h = sum(sc.Sz[i]*sc.Sz[i+1] for i in range(L-1))
        h = h + 0.5*sum(sc.Sx[i] for i in range(L))
    else:
        h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                + sc.Sz[i]*sc.Sz[i+1] for i in range(L-1))
    sc.set_hamiltonian(h)
    sc.maxm = params["maxm"][name]
    sc.nsweeps = params["nsweeps"]
    sc.cutoff = params["cutoff"]
    return sc


def fit_even(s):
    sel = (ell % 2 == 0) & (ell >= 4) & (ell <= L-4)
    slope, const = np.polyfit(x[sel], s[sel], 1)
    return 6*slope


# exact XX chain: free fermions with hopping 1/2 (Jordan-Wigner), half filling
hop = np.zeros((L, L))
for i in range(L-1):
    hop[i, i+1] = hop[i+1, i] = 0.5
ev, vec = np.linalg.eigh(hop)
occ = vec[:, :L//2]
cmat = occ @ occ.T
s_exact = []
for l in ell:
    nu = np.clip(np.linalg.eigvalsh(cmat[:l, :l]), 1e-15, 1-1e-15)
    s_exact.append(-np.sum(nu*np.log(nu) + (1-nu)*np.log(1-nu)))
s_exact = np.array(s_exact)
print("exact XX: E0 = %.10f, c(even fit) = %.4f"
      % (2*np.sum(ev[:L//2]) * 0.5, fit_even(s_exact)), flush=True)

entropy, c_lib, c_even, energy = [], [], [], []
for name in names:
    t0 = time.time()
    sc = chain(name)
    e = sc.gs_energy()
    wf = sc.get_gs()
    s = np.array([wf.get_bond_entropy(i, i+1) for i in range(L-1)])
    cl = wf.get_CFT_central_charge()
    ce = fit_even(s)
    entropy.append(s)
    c_lib.append(cl)
    c_even.append(ce)
    energy.append(e)
    print("%-10s E0=%.10f  c(lib)=%.4f  c(even)=%.4f  S(L/2)=%.6f  (%.1f s)"
          % (name, e, cl, ce, s[L//2-1], time.time()-t0), flush=True)
    if name == "xx":
        print("  XX: E0-E_exact = %.3e, max|S-S_exact| = %.3e"
              % (e - np.sum(ev[:L//2]), np.max(np.abs(s - s_exact))),
              flush=True)

common.save(__file__, params, names=np.array(names), ell=ell,
            entropy=np.array(entropy), s_exact=s_exact,
            c_lib=np.array(c_lib), c_even=np.array(c_even),
            energy=np.array(energy))
