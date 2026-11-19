"""Which site order keeps the bond dimension flat as baths are added?

A Heisenberg chain of ns spins 1/2, several of them Kondo coupled to their
own Wilson chain of n sites. Orders:
  linear : spins, then bath 0, bath 1, ... (every Kondo term long range)
  comb   : each contacted spin followed by its own bath
  ends   : two contacts only, bath 0 reversed, spins, bath 1
Usage: python order_probe.py c order maxm JK out.json [mode]
"""
import os, sys, json, time
sys.path.insert(0, os.environ.get("DMRGPY_SRC", "../../../../src"))
import numpy as np
from dmrgpy import fermionchain, wilsonchain
from dmrgpy.ysr import _hop

NS = 5
CONTACTS = {1: [2], 2: [1, 3], 3: [0, 2, 4], 4: [0, 1, 3, 4]}


def layout(ns, n, contacts, order):
    spin, bath = {}, {}
    if order == "linear":
        for i in range(ns):
            spin[i] = i
        for k, c in enumerate(contacts):
            bath[c] = [ns + k*n + m for m in range(n)]
    elif order == "comb":
        p = 0
        for i in range(ns):
            spin[i] = p
            p += 1
            if i in contacts:
                bath[i] = [p + m for m in range(n)]
                p += n
    elif order == "ends":
        assert len(contacts) == 2
        c0, c1 = contacts
        bath[c0] = [n - 1 - m for m in range(n)]
        for i in range(ns):
            spin[i] = n + i
        bath[c1] = [n + ns + m for m in range(n)]
    elif order == "sides":
        from dmrgpy.ysr import _metal_layout
        sp, bs = _metal_layout(ns, n, contacts, "sides")
        spin = {i: sp[i] for i in range(ns)}
        bath = {c: bs[k] for k, c in enumerate(contacts)}
    else:
        raise ValueError(order)
    return spin, bath


def build(ns, n, contacts, order, Js=0.5, JK=1.0, Lambda=2.0, D=1.0,
          maxm=64, mode="DMRG"):
    eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, band=(-D, D), z=0.0,
                                      mesh="mean")
    spin, bath = layout(ns, n, contacts, order)
    N = ns + n*len(contacts)
    fc = fermionchain.Spinful_Fermionic_Chain_Native(N)
    Upen = 2.0*(D + abs(JK) + abs(Js))
    S = lambda a, b: (fc.Sx[a]*fc.Sx[b] + fc.Sy[a]*fc.Sy[b]
                      + fc.Sz[a]*fc.Sz[b])
    h = 0
    for i in range(ns - 1):
        h = h + Js*S(spin[i], spin[i+1])
    for c in contacts:
        if JK != 0.0:
            h = h + JK*S(spin[c], bath[c][0])
    for i in range(ns):
        q = fc.Nup[spin[i]] + fc.Ndn[spin[i]] - 1
        h = h + Upen*q*q
    for c in contacts:
        b = bath[c]
        for m in range(n):
            if eps[m] != 0.0:
                h = h + eps[m]*(fc.Nup[b[m]] + fc.Ndn[b[m]])
            if m < n - 1:
                h = h + t[m]*_hop(fc, b[m], b[m+1])
    fc.set_hamiltonian(h)
    fc.maxm = maxm
    fc.nsweeps = 20
    if mode == "ED":
        fc.mode = "ED"
    return fc, spin, bath, (eps, t)


if __name__ == "__main__":
    c, order, maxm, JK, out = (int(sys.argv[1]), sys.argv[2],
                               int(sys.argv[3]), float(sys.argv[4]),
                               sys.argv[5])
    n = int(sys.argv[6]) if len(sys.argv) > 6 else 13
    t0 = time.time()
    fc, spin, bath, _ = build(NS, n, CONTACTS[c], order, JK=JK, maxm=maxm)
    e = float(np.real(fc.gs_energy()))
    wf = fc.get_gs()
    S = [float(wf.get_bond_entropy(b)) for b in range(fc.ns - 1)]
    res = dict(c=c, order=order, maxm=maxm, JK=JK, n=n, e=e, S=S,
               spin=spin, bath={str(k): v for k, v in bath.items()},
               time=time.time() - t0)
    json.dump(res, open(out, "w"))
    print(c, order, maxm, JK, "%.12f" % e, "Smax=%.3f" % max(S),
          "t=%.0fs" % res["time"])
