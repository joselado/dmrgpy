"""Plot figure 7, entanglement_central_charge, from data.npz.

Writes docs/guide/img/entanglement_central_charge.png and -dark.png.
    python3 plot.py [--theme light|dark]
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import common
import style

params, d = common.load(os.path.join(HERE, "data.npz"))
L = params["L"]
ell = d["ell"]
S = d["entropy"]
x = np.log(2*L/np.pi*np.sin(np.pi*ell/L))
even = (ell % 2 == 0) & (ell >= 4) & (ell <= L-4)
labels = {"xx": "XX", "heisenberg": "Heisenberg", "ising": "critical Ising",
          "spin1": "spin-1 Heisenberg"}
names = [str(a) for a in d["names"]]


def draw(theme):
    c = style.colors(theme)
    fig, (a, b) = style.figure(2)
    for k, name in enumerate(names):
        a.plot(ell, S[k], "-", color=c["series"][k], label=labels[name])
    a.plot(ell[::3], d["s_exact"][::3], "o", mfc="none", mec=c["ink"],
           ms=4, mew=0.9, label="XX, free fermions")
    a.set_xlabel(r"block length $\ell$")
    a.set_ylabel(r"entanglement entropy $S(\ell)$")
    a.set_xlim(0, L)
    a.set_ylim(bottom=0)
    a.legend(loc="lower center", ncol=2)

    xs = np.linspace(x[even].min(), x[even].max(), 2)
    for k, name in enumerate(names):
        if name == "spin1":
            lab = "%s (gapped)" % labels[name]
        else:
            lab = r"%s, fit $c=%.2f$" % (labels[name], d["c_even"][k])
        b.plot(x[even], S[k][even], "o", color=c["series"][k], ms=3.5,
               label=lab)
    for name, slope in (("heisenberg", 1/6), ("ising", 1/12)):
        k = names.index(name)
        s0 = np.mean(S[k][even] - slope*x[even])
        b.plot(xs, slope*xs + s0, **style.reference(theme),
               label=r"slope $c/6$ with $c=1$ and $c=1/2$"
               if name == "heisenberg" else None)
    b.set_xlabel(r"$\ln\left[\frac{2L}{\pi}\sin\frac{\pi\ell}{L}\right]$")
    b.set_ylabel(r"$S(\ell)$, even $\ell$")
    b.set_ylim(0.12, 1.1)
    b.legend(loc="lower right", fontsize=8)
    return fig


style.main(draw, "entanglement_central_charge")
