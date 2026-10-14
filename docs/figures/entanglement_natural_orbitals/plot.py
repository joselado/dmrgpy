"""Plot figure 8, entanglement_natural_orbitals, from data.npz.

Writes docs/guide/img/entanglement_natural_orbitals.png and -dark.png.
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
n = params["n"]
U = d["U"]


def draw(theme):
    c = style.colors(theme)
    fig, (a, b) = style.figure(2)
    alpha = np.arange(1, 2*n+1)
    for k, u in enumerate(U):
        a.plot(alpha, d["occ"][k], "o-", ms=3.5, color=c["series"][k],
               label=r"$U=%g$" % u)
    a.axhline(0.5, color=c["ink2"], lw=0.8, ls=":")
    a.set_xlabel(r"natural orbital $\alpha$")
    a.set_ylabel(r"occupation $n_\alpha$")
    a.set_ylim(-0.03, 1.03)
    a.legend(loc="lower left", ncol=2)

    b.plot(U, d["s"], "o-", color=c["series"][0], label=r"DMRG, $n=%d$" % n)
    b.plot(d["U_fine"], d["s_fine_ed"], **style.reference(theme),
           label=r"ED, $n=%d$" % params["n_small"])
    b.plot(U, d["s_small_dmrg"], "s", mfc="none", color=c["series"][1],
           mew=1.2, label=r"DMRG, $n=%d$" % params["n_small"])
    b.axhline(2*n*np.log(2), color=c["ink2"], lw=0.8, ls=":")
    b.text(U[-1], 2*n*np.log(2), r"$2n\ln 2$", color=c["ink2"],
           ha="right", va="bottom", fontsize=8)
    b.set_xlabel(r"$U$ [$t$]")
    b.set_ylabel(r"correlation entropy $S$")
    b.legend(loc="center right")
    return fig


style.main(draw, "entanglement_natural_orbitals")
