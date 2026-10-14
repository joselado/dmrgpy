"""Plot figure 12, mean-field_neel, from data.npz.

Writes docs/guide/img/mean-field_neel.png and -dark.png.
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
p = d["p"]
show = [0.0, 0.6, 0.8, 0.9, 1.0]


def draw(theme):
    c = style.colors(theme)
    fig, (a, b) = style.figure(2)
    sites = np.arange(n)
    k = 0
    for j, pj in enumerate(p):
        if not any(abs(pj - s) < 1e-9 for s in show):
            continue
        a.plot(sites, d["mz"][j], "o-", ms=3.5, color=c["series"][k],
               label=r"$p=%g$" % pj)
        k += 1
    a.set_xlabel(r"site $i$")
    a.set_ylabel(r"$\langle S^z_i\rangle$")
    a.set_ylim(-0.62, 0.75)
    a.legend(loc="upper center", ncol=5, fontsize=8, handlelength=1.2,
             columnspacing=0.8)

    b.plot(p, d["ms"], "o-", color=c["series"][0],
           label=r"staggered $m_s$ (left axis)")
    b.set_xlabel(r"fraction of exchange kept exactly, $p$")
    b.set_ylabel(r"staggered magnetization $m_s$", color=c["series"][0])
    b.set_ylim(-0.02, 0.55)
    b2 = b.twinx()
    b2.spines["right"].set_visible(True)
    b2.spines["top"].set_visible(False)
    b2.plot(p, d["e_var"]/n, "s", color=c["series"][1],
            label=r"energy $\langle H\rangle/n$ (right axis)")
    b2.axhline(-(n-1)/(4*n), **style.reference(theme))
    b2.axhline(float(d["e_exact"])/n, **style.reference(theme))
    b2.text(0.02, -(n-1)/(4*n), "Neel", va="bottom", fontsize=8,
            color=c["ink"])
    b2.text(0.02, float(d["e_exact"])/n, "exact chain", va="bottom",
            fontsize=8, color=c["ink"])
    b2.set_ylabel(r"energy per site [$J$]", color=c["series"][1])
    h1, l1 = b.get_legend_handles_labels()
    h2, l2 = b2.get_legend_handles_labels()
    b.legend(h1 + h2, l1 + l2, loc="lower left", bbox_to_anchor=(0.0, 0.08))
    return fig


style.main(draw, "mean-field_neel")
