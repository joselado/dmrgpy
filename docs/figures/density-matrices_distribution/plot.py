"""Plot figure 9, density-matrices_distribution, from data.npz.

Writes docs/guide/img/density-matrices_distribution.png and -dark.png.
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
h = d["h"]
ns, nl = params["n_small"], params["n_large"]
names = {0.2: "ordered", 0.5: "critical", 1.0: "paramagnet"}
# bulk magnetization per site of the ordered phase (Pfeuty), m = (1-(2h)^2)^(1/8)/2
m_pfeuty = 0.5*(1-(2*h[0])**2)**0.125


def draw(theme):
    c = style.colors(theme)
    fig, (a, b) = style.figure(2)
    for k, hk in enumerate(h):
        lab = r"$h=%g$, %s" % (hk, names.get(float(hk), ""))
        a.plot(d["m_small"], d["p_small"][k], color=c["series"][k], label=lab)
        b.plot(d["m_large"], d["p_large"][k], color=c["series"][k], label=lab)
    for k in range(len(h)):
        a.plot(d["m_small"], d["p_small_ed"][k], **style.reference(theme),
               label="ED" if k == 0 else None)
    b.axvline(m_pfeuty, color=c["ink2"], lw=0.8, ls=":")
    b.text(m_pfeuty, b.get_ylim()[1]*0.97, " bulk $m$\n (exact)",
           color=c["ink2"], ha="right", va="top", fontsize=8)
    for ax, n in ((a, ns), (b, nl)):
        ax.set_xlabel(r"magnetization per site $m$")
        ax.set_xlim(-0.55, 0.55)
        ax.set_ylim(bottom=0)
        ax.text(0.02, 0.97, r"$n=%d$" % n, transform=ax.transAxes,
                ha="left", va="top")
    a.set_ylabel(r"$P(m)$")
    a.legend(loc="upper left", bbox_to_anchor=(0.0, 0.9))
    return fig


style.main(draw, "density-matrices_distribution")
