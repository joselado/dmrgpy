"""Plot figure 2, ground-state_correlators, from data.npz.

Writes docs/guide/img/ground-state_correlators.png and -dark.png.
    python3 plot.py [--theme light|dark]
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import common
import matplotlib.ticker as ticker
import style

params, d = common.load(os.path.join(HERE, "data.npz"))
r = d["r"]
ok = r >= 2          # sqrt(ln r)/r vanishes at r=1


def draw(theme):
    c = style.colors(theme)
    fig, ax = style.figure(1)
    ax.plot(r, np.abs(d["czz_half"]), "o-", color=c["series"][0], ms=4,
            label=r"$S=1/2$")
    ax.plot(r, np.abs(d["czz_one"]), "s-", color=c["series"][1], ms=4,
            label=r"$S=1$")
    ax.plot(r[ok], d["ref_half"][ok], **style.reference(theme),
            label=r"$a\sqrt{\ln r}/r$")
    ax.plot(r, d["ref_one"], color=c["ink"], ls=":", lw=style.LW_REF,
            label=r"$b\,e^{-r/\xi}/\sqrt{r}$, $\xi=%.2f$" % params["xi_haldane"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([1, 2, 4, 8, 16, 32])
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel(r"distance $r$ [sites]")
    ax.set_ylabel(r"$|\langle S^z_iS^z_{i+r}\rangle|$")
    ax.legend()
    return fig


style.main(draw, "ground-state_correlators")
