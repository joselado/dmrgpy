"""Plot figure 6, excited-states_heff, from data.npz.

Writes docs/guide/img/excited-states_heff.png and -dark.png.
    python3 plot.py [--theme light|dark]
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import common
import matplotlib.ticker as ticker
import style

params, d = common.load(os.path.join(HERE, "data.npz"))
U = d["U"]


def draw(theme):
    c = style.colors(theme)
    fig, ax = style.figure(1)
    ax.plot(U, d["jfit"][:, 2], "o", color=c["series"][0],
            label=r"get_heff, fitted $J$")
    ax.plot(U, d["j_proj"], **style.reference(theme),
            label=r"exact, $(E_T-E_S)/(1/4+3w/4)$")
    ax.plot(U, d["j_exact"], color=c["ink"], ls=":", lw=style.LW_REF,
            label=r"splitting $E_T-E_S$")
    ax.plot(U, d["j_strong"], color=c["ink2"], ls="-.", lw=style.LW_REF,
            label=r"$4t^2/U$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([1, 2, 5, 10, 20])
    ax.set_yticks([0.2, 0.5, 1, 2, 4])
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        axis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel(r"$U$ [$t$]")
    ax.set_ylabel(r"exchange $J$ [$t$]")
    ax.legend()
    return fig


style.main(draw, "excited-states_heff")
