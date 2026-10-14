"""Plot figure 3, ground-state_sectors, from data.npz.

Writes docs/guide/img/ground-state_sectors.png and -dark.png.
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
M = d["M"]
HMAX = 2.4


def staircase(e):
    """m(h) as a step curve: M for h between h_{M-1} and h_M."""
    hm = np.diff(e)                       # h_M = E0(M+1) - E0(M)
    x = np.concatenate([[0.0], hm, [HMAX]])
    y = np.concatenate([M, [M[-1]]])/n
    return x, y


def draw(theme):
    c = style.colors(theme)
    fig, ax = style.figure(1)
    x, y = staircase(d["e_xx"])
    ax.step(x, y, where="post", color=c["series"][0], label="XX, DMRG")
    x, y = staircase(d["e_heis"])
    ax.step(x, y, where="post", color=c["series"][1], label="Heisenberg, DMRG")
    h = np.linspace(0, 1, 400)
    hh = np.concatenate([h, [HMAX]])
    ax.plot(hh, np.concatenate([np.arcsin(h)/np.pi, [0.5]]),
            **style.reference(theme), label=r"XX, exact $\arcsin(h/J)/\pi$")
    ax.axvline(2.0, color=c["ink"], ls=":", lw=style.LW_REF)
    ax.text(2.03, 0.05, r"$h_s=2J$", color=c["ink"])
    ax.set_xlim(0, HMAX)
    ax.set_ylim(0, 0.53)
    ax.set_xlabel(r"field $h$ [$J$]")
    ax.set_ylabel(r"magnetization per site $M/n$")
    ax.legend(loc="upper left")
    return fig


style.main(draw, "ground-state_sectors")
