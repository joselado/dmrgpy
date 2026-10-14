"""Plot figure 18, kondo_dIdV, from data.npz.

Writes docs/guide/img/kondo_dIdV.png and -dark.png.
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
mev = 1e3
eb, yb = d["eV_b"], d["dIdV_b"]
print("(b) B=0: zero bias %.4f (paper 1.13), tails %.4f %.4f (paper 0.886)"
      % (yb[0][len(eb)//2], yb[0][0], yb[0][-1]))
diff = np.max(np.abs(d["dIdV_dmrg"] - d["dIdV_ed"]))
print("(c) max |DMRG-ED| = %.4f, %.2f%% of max |ED|"
      % (diff, 100*diff/np.max(np.abs(d["dIdV_ed"]))))


def draw(theme):
    c = style.colors(theme)
    fig, (ax1, ax2, ax3) = style.figure(3)
    for e in d["eps_a"]:
        for s in (+1, -1):
            ax1.axvline(s*e*mev, color=c["ink2"], ls=":", lw=style.LW_REF)
    ax1.plot(d["eV_a"]*mev, d["dIdV_a"], color=c["series"][0], lw=style.LW)
    ax1.set_xlabel(r"$eV$ [meV]")
    for ax, t in zip((ax1, ax2, ax3), "abc"):
        ax.text(0.0, 1.02, "(%s)" % t, transform=ax.transAxes, va="bottom",
                color=c["ink"])
    ax1.set_ylabel(r"$dI/dV$ [arb. units]")

    for j, B in enumerate(d["B_b"]):
        ax2.plot(eb*mev, yb[j], color=c["series"][j], lw=style.LW,
                 label=r"$B=%g$ T" % B)
    for v in (1.13, 0.886):
        ax2.axhline(v, color=c["ink"], ls=":", lw=style.LW_REF)
    ax2.set_xlabel(r"$eV$ [meV]")
    ax2.set_ylabel(r"$dI/dV$ [$e^2T_0^2/h$]")
    ax2.legend(ncol=3, loc="upper center", handlelength=1.2, columnspacing=0.8)
    ax2.set_ylim(0.28, 1.42)

    ax3.plot(d["eV_c"]*mev, d["dIdV_ed"], **style.reference(theme), label="ED")
    ax3.plot(d["eV_c"]*mev, d["dIdV_dmrg"], "o", color=c["series"][0],
             ms=style.MS, label="DMRG")
    ax3.set_xlabel(r"$eV$ [meV]")
    ax3.set_ylabel(r"$dI/dV$ [arb. units]")
    ax3.legend()
    return fig


style.main(draw, "kondo_dIdV")
