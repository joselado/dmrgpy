"""Plot figure 22, post-processing_extrapolation, from data.npz.

Writes docs/guide/img/post-processing_extrapolation.png and -dark.png.
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
L = d["L"]
x = 1.0/L
y = d["energy"]/L
e_exact = float(d["e_exact"])
e1, e2 = float(d["e_inf_n1"]), float(d["e_inf_n2"])
# the same weighted quadratic fit as extrapolate.size_extrapolation(n=2)
c = np.polynomial.polynomial.polyfit(x, y, 2, w=L)
assert abs(c[0] - e2) < 1e-12
xf = np.linspace(0.0, x.max()*1.05, 200)
yf = np.polynomial.polynomial.polyval(xf, c)
c1 = np.polynomial.polynomial.polyfit(x, y, 1, w=L)
assert abs(c1[0] - e1) < 1e-12
yf1 = np.polynomial.polynomial.polyval(xf, c1)
print("n=2: %.6f  n=1: %.6f  exact: %.6f" % (e2, e1, e_exact))


def draw(theme):
    col = style.colors(theme)
    fig, ax = style.figure(1)
    ax.plot(xf, yf, **style.reference(theme), label="quadratic fit in $1/L$")
    ax.plot(x, y, "o", color=col["series"][0], label="DMRG, open chains")
    ax.plot([0], [e2], "s", color=col["series"][1], ms=6,
            label=r"extrapolated, $%.5f$" % e2)
    ax.plot([0], [e_exact], "x", color=col["ink"], ms=9, mew=1.6,
            label=r"Bethe ansatz $1/4-\ln 2=%.5f$" % e_exact)
    for Li, xi, yi in zip(L, x, y):
        if Li in (16, 32, 64, 128):
            ax.annotate("L=%d" % Li, (xi, yi), textcoords="offset points",
                        xytext=(4, -12), fontsize=8, color=col["ink2"])
    ax.plot(xf, yf1, color=col["ink2"], ls=":", lw=style.LW_REF,
            label=r"linear fit, $%.5f$" % e1)
    ax.set_xlim(-0.003, x.max()*1.08)
    # zoom on the extrapolation
    ins = ax.inset_axes([0.6, 0.08, 0.37, 0.37])
    ins.plot(xf, yf, **style.reference(theme))
    ins.plot(xf, yf1, color=col["ink2"], ls=":", lw=style.LW_REF)
    ins.plot(x, y, "o", color=col["series"][0])
    ins.plot([0], [e2], "s", color=col["series"][1], ms=6)
    ins.plot([0], [e_exact], "x", color=col["ink"], ms=9, mew=1.6)
    ins.set_xlim(-0.0002, 0.0012)
    ins.set_ylim(e_exact - 5e-5, e_exact + 2.5e-4)
    ins.set_xticks([0, 0.0005, 0.001])
    ins.set_xticklabels(["0", "0.0005", "0.001"])
    ins.set_yticks([e_exact, e_exact + 1e-4, e_exact + 2e-4])
    ins.set_yticklabels([r"$e_\infty$", r"$+10^{-4}$", r"$+2\cdot10^{-4}$"])
    ins.tick_params(labelsize=7)
    ins.set_facecolor("none")
    ax.set_xlabel(r"$1/L$")
    ax.set_ylabel(r"$E_0/L$ [$J$]")
    ax.legend(loc="upper left")
    return fig


style.main(draw, "post-processing_extrapolation")
