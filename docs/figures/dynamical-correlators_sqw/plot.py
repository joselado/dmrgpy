"""Plot figure 13, dynamical-correlators_sqw, from data.npz.

Writes docs/guide/img/dynamical-correlators_sqw.png and -dark.png.
    python3 plot.py [--theme light|dark]
"""
import os
import sys

import numpy as np
from matplotlib.colors import PowerNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import common
import style

params, d = common.load(os.path.join(HERE, "data.npz"))
qs, es, sqw = d["qs"], d["es"], d["sqw"]
sqw = np.clip(sqw, 0.0, None)          # Jackson ringing below zero, tiny
rel = np.abs(d["sumrule"] - d["sq"]) / d["sq"]
print("sum rule: max relative deviation %.2e" % rel.max())

qq = np.linspace(0, np.pi, 300)
lower = np.pi/2*np.abs(np.sin(qq))
upper = np.pi*np.abs(np.sin(qq/2))

cuts = [int(np.argmin(np.abs(qs - np.pi/2))), len(qs) - 1]


def draw(theme):
    c = style.colors(theme)
    fig, (a1, a2) = style.figure(2, gridspec_kw={"width_ratios": [1.15, 1]})
    vmax = np.percentile(sqw, 99.5)
    im = a1.pcolormesh(qs/np.pi, es, sqw.T, cmap=style.CMAP, shading="nearest",
                       norm=PowerNorm(0.5, vmin=0, vmax=vmax), rasterized=True)
    a1.plot(qq/np.pi, lower, **style.heatmap_line())
    a1.plot(qq/np.pi, upper, **style.heatmap_line())
    a1.set_xlim(0, 1)
    a1.set_ylim(es[0], es[-1])
    a1.set_xlabel(r"$q$ [$\pi$]")
    a1.set_ylabel(r"$\omega$ [$J$]")
    style.colorbar(fig, im, a1, r"$S^{zz}(q,\omega)$ [$1/J$]", theme)
    big = np.max(sqw[cuts[-1]])
    for k, idx in enumerate(cuts):
        q = qs[idx]
        col = c["series"][k]
        f = max(1, int(round(0.5*big/np.max(sqw[idx])))) if k == 0 else 1
        a2.plot(es, f*sqw[idx], color=col, lw=style.LW,
                label=r"$q=%.2f\pi$%s, $\int S\,d\omega$=%.4f (exact %.4f)"
                % (q/np.pi, r" ($\times%d$)" % f if f > 1 else "",
                   d["sumrule"][idx], d["sq"][idx]))
        a2.axvline(np.pi/2*np.sin(q), color=col, ls=":", lw=style.LW_REF)
    a2.set_xlim(es[0], es[-1])
    a2.set_ylim(0, None)
    a2.set_xlabel(r"$\omega$ [$J$]")
    a2.set_ylabel(r"$S^{zz}(q,\omega)$ [$1/J$]")
    a2.legend(loc="upper right")
    fig.tight_layout()
    return fig


style.main(draw, "dynamical-correlators_sqw")
