"""Plot figure 15, correlator-methods_comparison, from data.npz.

Writes docs/guide/img/correlator-methods_comparison.png and -dark.png.
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
methods = [str(m) for m in d["methods"]]
es, es_cvm = d["es"], d["es_cvm"]
offset = 1.05 * max(np.max(d["y_" + m]) for m in methods)


def reference(m):
    """ED reference of a method: KPM against ED with the same kernel."""
    if m == "CVM":
        return d["ed_cvm"]
    return d["ed_kpm"] if m == "KPM" else d["ed"]


def residual(m):
    ref = reference(m)
    return np.max(np.abs(d["y_" + m] - ref)) / np.max(np.abs(ref))


for m in methods:
    print("%-7s residual %.3f  time %.1f s" % (m, residual(m), d["time_" + m]))


def draw(theme):
    c = style.colors(theme)
    fig, ax = style.figure(1)
    fig.set_size_inches(6.4, 6.4)
    for k, m in enumerate(methods):
        sh = (len(methods) - 1 - k) * offset
        x = es_cvm if m == "CVM" else es
        if m == "KPM":
            ax.plot(es, d["ed"] + sh, color=c["ink2"], ls=":",
                    lw=style.LW_REF, label="ED, Lorentzian lines (KPM row)")
            ax.plot(es, d["ed_kpm"] + sh, **style.reference(theme), zorder=3)
        else:
            ax.plot(es, d["ed"] + sh, **style.reference(theme), zorder=3,
                    label="ED (for KPM: the same KPM kernel)" if k == 1 else None)
        ax.plot(x, d["y_" + m] + sh, color=c["series"][k], lw=style.LW,
                marker="o" if m == "CVM" else None, ms=2.5)
        ax.text(es[-1], sh + 0.12*offset,
                "%s\nmax. deviation %.1f%%" % (m, 100*residual(m)),
                ha="right", va="bottom", color=c["series"][k])
    ax.set_xlim(es[0], es[-1])
    ax.set_xlabel(r"$\omega$ [$J$]")
    ax.set_ylabel(r"$S^{zz}_{00}(\omega)$ [$1/J$], offset")
    ax.set_yticks([])
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2)
    return fig


style.main(draw, "correlator-methods_comparison")
