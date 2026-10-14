"""Plot figure 17, density-of-states_ldos, from data.npz.

Writes docs/guide/img/density-of-states_ldos.png and -dark.png.
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
es, A, Vs, sites = d["es"], d["A"], d["Vs"], d["sites"]
for j, i in enumerate(sites):
    r = np.max(np.abs(A[0, j] - d["A_free"][j])) / np.max(d["A_free"][j])
    print("site %d: V=0 max deviation from free fermions %.2e" % (i, r))
print("gaps:", d["gap"], " free:", d["gap_free"], " captured:", d["captured"])


def draw(theme):
    c = style.colors(theme)
    fig, axs = style.figure(2, sharey=True)
    for j, (ax, i) in enumerate(zip(axs, sites)):
        for k, V in enumerate(Vs):
            ax.plot(es, A[k, j], color=c["series"][k], lw=style.LW,
                    label=r"$V=%g\,t$" % V)
        ax.plot(es, d["A_free"][j], **style.reference(theme),
                label="exact, $V=0$")
        ax.axvline(0.0, color=c["ink2"], lw=0.8, ls=":")
        ax.set_xlim(es[0], es[-1])
        ax.set_ylim(0, None)
        ax.set_xlabel(r"$\omega-\mu$ [$t$]")
        where = "bulk" if i > 0 else "edge"
        ax.text(0.03, 0.95, "%s site, $i=%d$" % (where, i),
                transform=ax.transAxes, va="top", color=c["ink"])
    axs[0].set_ylabel(r"$A_i(\omega)$ [$1/t$]")
    axs[1].legend(loc="upper right")
    fig.tight_layout()
    return fig


style.main(draw, "density-of-states_ldos")
