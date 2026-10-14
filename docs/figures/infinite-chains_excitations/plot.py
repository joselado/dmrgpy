"""Plot figure 21, infinite-chains_excitations, from data.npz.

Writes docs/guide/img/infinite-chains_excitations.png and -dark.png.
    python3 plot.py [--theme light|dark]
"""
import os
import sys

import numpy as np
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import common
import style

params, d = common.load(os.path.join(HERE, "data.npz"))
ks_a, es_a, S_a = d["ks_a"], d["es_a"], d["S_a"]
ks_b, es_b, S_b = d["ks_b"], d["es_b"], d["S_b"]
band_b = d["e_b"][:, 0]
gap = float(d["gap_ref"])
print("(a) max |E-exact| = %.2e" % np.max(np.abs(d["band_a"] - d["exact_a"])))
ipi = int(np.argmin(np.abs(ks_b - np.pi)))
print("(b) E(pi) = %.6f, reference %.6f" % (band_b[ipi], gap))


def pi_ticks(ax, lo, hi):
    t = np.arange(np.ceil(lo/(np.pi/2) - 1e-9), np.floor(hi/(np.pi/2) + 1e-9) + 1)
    names = {0: "0", 1: r"$\pi/2$", -1: r"$-\pi/2$", 2: r"$\pi$",
             -2: r"$-\pi$", 3: r"$3\pi/2$", 4: r"$2\pi$"}
    ax.set_xticks(t*np.pi/2)
    ax.set_xticklabels([names[int(v)] for v in t])


def draw(theme):
    fig, (ax1, ax2) = style.figure(2)
    im = ax1.pcolormesh(ks_a, es_a, S_a.T, shading="auto", cmap=style.CMAP,
                        rasterized=True)
    ax1.plot(ks_a, d["band_a"], **style.heatmap_line(), label="VUMPS band")
    ax1.plot(ks_a, d["exact_a"], color="white", ls=":", lw=2.0,
             label="exact band")
    pi_ticks(ax1, ks_a[0], ks_a[-1])
    ax1.set_xlabel(r"momentum $k$")
    ax1.set_ylabel(r"$\omega$ [$J$]")
    ax1.legend(loc="lower center", labelcolor="white")
    style.colorbar(fig, im, ax1, r"$S^{xx}(k,\omega)$", theme)

    vmax = S_b.max()
    im = ax2.pcolormesh(ks_b, es_b, np.maximum(S_b.T, 1e-12), shading="auto",
                        cmap=style.CMAP, norm=LogNorm(vmin=1e-4*vmax, vmax=vmax),
                        rasterized=True)
    ax2.plot(ks_b, band_b, **style.heatmap_line(), label="magnon branch")
    ax2.axhline(gap, color="white", ls=":", lw=1.2,
                label=r"Haldane gap $\Delta=%.4f$" % gap)
    pi_ticks(ax2, ks_b[0], ks_b[-1])
    ax2.set_xlabel(r"momentum $k$")
    ax2.set_ylabel(r"$\omega$ [$J$]")
    ax2.legend(loc="upper center", labelcolor="white")
    style.colorbar(fig, im, ax2, r"$S^{zz}(k,\omega)$", theme)
    fig.tight_layout()
    return fig


style.main(draw, "infinite-chains_excitations")
