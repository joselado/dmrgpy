import numpy as np
from matplotlib.colors import LogNorm
import style

d = np.load("data.npz")
ks_a, es_a, S_a = d["ks_a"], d["es_a"], d["S_a"]
ks_b, es_b, S_b = d["ks_b"], d["es_b"], d["S_b"]
band_b, gap = d["e_b"][:, 0], float(d["gap_ref"])
PI_NAMES = {0: "0", 1: r"$\pi/2$", -1: r"$-\pi/2$", 2: r"$\pi$", -2: r"$-\pi$",
            3: r"$3\pi/2$", 4: r"$2\pi$"}


def pi_ticks(ax, k):
    q = np.pi/2
    t = np.arange(np.ceil(k[0]/q - 1e-9), np.floor(k[-1]/q + 1e-9) + 1)
    ax.set_xticks(t*q)
    ax.set_xticklabels([PI_NAMES[int(v)] for v in t])


def draw(theme):
    fig, (ax1, ax2) = style.figure(2)
    im = ax1.pcolormesh(ks_a, es_a, S_a.T, shading="auto", cmap=style.CMAP,
                        rasterized=True)
    ax1.plot(ks_a, d["band_a"], **style.heatmap_line(), label="VUMPS band")
    ax1.plot(ks_a, d["exact_a"], color="white", ls=":", lw=2.0,
             label="exact band")
    ax1.legend(loc="lower center", labelcolor="white")
    style.colorbar(fig, im, ax1, r"$S^{xx}(k,\omega)$", theme)

    vmax = S_b.max()
    im = ax2.pcolormesh(ks_b, es_b, np.maximum(S_b.T, 1e-12), shading="auto",
                        cmap=style.CMAP, rasterized=True,
                        norm=LogNorm(vmin=1e-4*vmax, vmax=vmax))
    ax2.plot(ks_b, band_b, **style.heatmap_line(), label="magnon branch")
    ax2.axhline(gap, color="white", ls=":", lw=1.2,
                label=r"Haldane gap $\Delta=%.4f$" % gap)
    ax2.legend(loc="upper center", labelcolor="white")
    style.colorbar(fig, im, ax2, r"$S^{zz}(k,\omega)$", theme)
    for ax, k in ((ax1, ks_a), (ax2, ks_b)):
        pi_ticks(ax, k)
        ax.set_xlabel(r"momentum $k$")
        ax.set_ylabel(r"$\omega$ [$J$]")
    fig.tight_layout()
    return fig


style.main(draw, "infinite-chains_excitations")
