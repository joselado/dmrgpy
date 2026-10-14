"""Plot figure 16, real-time_light_cone, from data.npz.

Writes docs/guide/img/real-time_light_cone.png and -dark.png.
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
ts, sz, sz_ex, x0 = d["ts"], d["sz"], d["sz_exact"], float(d["x0"])
x = np.arange(n)
dt_row = ts[1] - ts[0]
for t in (4.0, 8.0, 12.0):
    k = int(np.argmin(np.abs(ts - t)))
    print("t=%4.1f  max |TDVP-exact| = %.2e" % (ts[k], np.abs(sz[k] - sz_ex[k]).max()))
print("all times: %.2e" % np.abs(sz - sz_ex).max())


def draw(theme):
    c = style.colors(theme)
    fig, (ax1, ax2) = style.figure(2)
    # signed quantity: diverging map, symmetric about zero; the cone lines
    # are dark, since the middle of RdBu is white
    im = ax1.imshow(sz, origin="lower", aspect="auto", cmap="RdBu_r",
                    vmin=-0.5, vmax=0.5,
                    extent=[-0.5, n - 0.5, ts[0] - dt_row/2, ts[-1] + dt_row/2])
    tt = np.linspace(0, ts[-1], 50)
    for s in (+1, -1):
        ax1.plot(x0 + s*tt, tt, color="#0b0b0b", ls="--", lw=style.LW_REF)
    ax1.set_xlim(-0.5, n - 0.5)
    ax1.set_ylim(0, ts[-1])
    ax1.set_xlabel(r"site $x$")
    ax1.set_ylabel(r"time $t$ [$1/J$]")
    style.colorbar(fig, im, ax1, r"$\langle S^z_x(t)\rangle$", theme)

    for j, t in enumerate((4.0, 8.0, 12.0)):
        k = int(np.argmin(np.abs(ts - t)))
        ax2.plot(x - x0, sz[k], "o", color=c["series"][j], ms=style.MS,
                 label=r"TDVP, $t=%g$" % t)
        ax2.plot(x - x0, sz_ex[k], **style.reference(theme),
                 label="exact" if j == 0 else None)
    t = 12.0
    u = np.linspace(-t, t, 200)
    ax2.plot(u, -np.arcsin(u/t)/np.pi, color=c["ink2"], ls=":", lw=2.4, zorder=0,
             label=r"$-\frac{1}{\pi}\arcsin(x/t)$, $t=12$")
    ax2.set_xlabel(r"$x-x_0$")
    ax2.set_ylabel(r"$\langle S^z_x\rangle$")
    ax2.legend(loc="lower left")
    fig.subplots_adjust(wspace=0.4)
    return fig


style.main(draw, "real-time_light_cone")
