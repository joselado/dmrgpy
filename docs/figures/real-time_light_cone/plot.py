import numpy as np
import style

d = np.load("data.npz")
n, ts, sz, sz_ex = int(d["n"]), d["ts"], d["sz"], d["sz_exact"]
x0 = n/2 - 0.5  # the wall, between sites n/2-1 and n/2
x = np.arange(n)
dt_row = ts[1] - ts[0]


def draw(theme):
    c = style.colors(theme)
    fig, (ax1, ax2) = style.figure(2)
    # signed quantity: diverging map, dark cone lines on its white middle
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
    u = np.linspace(-12.0, 12.0, 200)
    ax2.plot(u, -np.arcsin(u/12.0)/np.pi, color=c["ink2"], ls=":", lw=2.4,
             zorder=0, label=r"$-\frac{1}{\pi}\arcsin(x/t)$, $t=12$")
    ax2.set_xlabel(r"$x-x_0$")
    ax2.set_ylabel(r"$\langle S^z_x\rangle$")
    ax2.legend(loc="lower left")
    fig.subplots_adjust(wspace=0.4)
    return fig


style.main(draw, "real-time_light_cone")
