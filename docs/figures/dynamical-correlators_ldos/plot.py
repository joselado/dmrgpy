import numpy as np
import style

d = np.load("data.npz")
es, A, A_free, Vs, sites = d["es"], d["A"], d["A_free"], d["Vs"], d["sites"]


def draw(theme):
    c = style.colors(theme)
    fig, axs = style.figure(2, sharey=True)
    for j, (ax, i) in enumerate(zip(axs, sites)):
        for k, V in enumerate(Vs):
            ax.plot(es, A[k, j], color=c["series"][k], lw=style.LW,
                    label=r"$V=%g\,t$" % V)
        ax.plot(es, A_free[j], **style.reference(theme), label="exact, $V=0$")
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


style.main(draw, "dynamical-correlators_ldos")
