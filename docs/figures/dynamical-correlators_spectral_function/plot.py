import numpy as np
import style

d = np.load("data.npz")
es, A, Us = d["es"], d["A"], d["Us_spec"]
off = 0.8*np.max(A)  # vertical offset between the spectra


def draw(theme):
    c = style.colors(theme)
    fig, (a1, a2) = style.figure(2)
    for k, u in enumerate(Us):
        a1.plot(es, A[k] + k*off, color=c["series"][k], lw=style.LW,
                label=r"$U=%g\,t$" % u)
        if k == 0:
            a1.plot(es, d["A_free"], **style.reference(theme),
                    label="exact, $U=0$")
    a1.axvline(0.0, color=c["ink2"], lw=0.8, ls=":")
    a1.set_xlim(es[0], es[-1])
    a1.set_yticks([])
    a1.set_xlabel(r"$\omega-\mu$ [$t$]")
    a1.set_ylabel(r"$A_{ii}(\omega)$ [$1/t$], offset")
    a1.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=3, fontsize=8)

    a2.plot(d["Ul"], d["gap_lw"], **style.reference(theme),
            label="Lieb-Wu gap, $n=\\infty$")
    a2.plot(d["Ul"], d["Ul"]/2, color=c["ink2"], ls=":", lw=style.LW_REF,
            label=r"$U/2$")
    for k, (L, m) in enumerate(zip(d["ns_gap"], "os")):
        a2.plot(d["Us_gap"], d["gap"][k], m, color=c["series"][0],
                mfc="none" if k else c["series"][0],
                label="charge gap, $n=%d$" % L)
    a2.plot(d["Us_gap"], d["mu"][0], "D", color=c["series"][1], ms=4,
            label=r"$\mu$, $n=%d$" % d["ns_gap"][0])
    a2.set_xlabel(r"$U$ [$t$]")
    a2.set_ylabel(r"energy [$t$]")
    a2.set_xlim(0, None)
    a2.set_ylim(0, None)
    a2.legend(loc="upper left")
    fig.tight_layout()
    return fig


style.main(draw, "dynamical-correlators_spectral_function")
