import numpy as np
import matplotlib.ticker as ticker
import style

d = np.load("data.npz")
U, t = d["U"], d["t"]
j_exact = 0.5*(np.sqrt(U**2 + 16*t**2) - U)  # triplet-singlet splitting


def draw(theme):
    c = style.colors(theme)
    fig, ax = style.figure(1)
    ax.plot(U, d["jfit"][:, 2], "o", color=c["series"][0],
            label=r"get_heff, singly occupied projector")
    ax.plot(U, d["jfit_noproj"][:, 2], "s", mfc="none",
            color=c["series"][1], label=r"get_heff, no projector")
    ax.plot(U, j_exact, **style.reference(theme), label=r"exact, $E_T-E_S$")
    ax.plot(U, 4*t**2/U, color=c["ink2"], ls="-.", lw=style.LW_REF,
            label=r"$4t^2/U$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([1, 2, 5, 10, 20])
    ax.set_yticks([0.2, 0.5, 1, 2, 4])
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        axis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel(r"$U$ [$t$]")
    ax.set_ylabel(r"exchange $J$ [$t$]")
    ax.legend()
    return fig


style.main(draw, "excited-states_heff")
