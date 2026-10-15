import numpy as np
import style

d = np.load("data.npz")
n = d["n"]
T, Te, Tm = d["T"], d["T_exact"], d["T_metts"]


def draw(theme):
    c = style.colors(theme)
    fig, (ax1, ax2) = style.figure(2)
    ax1.plot(Te, d["E_exact"]/n, **style.reference(theme), label="exact")
    ax1.plot(T, d["E_pur"]/n, "o", color=c["series"][0], ms=style.MS,
             label="purification")
    ax1.errorbar(Tm, d["E_metts"]/n, yerr=d["err_metts"]/n, fmt="s",
                 color=c["series"][1], ms=style.MS, capsize=3, label="METTS")
    ax1.set_ylabel(r"$\langle H\rangle_T/n$ [$J$]")
    ax2.plot(Te, d["C_exact"]/n, **style.reference(theme), label="exact")
    ax2.plot(T, d["C_pur"]/n, "o", color=c["series"][0], ms=style.MS,
             label="purification")
    ax2.set_ylabel(r"$C(T)/n$")
    for ax in (ax1, ax2):
        ax.set_xscale("log")
        ax.set_xlabel(r"$T$ [$J$]")
        ax.set_xlim(Te[0], Te[-1])
        ax.legend()
    return fig


style.main(draw, "finite-temperature_specific_heat")
