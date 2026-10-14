"""Plot figure 19, finite-temperature_specific_heat, from data.npz.

Writes docs/guide/img/finite-temperature_specific_heat.png and -dark.png.
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
T, Te, Tm = d["T"], d["T_exact"], d["T_metts"]


def exact_at(arr, t):
    return np.interp(np.log(t), np.log(Te), arr)


print("max |E_pur-E_exact|/n = %.2e (exact interpolated on T_exact; compute.log has pointwise values)" % np.max(np.abs(d["E_pur"] - exact_at(d["E_exact"], T))/n))
print("max |C_pur-C_exact|/n = %.2e (exact interpolated)" % np.max(np.abs(d["C_pur"] - exact_at(d["C_exact"], T))/n))
print("METTS deviations in error bars:",
      np.round(np.abs(d["E_metts"] - exact_at(d["E_exact"], Tm))/d["err_metts"], 2))
k = np.argmax(d["C_exact"])
print("exact C peak: C/n = %.4f at T = %.3f" % (d["C_exact"][k]/n, Te[k]))


def draw(theme):
    c = style.colors(theme)
    fig, (ax1, ax2) = style.figure(2)
    ax1.plot(Te, d["E_exact"]/n, **style.reference(theme), label="exact")
    ax1.plot(T, d["E_pur"]/n, "o", color=c["series"][0], ms=style.MS, label="purification")
    ax1.errorbar(Tm, d["E_metts"]/n, yerr=d["err_metts"]/n, fmt="s",
                 color=c["series"][1], ms=style.MS, capsize=3, label="METTS")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$T$ [$J$]")
    ax1.set_ylabel(r"$\langle H\rangle_T/n$ [$J$]")
    ax1.legend()
    ax2.plot(Te, d["C_exact"]/n, **style.reference(theme), label="exact")
    ax2.plot(T, d["C_pur"]/n, "o", color=c["series"][0], ms=style.MS, label="purification")
    ax2.set_xscale("log")
    ax2.set_xlabel(r"$T$ [$J$]")
    ax2.set_ylabel(r"$C(T)/n$")
    ax2.legend()
    for ax in (ax1, ax2):
        ax.set_xlim(Te[0], Te[-1])
    return fig


style.main(draw, "finite-temperature_specific_heat")
