"""Plot figure 4, excited-states_haldane, from data.npz.

Writes docs/guide/img/excited-states_haldane.png and -dark.png.
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
ns = d["ns"]
x = 1.0/ns**2
gap_inf, slope = float(d["gap_inf"]), float(d["gap_slope"])
xi = float(d["xi_ref"])
split = np.abs(d["split"])
# the prefactor of exp(-n/xi) is fitted in log space on the sizes whose
# splitting is above 1e-9; the energies of one sector moved by 5e-6 at
# n=48 between maxm=100 and 150, so the 2e-10 of n=128 is left out
ok = split > 1e-9
pref = np.exp(np.mean(np.log(split[ok]) + ns[ok]/xi))
prof = d["profile"]


def draw(theme):
    c = style.colors(theme)
    fig, (a1, a2) = style.figure(2)
    # the fit uses n >= fit_min_n; drawn a bit beyond them
    xx = np.linspace(0, 2.0/params["fit_min_n"]**2, 50)
    a1.plot(x, d["gap"], "o", color=c["series"][0], label="DMRG, sectors")
    a1.plot(xx, gap_inf + slope*xx, **style.reference(theme), label=r"fit $a+b/n^2$ ($n\geq%d$), $a=%.4f$"
            % (params["fit_min_n"], gap_inf))
    a1.axhline(float(d["gap_ref"]), color=c["ink"], ls=":", lw=style.LW_REF,
               label=r"$\Delta=0.4105\,J$")
    a1.set_xlabel(r"$1/n^2$")
    a1.set_ylabel(r"gap $E_0(S^z{=}2)-E_0(S^z{=}1)$ [$J$]")
    a1.set_xlim(left=0)
    a1.legend(loc="upper left")
    ins = a1.inset_axes([0.62, 0.22, 0.35, 0.35])
    ins.plot(ns[ok], split[ok], "s", color=c["series"][1], ms=4)
    nn = np.linspace(ns.min(), ns[ok].max(), 50)
    ins.plot(nn, pref*np.exp(-nn/xi), **style.reference(theme))
    ins.text(0.95, 0.92, r"$e^{-n/6.03}$", transform=ins.transAxes,
             ha="right", va="top", fontsize=9)
    ins.set_yscale("log")
    ins.set_xlabel(r"$n$", fontsize=9)
    ins.set_ylabel(r"edge splitting [$J$]", fontsize=9)
    ins.tick_params(labelsize=8)
    i = np.arange(len(prof))
    a2.plot(i, prof, "o-", color=c["series"][1], ms=3.5)
    a2.axhline(0, color=c["ink2"], lw=0.6)
    a2.set_xlabel(r"site $i$")
    a2.set_ylabel(r"$\langle S^z_i\rangle$ in $S^z_{\rm tot}=1$")
    a2.text(0.5, 0.95, r"$n=%d$" % len(prof), transform=a2.transAxes,
            ha="center", va="top")
    return fig


style.main(draw, "excited-states_haldane")
