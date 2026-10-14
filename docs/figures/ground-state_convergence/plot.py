"""Plot figure 1, ground-state_convergence, from data.npz.

Writes docs/guide/img/ground-state_convergence.png and -dark.png.
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
maxm = d["maxm"]
# per site; abs() because at the ED/DMRG floor the error can come out
# marginally negative
err = np.abs(d["energy"] - d["e_ed"]) / n
de = d["de"] / n

# dE^2 times a constant, fitted in log space on the points above the
# numerical floor of the energy, about 1e-13 in total (energy error per
# site > FLOOR); at maxm=64 and 96 both curves sit on their floors
FLOOR = 5e-14
fit = err > FLOOR
logc = np.mean(np.log(err[fit]) - np.log(de[fit]**2))
print("fit points maxm =", maxm[fit], " constant = %.3g" % np.exp(logc))


def draw(theme):
    c = style.colors(theme)
    fig, ax = style.figure(1)
    ax.plot(maxm, err, "o-", color=c["series"][0],
            label=r"energy error $(E_0-E_0^{\rm ED})/n$")
    ax.plot(maxm, de, "s-", color=c["series"][1],
            label=r"energy fluctuation $\delta E/n$")
    ax.plot(maxm[fit], np.exp(logc) * de[fit]**2, **style.reference(theme),
            label=r"$\propto(\delta E/n)^2$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(maxm)
    ax.set_xticklabels([str(m) for m in maxm])
    ax.minorticks_off()
    ax.set_xlabel(r"bond dimension $\chi$ (maxm)")
    ax.set_ylabel(r"per site [$J$]")
    ax.legend()
    return fig


style.main(draw, "ground-state_convergence")
