"""Plot figure 20, infinite-chains_energy, from data.npz.

Writes docs/guide/img/infinite-chains_energy.png and -dark.png.
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
maxm = d["maxm"]
err = np.abs(d["energy"] - d["e_exact"])
r = d["r"]
c = np.abs(d["czz"])
# Luttinger form a sqrt(ln r)/r, prefactor fitted in log space on r = 8..20
FIT = (r >= 8) & (r <= 20)
a = np.exp(np.mean(np.log(c[FIT]) - np.log(np.sqrt(np.log(r[FIT]))/r[FIT])))
rr = np.linspace(2, r.max(), 200)
print("fitted a = %.4f, field theory (2 pi)^(-3/2) = %.4f" % (a, (2*np.pi)**-1.5))
e_uc2 = abs(float(d["e_uc2"]) - float(d["e_exact"]))
print("energy error n_uc=1: %s" % err)
print("energy error n_uc=2 at maxm=%d: %.2e" % (maxm[-1], e_uc2))
print("3<Sz0 Sz1> averaged: %.9f, exact e/3*3: %.9f"
      % (3*d["czz"][0], float(d["e_exact"])))


def draw(theme):
    col = style.colors(theme)
    fig, (ax1, ax2) = style.figure(2)
    ax1.plot(maxm, err, "o-", color=col["series"][0],
             label="iDMRG, one-site cell")
    ax1.plot([maxm[-1]], [e_uc2], "s", color=col["series"][1],
             label="iDMRG, two-site cell")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xticks(maxm)
    ax1.set_xticklabels([str(m) for m in maxm])
    ax1.minorticks_off()
    ax1.set_xlabel(r"bond dimension $\chi$ (maxm)")
    ax1.set_ylabel(r"$|e(\chi)-e_\infty|$ [$J$]")
    ax1.legend()
    ax1.set_ylim(top=err.max()*3)

    ax2.plot(r, c, "o", color=col["series"][0], ms=4,
             label=r"iDMRG, $\chi=%d$, two-site cell" % maxm[-1])
    ax2.plot(rr, a*np.sqrt(np.log(rr))/rr, **style.reference(theme),
             label=r"$a\sqrt{\ln r}/r$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"distance $r$")
    ax2.set_ylabel(r"$|\langle S^z_0S^z_r\rangle|$")
    ax2.legend()
    fig.tight_layout()
    return fig


style.main(draw, "infinite-chains_energy")
