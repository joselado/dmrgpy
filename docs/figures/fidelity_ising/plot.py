"""Plot figure 11, fidelity_ising, from data.npz.

Writes docs/guide/img/fidelity_ising.png and -dark.png.
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
lam, ns, chi = d["lam"], d["ns"], d["chi"]


def bulk(l):
    return np.where(l < 1, 1/(16*(1-l**2)), 1/(16*l**2*(l**2-1)))


def draw(theme):
    c = style.colors(theme)
    fig, (a, b) = style.figure(2)
    for k, n in enumerate(ns):
        a.plot(lam, chi[k]/n, "-", color=c["series"][k], label=r"$n=%d$" % n)
    a.plot(d["lam_ed"], d["chi_ed"]/ns[0], **style.reference(theme),
           label=r"ED, $n=%d$" % ns[0])
    a.axvline(1.0, color=c["ink2"], lw=0.8, ls=":")
    a.set_xlabel(r"$\lambda$")
    a.set_ylabel(r"fidelity susceptibility $\chi/n$")
    a.set_xlim(lam[0], lam[-1])
    a.set_ylim(0, 1.1*np.max(chi[-1]/ns[-1]))
    a.legend(loc="upper right")

    # the edge contribution cancels in the difference of the two longest
    # chains, which leaves the bulk value per site
    n1, n2 = ns[-2], ns[-1]
    dens = (chi[-1] - chi[-2])/(n2 - n1)
    b.plot(lam, dens, "o", ms=3.5, color=c["series"][len(ns)-1],
           label=r"$(\chi_{%d}-\chi_{%d})/%d$" % (n2, n1, n2 - n1))
    for k, side in enumerate((lam < 0.97, lam > 1.03)):
        b.plot(lam[side], bulk(lam[side]), **style.reference(theme),
               label="infinite chain, exact" if k == 0 else None)
    b.axvline(1.0, color=c["ink2"], lw=0.8, ls=":")
    b.set_xlabel(r"$\lambda$")
    b.set_ylabel(r"bulk $\chi/n$")
    b.set_xlim(lam[0], lam[-1])
    b.set_ylim(0, 1.2*np.max(dens))
    b.legend(loc="upper right")
    return fig


style.main(draw, "fidelity_ising")
