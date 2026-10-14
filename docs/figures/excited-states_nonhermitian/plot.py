"""Plot figure 5, excited-states_nonhermitian, from data.npz.

Left: the right-vector density <n_i>_R of the n-site chain for a few g,
NH-DMRG against the exact density (Hermitian DMRG plus the similarity
transform). Right: the charge moved into the left half,
Q_L = sum_{i<n/2} <n_i>_R - n/4, against g, at n sites (NH-DMRG and the
exact transform) and at n_ed sites (NH-DMRG and ED).

Writes docs/guide/img/excited-states_nonhermitian.png and -dark.png.
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
gam = d["gamma"]
conv = d["converged"].astype(bool)   # a failed NH-DMRG solve is left out
SHOW = [0.0, 0.2, 0.4]               # the g of the left panel


def left_charge(dens):
    n = dens.shape[-1]
    return dens[..., :n//2].sum(axis=-1) - n/4


def draw(theme):
    c = style.colors(theme)
    fig, (a1, a2) = style.figure(2)
    shown = [k for k, g in enumerate(gam) if conv[k] and
             np.min(np.abs(g - np.array(SHOW))) < 1e-9]
    for j, k in enumerate(shown):
        i = np.arange(params["n"])
        a1.plot(i, d["dens"][k], "o-", color=c["series"][j], ms=4,
                label=r"$\gamma=%.1f$" % gam[k])
        a1.plot(i, d["dens_ref"][k], **style.reference(theme),
                label="exact" if j == len(shown)-1 else None)
    a1.set_xlabel(r"site $i$")
    a1.set_ylabel(r"$\langle n_i\rangle_R$")
    a1.legend(loc="lower left", ncol=2)
    a1.text(0.98, 0.98, r"$n=%d$" % params["n"], transform=a1.transAxes,
            ha="right", va="top")

    a2.plot(gam[conv], left_charge(d["dens"])[conv], "o-",
            color=c["series"][0], label=r"NH-DMRG, $n=%d$" % params["n"])
    a2.plot(gam, left_charge(d["dens_ref"]), **style.reference(theme),
            label="exact")
    a2.plot(gam, left_charge(d["dens_small"]), "s-", color=c["series"][1],
            label=r"NH-DMRG, $n=%d$" % params["n_ed"])
    a2.plot(gam, left_charge(d["dens_ed"]), color=c["ink"], ls=":",
            lw=style.LW_REF, label=r"ED, $n=%d$" % params["n_ed"])
    a2.set_xlabel(r"$\gamma$ [$t$]")
    a2.set_ylabel(r"charge in the left half, $Q_L-n/4$")
    a2.legend(loc="upper left")
    return fig


style.main(draw, "excited-states_nonhermitian")
