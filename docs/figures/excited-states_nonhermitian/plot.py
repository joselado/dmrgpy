import numpy as np
import style

d = np.load("data.npz")
gam, n, n_ed = d["gamma"], int(d["n"]), int(d["n_ed"])
conv = d["converged"].astype(bool)  # a failed NH-DMRG solve is left out
SHOW = [0.0, 0.2, 0.4]              # the g of the left panel


def left_charge(dens):  # charge moved into the left half
    m = dens.shape[-1]
    return dens[..., :m//2].sum(axis=-1) - m/4


def draw(theme):
    c = style.colors(theme)
    fig, (a1, a2) = style.figure(2)
    shown = [k for k, g in enumerate(gam) if conv[k] and
             np.min(np.abs(g - np.array(SHOW))) < 1e-9]
    i = np.arange(n)
    for j, k in enumerate(shown):
        a1.plot(i, d["dens"][k], "o-", color=c["series"][j], ms=4,
                label=r"$\gamma=%.1f$" % gam[k])
        a1.plot(i, d["dens_ref"][k], **style.reference(theme),
                label="exact" if j == len(shown)-1 else None)
    a1.set_xlabel(r"site $i$")
    a1.set_ylabel(r"$\langle n_i\rangle_R$")
    a1.legend(loc="lower left", ncol=2)
    a1.text(0.98, 0.98, r"$n=%d$" % n, transform=a1.transAxes,
            ha="right", va="top")

    a2.plot(gam[conv], left_charge(d["dens"])[conv], "o-",
            color=c["series"][0], label=r"NH-DMRG, $n=%d$" % n)
    a2.plot(gam, left_charge(d["dens_ref"]), **style.reference(theme),
            label="exact")
    a2.plot(gam, left_charge(d["dens_small"]), "s-", color=c["series"][1],
            label=r"NH-DMRG, $n=%d$" % n_ed)
    a2.plot(gam, left_charge(d["dens_ed"]), color=c["ink"], ls=":",
            lw=style.LW_REF, label=r"ED, $n=%d$" % n_ed)
    a2.set_xlabel(r"$\gamma$ [$t$]")
    a2.set_ylabel(r"charge in the left half, $Q_L-n/4$")
    a2.legend(loc="upper left")
    return fig


style.main(draw, "excited-states_nonhermitian")
