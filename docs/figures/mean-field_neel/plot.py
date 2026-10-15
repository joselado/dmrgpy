import numpy as np
import style

d = np.load("data.npz")
n, p, mz = int(d["n"]), d["p"], d["mz"]
e_neel, e_exact = -(n-1)/(4*n), float(d["e_exact"])/n
# staggered magnetization of the 12 central sites
centre = np.arange((n - 12)//2, (n + 12)//2)
ms = np.mean((-1.0)**centre*mz[:, centre], axis=1)
show = [0.0, 0.6, 0.8, 0.9, 1.0]


def draw(theme):
    c = style.colors(theme)
    fig, (a, b) = style.figure(2)
    shown = [j for j, pj in enumerate(p)
             if np.min(np.abs(pj - np.array(show))) < 1e-9]
    for k, j in enumerate(shown):
        a.plot(np.arange(n), mz[j], "o-", ms=3.5, color=c["series"][k],
               label=r"$p=%g$" % p[j])
    a.set_xlabel(r"site $i$")
    a.set_ylabel(r"$\langle S^z_i\rangle$")
    a.set_ylim(-0.62, 0.75)
    a.legend(loc="upper center", ncol=5, fontsize=8, handlelength=1.2,
             columnspacing=0.8)

    b.plot(p, ms, "o-", color=c["series"][0],
           label=r"staggered $m_s$ (left axis)")
    b.set_xlabel(r"fraction of exchange kept exactly, $p$")
    b.set_ylabel(r"staggered magnetization $m_s$", color=c["series"][0])
    b.set_ylim(-0.02, 0.55)
    b2 = b.twinx()
    b2.spines["right"].set_visible(True)
    b2.spines["top"].set_visible(False)
    b2.plot(p, d["e_var"]/n, "s", color=c["series"][1],
            label=r"energy $\langle H\rangle/n$ (right axis)")
    for e, name in ((e_neel, "Neel"), (e_exact, "exact chain")):
        b2.axhline(e, **style.reference(theme))
        b2.text(0.02, e, name, va="bottom", fontsize=8, color=c["ink"])
    b2.set_ylabel(r"energy per site [$J$]", color=c["series"][1])
    h1, l1 = b.get_legend_handles_labels()
    h2, l2 = b2.get_legend_handles_labels()
    b.legend(h1 + h2, l1 + l2, loc="lower left", bbox_to_anchor=(0.0, 0.08))
    return fig


style.main(draw, "mean-field_neel")
