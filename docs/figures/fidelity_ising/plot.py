import numpy as np
import style

d = np.load("data.npz")
lam, ns, chi = d["lam"], d["ns"], d["chi"]


def bulk(l):  # infinite chain, per site
    return np.where(l < 1, 1/(16*(1-l**2)), 1/(16*l**2*(l**2-1)))


# the edge contribution cancels in the difference of the two longest
# chains, which leaves the bulk value per site
n1, n2 = ns[-2], ns[-1]
dens = (chi[-1] - chi[-2])/(n2 - n1)


def draw(theme):
    c = style.colors(theme)
    fig, (a, b) = style.figure(2)
    for k, n in enumerate(ns):
        a.plot(lam, chi[k]/n, "-", color=c["series"][k], label=r"$n=%d$" % n)
    a.plot(lam, d["chi_ed"]/ns[0], **style.reference(theme),
           label=r"ED, $n=%d$" % ns[0])
    a.set_ylabel(r"fidelity susceptibility $\chi/n$")
    a.set_ylim(0, 1.1*np.max(chi[-1]/ns[-1]))
    b.plot(lam, dens, "o", ms=3.5, color=c["series"][len(ns)-1],
           label=r"$(\chi_{%d}-\chi_{%d})/%d$" % (n2, n1, n2 - n1))
    for k, side in enumerate((lam < 0.97, lam > 1.03)):
        b.plot(lam[side], bulk(lam[side]), **style.reference(theme),
               label="infinite chain, exact" if k == 0 else None)
    b.set_ylabel(r"bulk $\chi/n$")
    b.set_ylim(0, 1.2*np.max(dens))
    for ax in (a, b):
        ax.axvline(1.0, color=c["ink2"], lw=0.8, ls=":")
        ax.set_xlabel(r"$\lambda$")
        ax.set_xlim(lam[0], lam[-1])
        ax.legend(loc="upper right")
    return fig


style.main(draw, "fidelity_ising")
