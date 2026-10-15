import numpy as np
import style

d = np.load("data.npz")
ns, prof = d["ns"], d["profile"]
x = 1.0/ns**2
gap = d["e2"] - d["e1"]          # bulk gap, one magnon on the edge triplet
split = np.abs(d["e1"] - d["e0"])  # edge splitting
gap_ref, xi = 0.41048, 6.03        # infinite chain (White and Huse)
# fit gap = a + b/n^2, weighted by n^2, on n >= 48 (next correction is 1/n^3)
nmin = 48
fit = ns >= nmin
gap_inf, slope = np.polynomial.polynomial.polyfit(x[fit], gap[fit], 1,
                                                  w=ns[fit]**2)
# prefactor of exp(-n/xi) on the splittings above 1e-9 (below, the
# sector energies are not accurate enough)
ok = split > 1e-9
pref = np.exp(np.mean(np.log(split[ok]) + ns[ok]/xi))


def draw(theme):
    c = style.colors(theme)
    fig, (a1, a2) = style.figure(2)
    xx = np.linspace(0, 2.0/nmin**2, 50)
    a1.plot(x, gap, "o", color=c["series"][0], label="DMRG, sectors")
    a1.plot(xx, gap_inf + slope*xx, **style.reference(theme),
            label=r"fit $a+b/n^2$ ($n\geq%d$), $a=%.4f$" % (nmin, gap_inf))
    a1.axhline(gap_ref, color=c["ink"], ls=":", lw=style.LW_REF,
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
    a2.plot(np.arange(len(prof)), prof, "o-", color=c["series"][1], ms=3.5)
    a2.axhline(0, color=c["ink2"], lw=0.6)
    a2.set_xlabel(r"site $i$")
    a2.set_ylabel(r"$\langle S^z_i\rangle$ in $S^z_{\rm tot}=1$")
    a2.text(0.5, 0.95, r"$n=%d$" % len(prof), transform=a2.transAxes,
            ha="center", va="top")
    return fig


style.main(draw, "excited-states_haldane")
