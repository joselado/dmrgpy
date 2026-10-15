import numpy as np
import matplotlib.ticker as ticker
import style

d = np.load("data.npz")
r, czz_half, czz_one = d["r"], np.abs(d["czz_half"]), np.abs(d["czz_one"])
xi = 6.03  # Haldane correlation length
win = (r >= 8) & (r <= 32)  # fit window of the amplitudes a and b


def fit(c, f):
    """f times the amplitude fitted in log space on the window."""
    return np.exp(np.mean(np.log(c[win]) - np.log(f[win])))*f


ref_half = fit(czz_half, np.sqrt(np.log(r))/r)  # Luttinger liquid, log correction
ref_one = fit(czz_one, np.exp(-r/xi)/np.sqrt(r))  # Ornstein-Zernike, gapped
ok = r >= 2  # sqrt(ln r)/r vanishes at r=1
# xi fitted freely with the Ornstein-Zernike form, quoted in the guide
xi_free = -1/np.polyfit(r[win], np.log(czz_one[win]*np.sqrt(r[win])), 1)[0]


def draw(theme):
    c = style.colors(theme)
    fig, ax = style.figure(1)
    ax.plot(r, czz_half, "o-", color=c["series"][0], ms=4, label=r"$S=1/2$")
    ax.plot(r, czz_one, "s-", color=c["series"][1], ms=4, label=r"$S=1$")
    ax.plot(r[ok], ref_half[ok], **style.reference(theme),
            label=r"$a\sqrt{\ln r}/r$")
    ax.plot(r, ref_one, color=c["ink"], ls=":", lw=style.LW_REF,
            label=r"$b\,e^{-r/\xi}/\sqrt{r}$, $\xi=%.2f$" % xi)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([1, 2, 4, 8, 16, 32])
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel(r"distance $r$ [sites]")
    ax.set_ylabel(r"$|\langle S^z_iS^z_{i+r}\rangle|$")
    ax.legend()
    return fig


style.main(draw, "ground-state_correlators")
