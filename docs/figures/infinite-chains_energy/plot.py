import numpy as np
import style

d = np.load("data.npz")
maxm, r = d["maxm"], d["r"]
err = np.abs(d["energy"] - d["e_exact"])
err_uc2 = abs(d["e_uc2"] - d["e_exact"])
c = np.abs(d["czz"])
# Luttinger form a sqrt(ln r)/r, a fitted in log space on r = 8 to 20
fit = (r >= 8) & (r <= 20)
a = np.exp(np.mean(np.log(c[fit]) - np.log(np.sqrt(np.log(r[fit]))/r[fit])))
rr = np.linspace(2, r.max(), 200)


def draw(theme):
    col = style.colors(theme)
    fig, (ax1, ax2) = style.figure(2)
    ax1.plot(maxm, err, "o-", color=col["series"][0],
             label="iDMRG, one-site cell")
    ax1.plot([maxm[-1]], [err_uc2], "s", color=col["series"][1],
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
