import numpy as np
import style

d = np.load("data.npz")
maxm = d["maxm"]
err = np.abs(d["energy"] - d["e_ed"])/d["n"]  # abs: at the floor it can be < 0
de = d["de"]/d["n"]
# err = c dE^2, c fitted on the points above the numerical floor
fit = err > 5e-14
c_fit = np.exp(np.mean(np.log(err[fit]) - np.log(de[fit]**2)))


def draw(theme):
    c = style.colors(theme)
    fig, ax = style.figure(1)
    ax.plot(maxm, err, "o-", color=c["series"][0],
            label=r"energy error $(E_0-E_0^{\rm ED})/n$")
    ax.plot(maxm, de, "s-", color=c["series"][1],
            label=r"energy fluctuation $\delta E/n$")
    ax.plot(maxm[fit], c_fit*de[fit]**2, **style.reference(theme),
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
