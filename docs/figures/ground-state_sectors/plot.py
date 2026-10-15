import numpy as np
import style

d = np.load("data.npz")
n, M = d["n"], d["M"]
HMAX = 2.4


def staircase(e):
    """m(h) as a step curve: M for h between h_{M-1} and h_M."""
    x = np.concatenate([[0.0], np.diff(e), [HMAX]])  # h_M = E0(M+1) - E0(M)
    y = np.concatenate([M, [M[-1]]])/n
    return x, y


def draw(theme):
    c = style.colors(theme)
    fig, ax = style.figure(1)
    ax.step(*staircase(d["e_xx"]), where="post", color=c["series"][0],
            label="XX, DMRG")
    ax.step(*staircase(d["e_heis"]), where="post", color=c["series"][1],
            label="Heisenberg, DMRG")
    h = np.linspace(0, 1, 400)
    ax.plot(np.append(h, HMAX), np.append(np.arcsin(h)/np.pi, 0.5),
            **style.reference(theme), label=r"XX, exact $\arcsin(h/J)/\pi$")
    ax.axvline(2.0, color=c["ink"], ls=":", lw=style.LW_REF)
    ax.text(2.03, 0.05, r"$h_s=2J$", color=c["ink"])
    ax.set_xlim(0, HMAX)
    ax.set_ylim(0, 0.53)
    ax.set_xlabel(r"field $h$ [$J$]")
    ax.set_ylabel(r"magnetization per site $M/n$")
    ax.legend(loc="upper left")
    return fig


style.main(draw, "ground-state_sectors")
