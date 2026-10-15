import numpy as np
import style

d = np.load("data.npz")


def wrap(g):
    """A phase in units of pi, mapped to [-1/2, 3/2) so that 0 and 2 meet."""
    return (np.asarray(g) + 0.5) % 2 - 0.5


def draw(theme):
    c = style.colors(theme)
    fig, (a, b, e) = style.figure(3)
    # at delta=0 the gap closes at k=pi, which is on the k grid, and the
    # phase there is undefined: leave that point out
    ok = np.abs(d["delta_ssh"]) > 1e-9
    a.plot(d["delta_ssh"][ok], wrap(d["zak"][ok]), **style.reference(theme),
           label="SSH Zak phase")
    a.plot(d["delta"], wrap(d["gamma"]), "o", color=c["series"][0],
           label="spin ring, DMRG")
    a.set_xlabel(r"dimerization $\delta$")
    a.set_ylabel(r"Berry phase $\gamma/\pi$")
    a.set_yticks([0, 0.5, 1])
    a.set_ylim(-0.15, 1.25)
    a.legend(loc="upper right")

    bc = float(d["b_c"])
    b.plot(d["bz"], d["sz_ed"], **style.reference(theme), label="ED")
    b.plot(d["bz"], d["sz"], "o", ms=3.5, color=c["series"][1], label="DMRG")
    b.set_xlabel(r"impurity field $B$ [$t$]")
    b.set_ylabel(r"total $\langle S^z\rangle$")
    b.legend(loc="lower left")

    e.plot(d["bz"], d["parity_ed"], **style.reference(theme), label="ED")
    e.plot(d["bz"], d["parity"], "o", ms=3.5, color=c["series"][2],
           label="DMRG")
    e.set_xlabel(r"impurity field $B$ [$t$]")
    e.set_ylabel(r"fermion parity $P$")
    e.set_yticks([-1, 0, 1])
    e.legend(loc="lower left")
    for ax in (b, e):
        ax.axvline(bc, color=c["ink2"], lw=0.8, ls=":")
    return fig


style.main(draw, "topology_berry_phase")
