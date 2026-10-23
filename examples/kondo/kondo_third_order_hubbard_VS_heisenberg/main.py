# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Third-order STM/Kondo spectrum (Ternes, arXiv:1505.04430) at T=0 on the
# first site of a half-filled Hubbard chain, against the same spectrum on
# the Heisenberg chain it maps onto at large U. The tip couples to the
# local electron spin S_0 = c^dagger sigma c/2, so the Hubbard spectrum
# carries the charge fluctuations the spin model leaves out.
#
# Every energy is measured in units of the exchange J, and so are the
# bias, the field and the band cutoff omega0, so that curves at different
# U/t lie on one axis:
#   - dimer: J = (sqrt(U^2+16t^2)-U)/2, the exact singlet-triplet gap,
#     so the step sits at eV=J for every U/t and only its height moves
#   - chain: J = 4t^2/U, the leading order of the mapping, so the curves
#     approach the Heisenberg one as U/t grows; the deviation falls as
#     (t/U)^2 once U/t is a few tens, while at U/t below about 30 the
#     higher orders in t/U still matter, so the sweep of the deviation
#     goes up to U/t=256
# The Hubbard charge excitations sit at about U = (U/t)^2 J/4, outside the
# plotted bias window, and inside omega0 = 10 J for U/t < 6.
#
# The chain is solved with ED in the Nf=n sector, which is exact at T=0
# since S_0 conserves the electron number, and one U/t is cross-checked
# against mode="DMRG" (submode="TASK" on native spinful sites, exact here
# since maxm holds the whole sector).

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain, fermionchain

t = 1.0
kw = dict(site=0, Jrho_s=-0.05, T=0., kB=1., order=3)

def hubbard(n, U, b=0., itensor_version="python"):
    """Half-filled Hubbard chain, b*J Zeeman field on every site"""
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n,
            itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h - t*(fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
    h = h + h.get_dagger()
    for i in range(n):
        h = h + U*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5) + b*fc.Sz[i]
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=n) # half filling
    fc.maxm = 64 ; fc.nsweeps = 20 ; fc.cutoff = 0.
    return fc

def heisenberg(n, J, b=0.):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(n): h = h + b*sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc

def spectrum(chain, x, J, **extra):
    """dI/dV/2pi against eV/J, with omega0 and Gamma0 in units of J"""
    _, d = chain.get_kondo_spectrum(x*J, omega0=10*J, Gamma0=1e-2*J,
                                    **kw, **extra)
    return d/(2*np.pi)

x = np.linspace(-3., 3., 240) # eV/J, off the steps at eV=J of the dimer
Uts = [4., 8., 16., 32.] # the curves that are plotted
Uts_dev = Uts + [64., 128., 256.] # the deviation against t/U

# dimer, exact J
dimer = {}
for Ut in Uts:
    U = Ut*t
    J = (np.sqrt(U**2+16*t**2)-U)/2
    dimer[Ut] = spectrum(hubbard(2, U), x, J)
dimer_heis = spectrum(heisenberg(2, 1.), x, 1.)

# six-site chain, J = 4t^2/U, without and with a field b=J
n = 6
chain = {0.: {}, 1.: {}}
heis = {}
for b in chain:
    for Ut in Uts_dev:
        J = 4*t**2/Ut
        chain[b][Ut] = spectrum(hubbard(n, Ut*t, b=b*J), x, J)
    heis[b] = spectrum(heisenberg(n, 1., b=b), x, 1.)

# DMRG cross-check, one U/t (cutoff=0, maxm=64 holds the Nf=6 sector)
Ut = 8. ; J = 4*t**2/Ut
fc = hubbard(n, Ut*t, b=J)
fc.set_conserved_sector() # the TASK route solves the whole Fock space
_, d = fc.get_kondo_spectrum(x[::8]*J, omega0=10*J, Gamma0=1e-2*J,
                             mode="DMRG", submode="TASK", nkry=200, **kw)
err = np.max(np.abs(d/(2*np.pi)-chain[1.][Ut][::8]))
print("TASK against ED, U/t=8, b=J: max |difference| =", err)
assert err < 1e-6

# deviation from the Heisenberg curve against t/U, averaged over the bias
# window: at T=0 the steps are sharp, so a maximum would only measure
# whether a shifted step crossed a grid point
dev = {b: [np.mean(np.abs(chain[b][Ut]-heis[b])) for Ut in Uts_dev]
       for b in chain}
for b in chain:
    for Ut, dv in zip(Uts_dev, dev[b]):
        print("b=%.0f J, U/t=%4.0f: mean |Hubbard - Heisenberg| = %.5f" % (b, Ut, dv))

fig, axs = plt.subplots(2, 2, figsize=(10, 8))
colors = plt.cm.viridis(np.linspace(0., 0.85, len(Uts)))
panels = [(axs[0,0], dimer, dimer_heis, "dimer, B=0, J exact"),
          (axs[0,1], chain[0.], heis[0.], "n=6, B=0, J=4t$^2$/U"),
          (axs[1,0], chain[1.], heis[1.], "n=6, B=J, J=4t$^2$/U")]
for ax, curves, ref, title in panels:
    for c, Ut in zip(colors, Uts):
        ax.plot(x, curves[Ut], c=c, label="Hubbard U/t=%g" % Ut)
    ax.plot(x, ref, "k--", label="Heisenberg")
    ax.set_xlabel("eV / J") ; ax.set_ylabel(r"dI/dV  [$e^2T_0^2/\hbar$]")
    ax.set_title(title)
axs[0,0].legend(fontsize=8)
ax = axs[1,1]
tU = t/np.array(Uts_dev)
for b, m in zip(chain, ("o", "s")):
    ax.loglog(tU, dev[b], m+"-", label="n=6, B=%g J" % b)
# anchored at the largest U/t, where the deviation is asymptotic
ax.loglog(tU, dev[1.][-1]*(tU/tU[-1])**2, "k:", label=r"$\propto (t/U)^2$")
ax.set_xlabel("t / U") ; ax.set_ylabel("mean |Hubbard - Heisenberg|")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig("kondo_third_order_hubbard_VS_heisenberg.png", dpi=150)
plt.show()
