# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Conserved-sector targeting on the two backends that have it:
# itensor_version=3 (block-sparse QN tensors) and itensor_version="python"
# (dense storage plus a charge penalty). See docs/user_guide.md's
# "Targeting a quantum-number sector" and pyitensor/sector.py.
#
# set_conserved_sector(Nf=k) confines the whole calculation -- starting
# state, every sweep, the returned wavefunction -- to exactly k particles,
# which is the grand-canonical-free route to the addition spectrum E(N).
# The two backends get there by different means and must agree anyway;
# this script checks that they do, against ED restricted to the same
# sector, and then shows the one thing that is *not* shared.
#
# That one thing is the third panel. ITensor v3 cannot leave a sector: an
# amplitude outside it has nowhere to be stored. The pure-Python backend
# stores everything densely and uses the quantum numbers only as labels,
# so a truncating SVD leaks ~1e-16 of amplitude into neighbouring sectors
# on every bond (LAPACK's Householder bidiagonalization mixes rows across
# charge blocks), and a variational sweep amplifies that leak toward
# whichever sector is lower in energy. A charge penalty
# lambda*(N-k)^2 -- exactly zero on the target sector, so it changes no
# reported number -- is what keeps the sweep confined. Panel 3 sweeps
# lambda down to zero on a chain with *attractive* interactions, where the
# global ground state is the full band and the requested Nf=2 sector is as
# far from it as it gets: below a threshold set by the Hamiltonian's own
# energy scale, the solve simply falls out of the sector, and dmrgpy
# raises rather than reporting the wrong sector's energy.
#
# The model is a spinless t-V chain (nearest-neighbour hopping plus a
# nearest-neighbour interaction V).
import time

import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import cppext, fermionchain
from dmrgpy.multioperator import MO2matrix


def tV_chain(n, v, backend):
    """A t-V chain on the requested backend, plus its Hamiltonian and
    total-number operator."""
    fc = fermionchain.Fermionic_Chain(n)
    if backend == "python":
        fc.setup_python()
    else:
        fc.setup_cpp(version=backend)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + v * fc.N[i] * fc.N[i + 1]
    fc.maxm = 40
    fc.nsweeps = 10
    fc.set_hamiltonian(h)
    return fc, h, sum(fc.N)


def ed_sector_energy(chain, h, charge, target):
    """Lowest eigenvalue of h restricted to the sector where the (diagonal)
    operator `charge` equals `target` -- a submatrix in the ED product
    basis, not a filter on full-spectrum eigenvectors (those come back as
    arbitrary superpositions whenever a level is degenerate across
    sectors)."""
    ed = chain.get_ED_obj()
    H = np.array(MO2matrix(h, ed).todense())
    N = np.array(MO2matrix(charge, ed).todense())
    sel = np.abs(np.diag(N).real - target) < 1e-9
    return float(np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0])


# ---------------------------------------------------------------- E(N)
n, v = 10, 2.0
fillings = list(range(n + 1))
backends = ["python"] + ([3] if cppext.available(3) else [])
if 3 not in backends:
    print("mpscpp3 is not compiled here -- running the pure-Python backend alone")

energies, times = {}, {}
for backend in backends:
    fc, h, number = tV_chain(n, v, backend)
    es, t0 = [], time.time()
    for nf in fillings:
        fc.set_conserved_sector(Nf=nf)
        es.append(fc.gs_energy())
        # the state really is in the sector, not merely close in energy
        assert abs(fc.vev(number).real - nf) < 1e-6
    energies[backend] = np.array(es)
    times[backend] = time.time() - t0
    print("%-8s addition spectrum in %.1f s" % (str(backend), times[backend]))

fc_ed, h_ed, number_ed = tV_chain(n, v, "python")
e_ed = np.array([ed_sector_energy(fc_ed, h_ed, number_ed, nf) for nf in fillings])

for backend in backends:
    err = np.max(np.abs(energies[backend] - e_ed))
    print("%-8s max |E - E_ED| over all sectors: %.2e" % (str(backend), err))
    assert err < 1e-6

# ------------------------------------------------- the charge penalty
# An attractive interaction puts the global ground state at full filling,
# so Nf=2 is a sector the sweep actively wants to leave.
n_pen, nf_pen = 12, 2
fc_pen, h_pen, number_pen = tV_chain(n_pen, -1.0, "python")
fc_pen.nsweeps = 20
e_pen_ed = ed_sector_energy(fc_pen, h_pen, number_pen, nf_pen)
default_scale = sum(abs(complex(c)) for c, _ in h_pen.to_terms())

lambdas = np.array([0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]) * default_scale
errors, confined = [], []
for lam in lambdas:
    fc_pen.restart()
    fc_pen.set_conserved_sector(Nf=nf_pen)
    fc_pen.set_sector_penalty(lam)
    try:
        e = fc_pen.gs_energy()
        errors.append(abs(e - e_pen_ed))
        confined.append(True)
    except RuntimeError as err:  # the solve left the sector, and says so
        errors.append(np.nan)
        confined.append(False)
    print("lambda = %8.2f  ->  %s" % (lam, "in sector" if confined[-1] else "LEFT the sector"))
errors = np.array(errors)
assert confined[-1] and not confined[0], "the penalty must be what confines the sweep"

# ------------------------------------------------------------- plots
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

ax[0].plot(fillings, e_ed, "k-", lw=3, alpha=0.35, label="ED (sector-restricted)")
for backend, style in zip(backends, ["o--", "s:"]):
    ax[0].plot(fillings, energies[backend], style, label="DMRG, itensor_version=%s"
               % repr(backend))
ax[0].set_xlabel("particle number $N$")
ax[0].set_ylabel("$E_0(N)$")
ax[0].set_title("addition spectrum, $t$-$V$ chain ($n=%d$, $V=%.1f$)" % (n, v))
ax[0].legend()

for backend, style in zip(backends, ["o-", "s-"]):
    ax[1].semilogy(fillings, np.maximum(np.abs(energies[backend] - e_ed), 1e-16),
                   style, label="itensor_version=%s" % repr(backend))
ax[1].axhline(1e-6, color="k", ls=":", lw=1, label="test tolerance")
ax[1].set_xlabel("particle number $N$")
ax[1].set_ylabel(r"$|E_\mathrm{DMRG}-E_\mathrm{ED}|$")
ax[1].set_title("both backends against sector-restricted ED")
ax[1].legend()

left = np.array(lambdas)[~np.array(confined)]
ax[2].semilogy(lambdas, np.where(np.isnan(errors), 1e0, errors), "o-",
               label=r"$|E-E_\mathrm{ED}|$ (in sector)")
for lam in left:
    ax[2].axvline(lam, color="crimson", lw=6, alpha=0.25)
ax[2].plot([], [], color="crimson", lw=6, alpha=0.25, label="solve left the sector")
ax[2].axvline(default_scale, color="k", ls="--", lw=1,
              label=r"dmrgpy's default $\lambda$")
ax[2].set_xlabel(r"charge-penalty strength $\lambda$")
ax[2].set_ylabel(r"$|E-E_\mathrm{ED}|$")
ax[2].set_title(r'itensor_version="python", $n=%d$, $N_f=%d$, attractive $V$'
                % (n_pen, nf_pen))
ax[2].legend()

plt.tight_layout()
plt.show()
