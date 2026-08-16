# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Cost of the two workhorse dynamical-correlator submodes -- "KPM" (the
# Chebyshev-moment recursion, kpmdmrg.py) and "TD" (real-time TDVP plus a
# windowed Fourier transform, timedependent.py) -- as a function of chain
# length, on both always-available DMRG backends, together with the
# spectral functions themselves so the timings can be read next to the
# result they buy.
#
# This is the script behind docs/documentation.md's section 4.8a ("Cost of
# the dynamical-correlator submodes, and what makes them fast"), whose
# reference point is L=30; it sweeps a range of smaller sizes so it stays
# runnable by hand in a couple of minutes.
#
# NOTE ON TIMING: pin BLAS to one thread before running this, e.g.
#     MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 main.py
# DMRG is a great many *small* dense linear-algebra calls, and under
# thread oversubscription the numbers below stop being comparable at all
# -- see CLAUDE.md and src/dmrgpy/blasthreads.py.

import time
import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import spinchain

Ls = [8, 12, 16, 20]
SUBMODES = ["KPM", "TD"]
BACKENDS = [3, "python"]
es = np.linspace(-1.0, 6.0, 300)
DELTA = 0.3


def heisenberg(L, backend):
    """A uniform S=1/2 Heisenberg chain on the requested backend. Written
    the textbook way (Sx.Sx + Sy.Sy + Sz.Sz) on purpose: that is exactly
    the form mpscpp3/mo_terms.h's Sy -> (S+ - S-)/(2i) rewrite turns into
    a real-valued MPO, and hence the form the timings below are about."""
    sc = spinchain.Spin_Chain(["S=1/2"] * L)
    if backend == "python":
        sc.setup_python()
    h = 0
    for i in range(L - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc


times = {(b, s): [] for b in BACKENDS for s in SUBMODES}
spectra = {}  # (backend, submode) -> the L=max(Ls) spectral function

for L in Ls:
    for backend in BACKENDS:
        sc = heisenberg(L, backend)
        sc.gs_energy()  # ground state once, shared by both submodes
        name = (sc.Sz[L // 2], sc.Sz[L // 2])
        for submode in SUBMODES:
            t0 = time.time()
            _es, y = sc.get_dynamical_correlator(name=name, submode=submode,
                                                  es=es, delta=DELTA)
            dt = time.time() - t0
            times[(backend, submode)].append(dt)
            if L == Ls[-1]:
                spectra[(backend, submode)] = y
            print("L=%2d  %-8s %-4s  %7.2f s" % (L, str(backend), submode, dt),
                  flush=True)

fig, (ax_t, ax_s) = plt.subplots(1, 2, figsize=(11, 4.2))

styles = {(3, "KPM"): ("o-", "C0"), (3, "TD"): ("s--", "C1"),
          ("python", "KPM"): ("o-", "C2"), ("python", "TD"): ("s--", "C3")}
for key, ts in times.items():
    marker, color = styles[key]
    ax_t.plot(Ls, ts, marker, c=color,
              label="itensor_version=%s, %s" % (repr(key[0]), key[1]))
ax_t.set_xlabel("chain length L")
ax_t.set_ylabel("wall time (s)")
ax_t.set_yscale("log")
ax_t.set_title("cost of one dynamical correlator")
ax_t.legend(fontsize=8)

for key, y in spectra.items():
    marker, color = styles[key]
    ax_s.plot(es, -y.imag if key[1] == "TD" else y.real, c=color,
              ls="--" if key[1] == "TD" else "-",
              label="itensor_version=%s, %s" % (repr(key[0]), key[1]))
ax_s.set_xlabel(r"$\omega$")
ax_s.set_ylabel(r"$S^{zz}(\omega)$  (L=%d, middle site)" % Ls[-1])
ax_s.set_title("what those timings buy")
ax_s.legend(fontsize=8)

plt.tight_layout()
plt.show()
