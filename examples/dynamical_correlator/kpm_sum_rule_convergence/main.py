# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# KPM dynamical correlator: how the Chebyshev reconstruction converges as
# the bond dimension of the moment recursion (kpmmaxm) grows, checked
# against the exact zeroth-moment sum rule.
#
# For A = Sz on a spin-1/2 site, A*A = 1/4 identically, so
#
#     int dw S_AA(w) = mu_0 = <GS|A A|GS> = 1/4
#
# exactly -- independent of chain length, coupling or bond dimension. That
# makes the integrated weight a parameter-free accuracy measure for the
# whole moment recursion, not just for the dominant peak, and it is the
# invariant tests/test_kpm_dynamical_correlator_python.py asserts and that
# any port of this path (e.g. to GPU, see
# docs/pyitensor_gpu_port_plan.md) has to preserve.
#
# The recursion truncates the MPS after every Chebyshev step
# (pyitensor/chain.py's _kpm_moments_accelerated), so kpmmaxm controls how
# much of the spectral weight survives that truncation -- which is exactly
# what the right-hand panel shows converging.

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain

n = 8                                     # sites
kpmmaxms = [10, 20, 40, 80]               # KPM bond dimensions to compare
es = np.linspace(-1.0, 6.0, 500)          # frequency window for the lineshape
es_wide = np.linspace(-9.0, 9.0, 2500)    # wide enough to hold all the weight
delta = 0.1

def chain(kpmmaxm):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)],
                              itensor_version="python")
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1]
        h = h + sc.Sy[i]*sc.Sy[i+1]
        h = h + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 30
    sc.kpmmaxm = kpmmaxm
    return sc

spectra, errors = [], []
for kpmmaxm in kpmmaxms:
    sc = chain(kpmmaxm)
    name = (sc.Sz[0], sc.Sz[0])
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                       name=name, es=es, delta=delta)
    spectra.append((np.array(x), np.array(y).real))
    # the sum rule needs the whole spectrum, so integrate over the wide window
    xw, yw = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                         name=name, es=es_wide, delta=delta)
    weight = np.trapezoid(np.real(yw), xw)
    errors.append(abs(weight - 0.25))
    print("kpmmaxm = %3d   integrated weight = %.6f   |error| = %.2e"
          % (kpmmaxm, weight, errors[-1]))

# exact reference lineshape from ED, for the same window
sc = chain(kpmmaxms[-1])
xe, ye = sc.get_dynamical_correlator(mode="ED", submode="INV",
                                     name=(sc.Sz[0], sc.Sz[0]),
                                     es=es, delta=delta)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

for kpmmaxm, (x, y) in zip(kpmmaxms, spectra):
    ax1.plot(x, y, label="KPM, kpmmaxm = %d" % kpmmaxm)
ax1.plot(xe, np.array(ye).real, "k--", lw=1.2, label="ED (exact)")
ax1.set_xlabel("Energy $\\omega$")
ax1.set_ylabel("$S_{zz}(\\omega)$")
ax1.set_title("spectral function, %d-site Heisenberg chain" % n)
ax1.legend(fontsize=8)

ax2.plot(kpmmaxms, errors, "o-")
ax2.set_xscale("log", base=2)
ax2.set_yscale("log")
ax2.set_xlabel("kpmmaxm")
ax2.set_ylabel(r"$|\int S_{zz}\,d\omega - 1/4|$")
ax2.set_title("sum-rule error (exact answer is 1/4)")
ax2.grid(alpha=.3)

fig.tight_layout()
plt.show()
