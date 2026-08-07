# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

# Compare submode="TD"'s spectral tail against KPM, and show how the two
# tools this feature adds -- a selectable damping/window taper
# (timedependent.py::_damping_window) and linear-prediction extrapolation
# (dynamicstk/linearprediction.py) -- sharpen it. Empirically,
# damping="exp"+predict=True gave the narrowest, best-centered peak among
# every combination tried (see docs/td_dynamical_correlator_sharpening_plan.md),
# so that combination -- damping="exp" (unchanged) and predict=True -- is
# now dynamical_correlator's own default; pass predict=False to recover
# the old, pre-this-feature behavior exactly.

def heisenberg_chain(n):
    spins = ["S=1/2" for i in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 30
    return sc

# Panel 1: a small chain (4 sites) whose lowest excitation is a single,
# well-isolated resonance -- deliberately at a broad delta, so the *shape*
# differences between damping choices and linear prediction (all designed
# to sharpen a broadening-limited peak) are directly visible. A narrow
# delta already resolves the peak well on its own, which would hide the
# effect this feature is about.
sc4 = heisenberg_chain(4)
name4 = (sc4.Sz[0],sc4.Sz[0])
delta_peak = 1.5e-1
es_peak = np.linspace(0.3,1.0,600)

peak_results = {}
peak_results["TD old default (exp, no predict)"] = sc4.get_dynamical_correlator(
        es=es_peak,name=name4,submode="TD",delta=delta_peak,
        damping="exp",predict=False)
peak_results["TD new default (exp + predict)"] = sc4.get_dynamical_correlator(
        es=es_peak,name=name4,submode="TD",delta=delta_peak)
peak_results["TD gaussian + predict"] = sc4.get_dynamical_correlator(
        es=es_peak,name=name4,submode="TD",delta=delta_peak,damping="gaussian")
peak_results["TD parzen + predict"] = sc4.get_dynamical_correlator(
        es=es_peak,name=name4,submode="TD",delta=delta_peak,damping="parzen")

# Panel 2: a longer chain (6 sites, several excitations) at "TD"'s own
# default delta, comparing its far tail against KPM -- this is the "long
# tail" complaint this feature addresses: TD's tail (even at the new
# default) is heavy (Lorentzian-derived) compared to KPM's Jackson-kernel
# reconstruction, which is designed specifically to suppress it. Neither
# `damping` nor `predict` change the *algebraic* nature of the tail far
# from any resonance (they reshape the peak/near-tail, not the
# asymptotic decay), so this gap is expected to remain.
sc6 = heisenberg_chain(6)
name6 = (sc6.Sz[0],sc6.Sz[0])
delta_tail = 5e-2
es_wide = np.linspace(-1.0,6.0,3000)

tail_results = {}
tail_results["TD (new default: exp + predict)"] = sc6.get_dynamical_correlator(
        es=es_wide,name=name6,submode="TD",delta=delta_tail)
tail_results["KPM (Jackson kernel)"] = sc6.get_dynamical_correlator(
        es=es_wide,name=name6,submode="KPM",delta=delta_tail)

# plot the results
import matplotlib.pyplot as plt

fig,axs = plt.subplots(1,2,figsize=(11,4.5))

ax = axs[0]
for label,(x,y) in peak_results.items():
    ax.plot(x,np.abs(y),label=label)
ax.set_xlabel("frequency [J]")
ax.set_ylabel("|S(omega)|")
ax.set_title("Peak shape: old vs new default vs other choices\n(delta=%.2f)"%delta_peak)
ax.legend(fontsize=8)

ax = axs[1]
for label,(x,y) in tail_results.items():
    ax.semilogy(x,np.abs(y)+1e-12,label=label)
ax.set_xlabel("frequency [J]")
ax.set_ylabel("|S(omega)|  (log scale)")
ax.set_xlim([-1.0,6.0])
ax.set_title("Tail away from resonance: TD vs KPM\n(delta=%.2f)"%delta_tail)
ax.legend(fontsize=8)

fig.tight_layout()
plt.savefig("dynamical_correlator_td_sharpen.png",dpi=150)
plt.show()
