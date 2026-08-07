# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Demonstrates KPM energy truncation (Holzner, Weichselbaum, McCulloch &
# von Delft, "Chebyshev matrix product state approach for spectral
# functions", PRB 83, 195115 (2011), Sec. III-B; pyitensor-only, see
# src/dmrgpy/pyitensor/kpm_energy_truncation.py).
#
# By default, kpmdmrg.py's KPM rescales H into a window Ws about as wide
# as the full many-body bandwidth W (kpm_scale~0.7, a safety margin, not
# a resolution knob) -- this always works, but wastes the KPM expansion's
# resolution on parts of the spectrum a given correlator has no weight
# in. Choosing a narrower window (smaller kpm_scale) concentrates that
# same expansion order onto the physically relevant region instead -- but
# without a safeguard, numerical noise in the Chebyshev recursion lets
# "leaked" high-energy components blow up exponentially (Chebyshev
# polynomials are unbounded outside [-1,1]), which chain.py's own
# _check_kpm_moment only detects (raises), never fixes. Energy truncation
# is that fix: after every Chebyshev vector is formed, it projects out
# any locally-high-energy component (in a small per-site Krylov
# subspace) before the vector is used or fed into the next recursion
# step.
#
# Enabling kpm_energy_truncate also switches the rescaling itself from
# _scaled_hamiltonian() (centered on the full bandwidth's midpoint) to
# _scaled_hamiltonian_gs_anchored() (anchored at the ground state, per
# the paper's own Eq. 21b) -- see that method's docstring for why: a
# dynamical correlator's Chebyshev vectors are built by acting operators
# on the ground state, so the relevant window sits just above it, not
# around the spectrum's geometric middle.
import numpy as np
from dmrgpy import spinchain

n = 4
delta = 0.05
es = np.linspace(0.3, 1.0, 71)
gap = 0.658919  # exact E1-E0 gap for this chain (see test_excited_states.py)


def make_chain():
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    sc.setup_python()  # energy truncation is pyitensor-only
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc


name = lambda sc: (sc.Sz[0], sc.Sz[0])  # <Sz_0(t) Sz_0(0)>

# --- exact ED ground truth ---
sc_ed = make_chain()
x_ed, y_ed = sc_ed.get_dynamical_correlator(mode="ED", submode="INV", name=name(sc_ed),
                                             es=es, delta=delta)
y_ed = np.array(y_ed).real
peak_ed = x_ed[np.argmax(y_ed)]
print("ED peak: %.3f (exact gap: %.6f)" % (peak_ed, gap))
assert abs(peak_ed - gap) < 0.02

# --- default (safe) KPM window: works fine, no truncation needed ---
sc_default = make_chain()
x_def, y_def = sc_default.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name(sc_default),
                                                     es=es, delta=delta)
y_def = np.array(y_def).real
peak_def = x_def[np.argmax(y_def)]
print("KPM (default kpm_scale=%.2f) peak: %.3f" % (sc_default.kpm_scale, peak_def))
assert abs(peak_def - gap) < 0.05

# --- narrowed window, energy truncation OFF: reproduces the known failure ---
kpm_scale = 0.65
sc_narrow = make_chain()
sc_narrow.kpm_scale = kpm_scale
sc_narrow.kpm_energy_truncate = True   # switches to the gs-anchored window...
sc_narrow.kpm_truncate_nsweeps = 0     # ...but disable the truncation sweeps themselves
try:
    sc_narrow.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name(sc_narrow),
                                        es=es, delta=delta)
    raise SystemExit("expected a RuntimeError here -- did chain.py's KPM change?")
except RuntimeError as e:
    print("kpm_scale=%.2f without energy truncation: raised as expected (%s)" % (kpm_scale, e))

# --- same narrow window, energy truncation ON: computes cleanly, still ---
# --- finds the right peak ---
sc_trunc = make_chain()
sc_trunc.kpm_scale = kpm_scale
sc_trunc.kpm_energy_truncate = True
sc_trunc.kpm_truncate_dK = 30       # Table I's recommended values
sc_trunc.kpm_truncate_nsweeps = 10
x_tr, y_tr = sc_trunc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name(sc_trunc),
                                                es=es, delta=delta)
y_tr = np.array(y_tr).real
peak_tr = x_tr[np.argmax(y_tr)]
print("KPM (kpm_scale=%.2f, energy truncation ON) peak: %.3f" % (kpm_scale, peak_tr))
assert abs(peak_tr - gap) < 0.05

print("OK: energy truncation lets a narrower KPM window compute safely "
      "and still locate the exact excitation gap.")

# --- plot: ED ground truth vs default-window KPM vs narrow-window KPM
# with energy truncation, overlaid on the same frequency axis ---
import matplotlib.pyplot as plt

fig = plt.figure()
plt.plot(x_ed, y_ed, c="black", lw=2, label="ED (exact)")
plt.plot(x_def, y_def, c="blue", ls="--",
         label="KPM, default window (kpm_scale=%.2f)" % sc_default.kpm_scale)
plt.plot(x_tr, y_tr, c="red", ls="-.",
         label="KPM, narrow window + energy truncation (kpm_scale=%.2f)" % kpm_scale)
plt.axvline(gap, c="gray", ls=":", label="exact gap")
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.legend(fontsize=8)
plt.title("KPM energy truncation: narrow-window accuracy vs ED")
plt.savefig("kpm_energy_truncation.png", dpi=150)
print("\nPlot saved to kpm_energy_truncation.png")
plt.show()
