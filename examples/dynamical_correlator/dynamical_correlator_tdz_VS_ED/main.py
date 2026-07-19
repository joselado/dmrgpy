# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test for submode="TDZ" (tdz.py): complex-time evolution +
# perturbative real-axis reconstruction, following Cao, Lu, Stoudenmire &
# Parcollet, "Dynamical correlation functions from complex time
# evolution", arXiv:2311.10909.
#
# Checked two ways:
#
#   - Peak position vs exact ED (like KPM in dynamical_correlator_VS_ED):
#     TDZ reconstructs the real-time correlator only up to the requested
#     Taylor order in alpha0 plus the usual dt time-discretization error,
#     so its dominant peak should land close to (not necessarily exactly
#     on) the exact excitation energy -- checked on all three DMRG
#     backends (2, 3, "python"). itensor_version=2 has no TDVP, so it
#     exercises tdz.py's MPO-Taylor fallback path (evolve_taylor_step())
#     rather than TDVP (tdvp_step()); both must agree with each other and
#     with ED.
#
#   - Convergence of the alpha0-Taylor reconstruction itself: at a fixed,
#     deliberately large alpha0, the raw complex-time correlator (n_max=0,
#     no reconstruction at all) is compared to the reconstructed one
#     (n_max=4) against the real-time ("TD" submode) reference, directly
#     in the time domain (not just peak positions) -- confirming that
#     each order of the Taylor-in-alpha0 reconstruction (Appendix B of
#     the paper) actually improves the agreement, not just that the final
#     answer happens to look reasonable. This is the check that caught a
#     real sign/prefactor bug during development (an erroneous extra
#     factor of "i" in assembling phi^(0)+sum(g^(n)) that made n_max=1
#     *worse* than n_max=0): with the bug, increasing n_max did not
#     monotonically improve agreement; with the fix, each order improves
#     it by roughly a factor of alpha0.
import numpy as np
from scipy.signal import find_peaks
from dmrgpy import spinchain
from dmrgpy import tdz as tdzmod
from dmrgpy import timedependent as tdmod

n = 4
es = np.linspace(-0.5, 4.0, 200)  # fine enough to resolve peak positions
delta = 0.15


def make_chain(itensor_version=None):
    sc = spinchain.Spin_Chain([2 for i in range(n)])
    if itensor_version is not None:
        if itensor_version != "python": sc.setup_cpp(itensor_version)
        else: sc.setup_python()
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc


def peak_energies(x, y):
    """Energies of the resolvable peaks (>=5% of the tallest one)"""
    y = np.array(y).real
    x = np.array(x)
    idx, _ = find_peaks(y, height=0.05*np.max(y))
    return x[idx]


# --- exact ED ground truth ---
sc_ed = make_chain()
name = (sc_ed.Sz[0], sc_ed.Sz[0])  # <Sz_0(t) Sz_0(0)>
x_ed, y_ed = sc_ed.get_dynamical_correlator(mode="ED", submode="ED", name=name, es=es, delta=delta)
y_ed = np.array(y_ed)
peaks_ed = peak_energies(x_ed, y_ed)
print("ED peak energies =", peaks_ed)
assert len(peaks_ed) >= 2, "expected at least 2 resolvable peaks in this window"

# --- TDZ on all three DMRG backends ---
tol_peak = 0.15  # generous: TDZ's default dt=0.1 time-discretization + n_max=4 truncation
for v in [2, 3, "python"]:
    sc = make_chain(v)
    name = (sc.Sz[0], sc.Sz[0])
    x_tdz, y_tdz = sc.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=name,
            es=es, alpha0=0.1, n_max=4, dt=0.05, delta=delta)
    peaks_tdz = peak_energies(x_tdz, y_tdz)
    print("TDZ(v%s) peak energies =" % v, peaks_tdz)
    assert len(peaks_tdz) >= len(peaks_ed), \
        "v%s TDZ found %d peaks, expected at least the %d ED found" % (v, len(peaks_tdz), len(peaks_ed))
    for e_ed in peaks_ed:
        closest = peaks_tdz[np.argmin(np.abs(peaks_tdz-e_ed))]
        diff = abs(closest-e_ed)
        print("  v%s ED peak at %.4f -> nearest TDZ peak at %.4f (diff=%.3f)" % (v, e_ed, closest, diff))
        assert diff < tol_peak, \
            "v%s TDZ has no peak within %g of the ED peak at %.4f (closest: %.4f)" % (v, tol_peak, e_ed, closest)

# --- alpha0-Taylor reconstruction: each order should improve agreement
# with the real-time ("TD") reference at a fixed, deliberately large
# alpha0, confirming the reconstruction itself (not just the final
# peak-matching result) is implemented correctly ---
sc = make_chain("python")
sc.get_gs()
name = (sc.Sz[0], sc.Sz[0])
A = name[0].get_dagger()
B = name[1]
dt, nt = 0.05, 40
omega0 = 2.*np.pi/(nt*dt)

ts_ref, cs_ref = tdmod.evolution_dmrg_DC(sc, name=name, nt=nt, dt=dt)
cs_ref_raw = cs_ref.real - 1j*cs_ref.imag  # undo evolution_dmrg_DC's own conjugation

alpha0 = 0.2
diffs = []
for n_max in range(5):
    _, cs_z = tdzmod._complex_time_correlator(sc, A, B, alpha0, n_max, dt, nt, omega0)
    # time-aligned comparison: "TD"'s correlator[k] is measured after the
    # (k+1)-th step (true time (k+1)*dt) but returned under label k*dt --
    # see tdz.py's own note on this pre-existing timedependent.py quirk.
    diff = np.max(np.abs(cs_z[1:] - cs_ref_raw[:len(cs_z)-1]))
    diffs.append(diff)
    print("n_max=%d  max|TDZ - TD| (time-aligned) = %.3e" % (n_max, diff))

for k in range(1, len(diffs)):
    assert diffs[k] < diffs[k-1], \
        "increasing n_max from %d to %d did not improve agreement with the real-time reference (%.3e -> %.3e) -- the alpha0-Taylor reconstruction is broken" % (k-1, k, diffs[k-1], diffs[k])
assert diffs[-1] < 1e-6, "n_max=4 reconstruction should match the real-time reference very closely at alpha0=0.2, got %.3e" % diffs[-1]

print("TEST PASSED")
