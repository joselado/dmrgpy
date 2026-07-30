# Timing comparison: ITensor v3's native KPM (mpscpp3/chain_session.h)
# with and without energy truncation (Holzner, Weichselbaum, McCulloch &
# von Delft, "Chebyshev matrix product state approach for spectral
# functions", PRB 83, 195115 (2011), Sec. III-B).
#
# Energy truncation (kpm_energy_truncate=True) lets the KPM window be
# narrower than the full many-body bandwidth (see
# examples/dynamical_correlator/dynamical_correlator_kpm_energy_truncation
# for the accuracy-focused demo, and
# src/dmrgpy/mpscpp3/chain_session.h's kpm_energy_truncate()/
# scaled_hamiltonian_gs_anchored() for the implementation) -- but the
# extra truncation sweeps it runs (kpm_truncate_dK Krylov vectors x
# kpm_truncate_nsweeps sweeps, once per Chebyshev moment) are themselves
# real cost. This script turns what was previously an ad hoc, unsaved
# measurement into a saved, rerunnable timing regression: it is a net
# *slowdown* on small test systems like the ones below (the narrower
# window's resolution gain doesn't pay for itself until a system is
# large enough for a correlator's spectral weight to genuinely separate
# from the full bandwidth) -- this script exists to quantify that
# slowdown concretely, not to claim truncation is faster here.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time
import numpy as np
from dmrgpy import spinchain

DELTA = 0.08
ES = np.linspace(0.3, 1.2, 41)
NARROW_KPM_SCALE = 0.65  # matches tests/test_kpm_energy_truncation_v3_accuracy.py
TRUNC_DK = 30       # paper's Table I recommended values
TRUNC_NSWEEPS = 10


def make_chain(n):
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    sc.setup_cpp(3)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc


def run_untruncated(n):
    sc = make_chain(n)
    name = (sc.Sz[0], sc.Sz[0])
    t0 = time.time()
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name,
                                        es=ES, delta=DELTA)
    dt = time.time() - t0
    y = np.array(y).real
    peak = np.array(x)[np.argmax(y)]
    return dt, peak, y


def run_truncated(n):
    sc = make_chain(n)
    sc.kpm_scale = NARROW_KPM_SCALE
    sc.kpm_energy_truncate = True
    sc.kpm_truncate_dK = TRUNC_DK
    sc.kpm_truncate_nsweeps = TRUNC_NSWEEPS
    name = (sc.Sz[0], sc.Sz[0])
    t0 = time.time()
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name,
                                        es=ES, delta=DELTA)
    dt = time.time() - t0
    y = np.array(y).real
    peak = np.array(x)[np.argmax(y)]
    return dt, peak, y


print(f"{'n':>3} {'untrunc [s]':>12} {'trunc [s]':>10} {'slowdown':>9} "
      f"{'peak untrunc':>13} {'peak trunc':>11}")
for n in (4, 6):
    dt_u, peak_u, y_u = run_untruncated(n)
    dt_t, peak_t, y_t = run_truncated(n)
    slowdown = dt_t / dt_u
    print(f"{n:3d} {dt_u:12.3f} {dt_t:10.3f} {slowdown:8.1f}x "
          f"{peak_u:13.3f} {peak_t:11.3f}")

    assert np.all(np.isfinite(y_u)) and np.all(np.isfinite(y_t))
    # both must locate the same resonance -- energy truncation (and any
    # optimization of its own internals) must not move the physics,
    # only the wall-clock time
    assert abs(peak_u - peak_t) <= 2 * (ES[1] - ES[0]), \
        f"n={n}: truncated/untruncated peaks disagree ({peak_t} vs {peak_u})"

print("OK: truncated and untruncated v3 KPM agree on the resonance location.")
