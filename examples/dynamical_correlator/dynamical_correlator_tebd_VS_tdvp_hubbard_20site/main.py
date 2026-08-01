# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Benchmark: dynamical correlator (submode="TD", the real-time-evolution
# + Fourier-transform method, see timedependent.py) computed with
# tevol_method="TEBD" vs the default tevol_method="TDVP", on a 20-orbital
# native Hubbard chain (fermionchain.Spinful_Fermionic_Chain_Native --
# one dimension-4 (Empty/Up/Down/UpDn) site per orbital, itensor_version=
# "python" only). Unlike the earlier TEBD-vs-TDVP benchmark for this same
# model (single-operator evolve_and_measure_dmrg -- see CLAUDE.md/PR #107
# follow-up), this exercises the actual two-operator quench dynamical
# correlator (Chain.quench_tebd()/quench_tdvp()): apply A to the ground
# state, evolve, and overlap against B at every step, giving C(t) =
# <GS|A(t)B(0)|GS> and, after windowing + FFT, the frequency-domain
# S_AB(omega).
#
# The standard interleaved Spinful_Fermionic_Chain isn't nearest-neighbor
# after Jordan-Wigner threading (hopping skips over the other spin
# flavor's site), so TEBD would reject it outright -- the native class is
# what makes both hopping and the onsite U term genuinely 2-site-local
# (see tebd.py's bond_hamiltonians()/_true_span()).
#
# TDVP and TEBD are run for a *different* number of steps here, not the
# same nt like the smaller-system example
# (dynamical_correlator_time_evolution_tebd/main.py): measured directly,
# TDVP's per-step cost at this system size/local dimension (4, i.e. 16
# for a 2-site block, at maxm=80) is ~200 s/step -- a full 30-step run
# would take on the order of two hours, impractical for a benchmark
# script. TEBD's own gates are built once up front from the bare bond
# Hamiltonians and reused unchanged every step, so its cost barely
# depends on nt in comparison, and is run for the full NT_TEBD steps.
# TDVP is only run for the much smaller NT_TDVP steps -- enough to (a)
# directly cross-check C(t) against TEBD's own trajectory over the steps
# they share (a single real-time trajectory is deterministic, so TEBD's
# first NT_TDVP entries from its NT_TEBD-step run are exactly what a
# standalone NT_TDVP-step TEBD run would have produced -- no separate run
# needed) and (b) get a real measured per-step cost to extrapolate the
# total NT_TEBD-step cost from, rather than actually waiting for it.
#
# IMPORTANT CAVEAT, found while building this benchmark: at this system
# size, TEBD and TDVP do NOT agree tightly (~0.1, on values of order
# ~0.5) the way they do at small n (see
# tests/test_time_evolution.py's ED cross-checks, which confirm both
# individually match *exact diagonalization* to ~3e-4 at n=4/n=6) or at
# moderate n=12 (~3e-3 over 6 steps at the same maxm=80). Diagnosed
# directly (not assumed) by inspecting bond dimensions after each TEBD
# step: at maxm=80, essentially the entire bulk of the chain is already
# at the maxm cap from the very first step, for *both* n=12 and n=20 --
# so this isn't "n=20 truncates harder" in a simple sense. What differs
# is the *number* of saturated bonds (~16 at n=20 vs ~7 at n=12): TDVP's
# per-bond Krylov-approximated local update and TEBD's per-bond exact
# bond-Hamiltonian-exponential update are two genuinely different
# approximations once a bond is truncation-saturated (neither reproduces
# the true local propagator exactly there), and those per-bond
# differences compound across more bonds at n=20. This is a real
# numerical-methods characteristic of running two different truncated
# real-time MPS propagators at a hard, saturated bond-dimension cap on a
# strongly-interacting model -- not a bug in either implementation (ruled
# out by the ED cross-checks above) -- so this script reports the
# discrepancy rather than asserting a tolerance that heavy truncation
# can't be expected to satisfy.
import time
import numpy as np
from dmrgpy import fermionchain
from dmrgpy import timedependent

n = 20
t = 1.0
U = 4.0

fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
h = 0
for i in range(n - 1):
    h = h + t*(fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
for i in range(n):
    h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
h = h + h.get_dagger()
fc.setup_python() # TEBD only exists on the pure-Python backend
fc.maxm = 80
fc.nsweeps = 16
fc.set_hamiltonian(h)

print(f"### Ground state (n={n} native Hubbard orbitals, t={t}, U={U}) ###")
t0 = time.time()
e0 = fc.gs_energy(mode="DMRG")
t_gs = time.time() - t0
print(f"E0 = {e0:.6f}   wall = {t_gs:.1f} s")
print()

name = (fc.Cup[0], fc.Cdagup[0]) # <Cup_0(t) Cdagup_0(0)>, single-particle Green's function
dt = 0.05
NT_TEBD = 30
NT_TDVP = 4

print(f"### Dynamical correlator: quench + measure, dt={dt} ###")

fc.tevol_method = "TEBD"
t0 = time.time()
(ts_tebd, cs_tebd) = timedependent.evolution_DC(fc, mode="DMRG", name=name, nt=NT_TEBD, dt=dt)
t_tebd = time.time() - t0
cs_tebd = np.array(cs_tebd)
print(f"TEBD: {NT_TEBD} steps, wall = {t_tebd:.2f} s ({t_tebd/NT_TEBD:.3f} s/step)")

fc.tevol_method = "TDVP"
t0 = time.time()
(ts_tdvp, cs_tdvp) = timedependent.evolution_DC(fc, mode="DMRG", name=name, nt=NT_TDVP, dt=dt)
t_tdvp = time.time() - t0
cs_tdvp = np.array(cs_tdvp)
t_tdvp_per_step = t_tdvp / NT_TDVP
print(f"TDVP: {NT_TDVP} steps, wall = {t_tdvp:.2f} s ({t_tdvp_per_step:.3f} s/step)")
print()

diff = np.max(np.abs(cs_tebd[:NT_TDVP] - cs_tdvp))
print(f"max |C_TEBD(t) - C_TDVP(t)| over the shared {NT_TDVP} steps = {diff:.3e}")
print("(large at this maxm/system size due to bond-dimension-saturation "
      "sensitivity -- see this file's module docstring; both propagators "
      "are separately validated against exact diagonalization in "
      "tests/test_time_evolution.py)")
print(f"TEBD speedup over TDVP, per step: {t_tdvp_per_step/(t_tebd/NT_TEBD):.1f}x")
print(f"Extrapolated TDVP wall for {NT_TEBD} steps: {t_tdvp_per_step*NT_TEBD:.0f} s "
      f"(~{t_tdvp_per_step*NT_TEBD/60:.0f} min), vs TEBD's measured {t_tebd:.1f} s")
