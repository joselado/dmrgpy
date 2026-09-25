# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test for two real physics bugs in thermal.py's anneal():
#
# 1. It used to dose the purified wavefunction with the FULL beta=1/T
#    instead of beta/2, i.e. it applied |Psi(beta)> = e^{-beta H}|Psi(0)>
#    instead of the standard purification convention
#    |Psi(beta)> = e^{-beta H/2}|Psi(0)> (needed so that tracing out the
#    ancilla gives rho ~ e^{-beta H}, the actual thermal state of T).
#    Thermal_Spin_Chain(T=1) reported the thermal energy of the exact T=0.5
#    state, off by nearly 2x.
# 2. Its steps were first-order, (1 - 0.1*H) on the unshifted H, so the
#    error was set by 0.1 times the absolute, extensive energy: 4.2 per
#    cent at n=4 (where this example used to run, at an absolute tolerance
#    of 0.05 that it passed by 0.003), 7.0 per cent at n=6 and 12.3 per
#    cent at n=10, and it moved with a constant energy offset (2026-09-25b
#    audit, finding 18). The steps are now an order-8 Taylor polynomial of
#    exp(-dtau*(H-E_ref)) with dtau*W fixed, W the spectral width, so the
#    error no longer depends on the chain length, the units or an offset.
#
# This checks the purification-based thermal energy directly against an
# exact Boltzmann-weighted average over ED eigenvalues, at a RELATIVE
# tolerance, on a chain long enough for the old error to have shown, and
# with the Hamiltonian shifted by a constant, which must not change
# <H+c>-c.
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
from dmrgpy import thermal

n = 6 # small enough for exact diagonalization of the doubled chain
spins = ["S=1/2" for i in range(n)]
offset = 20. # a constant energy offset, which the Gibbs state ignores

# exact reference: Tr[H exp(-H/T)] / Tr[exp(-H/T)] via full ED spectrum
sc_ed = spinchain.Spin_Chain(spins)
h_ed = 0
for i in range(n-1):
    h_ed = h_ed + sc_ed.Sx[i]*sc_ed.Sx[i+1] + sc_ed.Sy[i]*sc_ed.Sy[i+1] + sc_ed.Sz[i]*sc_ed.Sz[i+1]
sc_ed.set_hamiltonian(h_ed)
Hmat = sc_ed.get_ED_obj().get_hamiltonian()
Hmat = Hmat.toarray() if hasattr(Hmat,"toarray") else Hmat
evals = np.linalg.eigvalsh(Hmat)

def exact_thermal_energy(T):
    w = np.exp(-(evals-evals.min())/T)
    return np.sum(evals*w)/np.sum(w)

def purification_energy(T,c=0.):
    """<H+c>-c in the purified state of H+c at temperature T (mode="ED",
    so no compiled extension is needed for the ancilla chain's ground
    state or the annealing steps)"""
    tc = thermal.Thermal_Spin_Chain(spins)
    h = 0
    for i in range(n-1):
        h = h + tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
    h = h + c
    tc.set_hamiltonian(h)
    tc.T = T
    tc.mode = "ED"
    wf = tc.get_gs()
    return wf.dot(h*wf).real - c

Ts = [0.1, 0.25, 0.5, 1.0, 2.0, 10.0]
rtol = 1e-5 # relative; the stepper's own error is below 1e-6 at its defaults
e_exacts, e_purifications, e_shifted = [], [], []
for T in Ts:
    e_exact = exact_thermal_energy(T)
    e_purification = purification_energy(T)
    e_shift = purification_energy(T,c=offset)
    print("T=%.2f  exact %.8f  purification %.8f  with offset %.0f: %.8f"
          %(T,e_exact,e_purification,offset,e_shift))
    e_exacts.append(e_exact) ; e_purifications.append(e_purification)
    e_shifted.append(e_shift)
    for e in (e_purification,e_shift):
        diff = abs(e-e_exact)/abs(e_exact)
        assert diff<rtol, "purification thermal energy disagrees with exact ED by %g (relative, rtol=%g) at T=%g -- check thermal.py's anneal()"%(diff,rtol,T)

print("TEST PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(Ts, e_exacts, "k-o", label="Exact ED")
ax1.plot(Ts, e_purifications, "s--", label="Purification (ED mode)")
ax1.plot(Ts, e_shifted, "x:", label=r"Purification of $H+%.0f$, minus %.0f"%(offset,offset))
ax1.set_xscale("log")
ax1.set_xlabel("Temperature T")
ax1.set_ylabel(r"Thermal energy $\langle H\rangle_T$")
ax1.set_title("n=%d Heisenberg chain"%n)
ax1.legend()
ax1.grid(alpha=0.3)
ax2.semilogy(Ts, np.abs(np.array(e_purifications)-e_exacts)/np.abs(e_exacts), "s--", label=r"$H$")
ax2.semilogy(Ts, np.abs(np.array(e_shifted)-e_exacts)/np.abs(e_exacts), "x:", label=r"$H+%.0f$"%offset)
ax2.axhline(rtol, color="gray", lw=0.8, label="tolerance")
ax2.set_xscale("log")
ax2.set_xlabel("Temperature T")
ax2.set_ylabel("relative error")
ax2.legend()
ax2.grid(alpha=0.3)
plt.tight_layout()
plt.show()
