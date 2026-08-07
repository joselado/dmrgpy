# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test for a real physics bug: thermal.py's anneal() used to
# dose the purified wavefunction with the FULL beta=1/T instead of
# beta/2, i.e. it applied |Psi(beta)> = e^{-beta H}|Psi(0)> instead of
# the standard purification convention |Psi(beta)> = e^{-beta H/2}|Psi(0)>
# (needed so that tracing out the ancilla gives rho ~ e^{-beta H}, the
# actual thermal state of T). Confirmed directly: Thermal_Spin_Chain(T=1)
# used to report the thermal energy of the exact T=0.5 state, not T=1 --
# i.e. every finite-T calculation through this class was silently
# computing at T/2. Fixed by using beta/2 as the imaginary-time dose.
# This test checks the purification-based thermal energy directly
# against an exact Boltzmann-weighted average over ED eigenvalues.
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
from dmrgpy import thermal

n = 4 # small enough for exact diagonalization
spins = ["S=1/2" for i in range(n)]

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

# purification-based thermal energy (mode="ED" so no compiled extension
# is needed for the ancilla chain's ground state / annealing steps)
h = 0
Ts = [0.25, 0.5, 1.0, 2.0]
tol = 0.05
e_exacts, e_purifications = [], []
for T in Ts:
    e_exact = exact_thermal_energy(T)
    print("Exact thermal energy (T=%.2f) ="%T, e_exact)

    tc = thermal.Thermal_Spin_Chain(spins)
    h = 0
    for i in range(n-1):
        h = h + tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
    tc.set_hamiltonian(h)
    tc.T = T
    tc.mode = "ED"
    wf = tc.get_gs()
    e_purification = wf.dot(h*wf).real
    print("Purification thermal energy (T=%.2f) ="%T, e_purification)

    diff = abs(e_purification-e_exact)
    print("Difference =",diff)

    e_exacts.append(e_exact) ; e_purifications.append(e_purification)

    # Generous tolerance: anneal() uses first-order Euler steps (dbeta=0.1
    # by default), so some discretization error versus the exact Boltzmann
    # average is expected -- this tolerance is meant to catch a real
    # regression (e.g. the old missing-factor-of-2 bug, which was off by
    # nearly 2x in energy, not this), not to pin down the discretization.
    assert diff<tol, "purification thermal energy disagrees with exact ED by %g (tol=%g) at T=%g -- check thermal.py's anneal() beta schedule"%(diff,tol,T)

print("TEST PASSED")

plt.plot(Ts, e_exacts, "k-o", label="Exact ED")
plt.plot(Ts, e_purifications, "s--", label="Purification DMRG (ED mode)")
plt.xlabel("Temperature T")
plt.ylabel(r"Thermal energy $\langle H\rangle_T$")
plt.title("Purification thermal energy vs exact ED, n=%d Heisenberg chain"%n)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
