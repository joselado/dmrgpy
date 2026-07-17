# Self-check for dmrgpy.pyitensor Phase 6 (two-site TDVP real-time
# evolution). Cross-checked against exact dense time evolution
# (scipy.linalg.expm(-1j*dt*H) applied directly to the dense state vector)
# for small systems where maxdim is generous enough that SVD truncation
# isn't the limiting factor -- so this isolates TDVP's own algorithm
# (forward + backward Krylov steps) from any truncation error.
import numpy as np
from scipy.linalg import expm

from dmrgpy.pyitensor import AutoMPO, to_mpo
from dmrgpy.pyitensor.sites import SiteX
from dmrgpy.pyitensor.tdvp import tdvp_step
from testutil import dense_to_mps, mps_to_dense


def check(name, cond):
    status = "OK" if cond else "FAIL"
    print("[{}] {}".format(status, name))
    if not cond:
        raise SystemExit("selftest failed: " + name)


def test_heisenberg_time_evolution():
    np.random.seed(0)
    n = 5
    sites = SiteX([2] * n)
    ampo = AutoMPO(sites)
    for i in range(1, n):
        ampo.add(1.0, "Sz", i, "Sz", i + 1)
        ampo.add(0.5, "Sp", i, "Sm", i + 1)
        ampo.add(0.5, "Sm", i, "Sp", i + 1)
    Hdense = ampo.dense_matrix()
    H = to_mpo(ampo, cutoff=0.0, maxdim=2000)

    v0 = np.random.randn(2 ** n) + 1j * np.random.randn(2 ** n)
    v0 /= np.linalg.norm(v0)
    psi = dense_to_mps(v0, sites, maxdim=2000)
    check("dense_to_mps reproduces the starting state",
          np.allclose(mps_to_dense(psi), v0, atol=1e-10))

    # Two-site TDVP's forward/backward split is a 2nd-order-in-dt
    # integrator (confirmed directly: halving dt quarters the infidelity
    # 1-|<exact|tdvp>| -- 2.5e-5 -> 6.3e-6 -> 1.6e-6 for dt=0.05/0.025/
    # 0.0125), so a real, expected O(dt^2) integration error is baked into
    # this comparison -- this test only checks that error is small and
    # shrinks correctly with dt, not that it vanishes at a fixed dt.
    def evolve_and_compare(dt, nsteps):
        psi_t = dense_to_mps(v0, sites, maxdim=2000)
        v_exact = v0.copy()
        Udt = expm(-1j * dt * Hdense)
        for _ in range(nsteps):
            tdvp_step(psi_t, H, dt, cutoff=1e-13, maxdim=64, niter=40)
            v_exact = Udt @ v_exact
        return mps_to_dense(psi_t), v_exact

    v_tdvp, v_exact = evolve_and_compare(dt=0.05, nsteps=6)
    check("TDVP preserves the norm", np.isclose(np.vdot(v_tdvp, v_tdvp), 1.0, atol=1e-6))
    infidelity_coarse = 1 - abs(np.vdot(v_exact, v_tdvp))
    check("TDVP matches exact time evolution to within 2nd-order Trotter "
          "error at dt=0.05 (infidelity={:.2e})".format(infidelity_coarse),
          infidelity_coarse < 1e-4)

    v_tdvp_fine, v_exact_fine = evolve_and_compare(dt=0.0125, nsteps=24)
    infidelity_fine = 1 - abs(np.vdot(v_exact_fine, v_tdvp_fine))
    check("halving dt twice (0.05 -> 0.0125) shrinks the infidelity by "
          "~16x as expected for a 2nd-order integrator ({:.2e} -> {:.2e})".format(
              infidelity_coarse, infidelity_fine),
          infidelity_fine < infidelity_coarse / 8)


def test_fermion_hopping_time_evolution():
    np.random.seed(1)
    n = 4
    sites = SiteX([0] * n)
    ampo = AutoMPO(sites)
    for i in range(1, n):
        ampo.add(1.0, "Cdag", i, "C", i + 1)
        ampo.add(1.0, "Cdag", i + 1, "C", i)
    ampo.add(0.3, "Cdag", 1, "C", n)
    ampo.add(0.3, "Cdag", n, "C", 1)
    Hdense = ampo.dense_matrix()
    H = to_mpo(ampo, cutoff=0.0, maxdim=2000)

    v0 = np.random.randn(2 ** n) + 1j * np.random.randn(2 ** n)
    v0 /= np.linalg.norm(v0)
    psi = dense_to_mps(v0, sites, maxdim=2000)

    dt = 0.05
    nsteps = 6
    v_exact = v0.copy()
    Udt = expm(-1j * dt * Hdense)
    for _ in range(nsteps):
        tdvp_step(psi, H, dt, cutoff=1e-13, maxdim=64, niter=40)
        v_exact = Udt @ v_exact

    v_tdvp = mps_to_dense(psi)
    infidelity = 1 - abs(np.vdot(v_exact, v_tdvp))
    # Same 2nd-order-in-dt Trotter error as the spin case (see its comment);
    # tolerance set the same way, not tightened to exact match.
    check("TDVP matches exact time evolution for a fermionic (JW-threaded) "
          "Hamiltonian (infidelity={:.2e})".format(infidelity), infidelity < 1e-4)


if __name__ == "__main__":
    test_heisenberg_time_evolution()
    test_fermion_hopping_time_evolution()
    print("All pyitensor Phase 6 checks passed.")
