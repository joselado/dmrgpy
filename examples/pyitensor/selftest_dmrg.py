# Self-check for dmrgpy.pyitensor Phase 5 (two-site ground-state DMRG).
# Cross-checked against exact diagonalization: AutoMPO.dense_matrix() for
# the Heisenberg case, dmrgpy's independent ED backend (pyfermion) for the
# fermionic case.
import numpy as np

from dmrgpy.pyitensor import AutoMPO, Sweeps, inner, randomMPS, to_mpo
from dmrgpy.pyitensor.dmrg import dmrg, dmrg_excited
from dmrgpy.pyitensor.sites import SiteX
from dmrgpy.pyfermion.mbfermion import MBFermion


def check(name, cond):
    status = "OK" if cond else "FAIL"
    print("[{}] {}".format(status, name))
    if not cond:
        raise SystemExit("selftest failed: " + name)


def run_dmrg(sites, ampo, maxdim=40, nsweep=12):
    H = to_mpo(ampo, cutoff=1e-14, maxdim=5000)
    psi = randomMPS(sites, maxdim)
    sweeps = Sweeps(nsweep)
    sweeps.maxdim = maxdim
    sweeps.cutoff = 1e-12
    sweeps.niter = 30
    energy = dmrg(psi, H, sweeps, quiet=True)
    return energy, psi


def test_heisenberg_chain():
    np.random.seed(0)
    n = 8
    sites = SiteX([2] * n)
    ampo = AutoMPO(sites)
    for i in range(1, n):
        ampo.add(1.0, "Sz", i, "Sz", i + 1)
        ampo.add(0.5, "Sp", i, "Sm", i + 1)
        ampo.add(0.5, "Sm", i, "Sp", i + 1)
    ref = np.linalg.eigvalsh(ampo.dense_matrix())[0].real

    energy, psi = run_dmrg(sites, ampo)
    check("Heisenberg-8 DMRG ground energy matches exact diagonalization "
          "(ref={:.8f}, dmrg={:.8f})".format(ref, energy),
          abs(energy - ref) < 1e-6)


def test_heisenberg_chain_open_random_restart():
    # Same Hamiltonian, different random seed/starting MPS -- DMRG should
    # converge to the same ground energy regardless of the random start.
    np.random.seed(42)
    n = 8
    sites = SiteX([2] * n)
    ampo = AutoMPO(sites)
    for i in range(1, n):
        ampo.add(1.0, "Sz", i, "Sz", i + 1)
        ampo.add(0.5, "Sp", i, "Sm", i + 1)
        ampo.add(0.5, "Sm", i, "Sp", i + 1)
    ref = np.linalg.eigvalsh(ampo.dense_matrix())[0].real
    energy, psi = run_dmrg(sites, ampo)
    check("Heisenberg-8 DMRG is seed-independent (ref={:.8f}, dmrg={:.8f})".format(ref, energy),
          abs(energy - ref) < 1e-6)


def test_fermion_hopping_chain():
    np.random.seed(1)
    n = 6
    sites = SiteX([0] * n)  # spinless fermion
    ampo = AutoMPO(sites)
    for i in range(1, n):
        ampo.add(1.0, "Cdag", i, "C", i + 1)
        ampo.add(1.0, "Cdag", i + 1, "C", i)
    # a small interaction term to make this more than free fermions
    for i in range(1, n):
        ampo.add(0.5, "N", i, "N", i + 1)

    mb = MBFermion(n)
    ref_h = mb.get_zero()
    for i in range(n - 1):
        ref_h = ref_h + mb.Cdag[i] @ mb.C[i + 1] + mb.Cdag[i + 1] @ mb.C[i]
        ref_h = ref_h + 0.5 * (mb.N[i] @ mb.N[i + 1])
    ref = np.linalg.eigvalsh(ref_h.toarray())[0].real

    energy, psi = run_dmrg(sites, ampo, maxdim=40, nsweep=14)
    check("Interacting fermion chain DMRG ground energy matches independent "
          "ED backend (ref={:.8f}, dmrg={:.8f})".format(ref, energy),
          abs(energy - ref) < 1e-6)


def test_excited_states_penalty_method():
    # A non-uniform field breaks the Heisenberg chain's degeneracies so the
    # first excited state is unambiguous, then dmrg_excited() (the
    # overlap-penalty method chain_session.h's excited_states() uses) is
    # cross-checked against the two lowest eigenvalues of exact
    # diagonalization, and the found states' mutual overlap is checked to
    # be small (the penalty actually did its job).
    np.random.seed(2)
    n = 6
    sites = SiteX([2] * n)
    ampo = AutoMPO(sites)
    for i in range(1, n):
        ampo.add(1.0, "Sz", i, "Sz", i + 1)
        ampo.add(0.5, "Sp", i, "Sm", i + 1)
        ampo.add(0.5, "Sm", i, "Sp", i + 1)
    for i in range(1, n + 1):
        ampo.add(0.1 * i, "Sz", i)
    ref = np.linalg.eigvalsh(ampo.dense_matrix())
    H = to_mpo(ampo, cutoff=1e-14, maxdim=5000)

    # The penalty method's own local minima mean it can wander for several
    # sweeps before settling (confirmed directly: 12 sweeps sometimes still
    # oscillating, 18+ reliably converged) -- generous sweep/niter counts
    # here are about giving the *algorithm* room to converge, not a
    # precision knob. Re-confirmed again after the contract_many()
    # contraction-order fix (same energies, different floating-point
    # summation path -- see tensor.py's contract_many() docstring): 18
    # sweeps at weight=1x bandwidth was right at the edge and flipped to
    # not-converged under the new (still mathematically equivalent) path;
    # 25 sweeps at weight=2x bandwidth is reliably converged, matching the
    # same margin chain.py's own excited_states() test needed. Since DMRG
    # is now dramatically faster (see dmrg.py/mpsalgebra.py's
    # contract_many() migration), a few extra sweeps costs nothing.
    sweeps = Sweeps(25)
    sweeps.maxdim = 60
    sweeps.cutoff = 1e-13
    sweeps.niter = 50

    psi0 = randomMPS(sites, 60)
    e0 = dmrg(psi0, H, sweeps, quiet=True)
    check("excited-states test: ground energy matches exact diagonalization",
          abs(e0 - ref[0].real) < 1e-6)

    bandwidth = ref[-1].real - ref[0].real
    weight = 2.0 * bandwidth  # mirrors chain_session.h's weight = bandwidth()*scale_lagrange
    psi1 = randomMPS(sites, 60)
    dmrg_excited(psi1, H, [psi0], weight, sweeps, quiet=True)
    psi1.normalize()
    e1 = inner(psi1, H, psi1).real

    check("first excited-state energy matches exact diagonalization "
          "(ref={:.8f}, dmrg={:.8f})".format(ref[1].real, e1), abs(e1 - ref[1].real) < 1e-5)
    check("the penalty method actually suppressed overlap with the ground state",
          abs(inner(psi0, psi1)) < 1e-4)


if __name__ == "__main__":
    test_heisenberg_chain()
    test_heisenberg_chain_open_random_restart()
    test_fermion_hopping_chain()
    test_excited_states_penalty_method()
    print("All pyitensor Phase 5 checks passed.")
