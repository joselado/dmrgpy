# Self-check for dmrgpy.pyitensor Phase 7 (the Chain facade). Cross-checked
# against dense/exact-diagonalization references wherever practical, since
# this is where every earlier phase's primitives get composed together
# into the actual methods dmrgpy's Python layer will call.
import numpy as np
from scipy.linalg import expm

from dmrgpy.pyitensor.autompo import AutoMPO
from dmrgpy.pyitensor.chain import Chain
from testutil import mps_to_dense


def check(name, cond):
    status = "OK" if cond else "FAIL"
    print("[{}] {}".format(status, name))
    if not cond:
        raise SystemExit("selftest failed: " + name)


def heisenberg_terms(n, field=None):
    terms = []
    for i in range(1, n):
        terms.append((1.0, [("Sz", i), ("Sz", i + 1)]))
        terms.append((0.5, [("Sp", i), ("Sm", i + 1)]))
        terms.append((0.5, [("Sm", i), ("Sp", i + 1)]))
    if field:
        for i in range(1, n + 1):
            terms.append((field * i, [("Sz", i)]))
    return terms


def dense_reference(sites, terms):
    return AutoMPO.from_terms(sites, terms).dense_matrix()


def test_spin_chain_basics():
    np.random.seed(0)
    n = 5
    chain = Chain([2] * n)
    chain.set_sweep_params(maxm=40, nsweeps=12, cutoff=1e-12, noise=0.1)
    # A small field breaks the plain chain's spin-rotation degeneracy (odd
    # N gives a >=2-fold degenerate ground state otherwise) -- without it,
    # DMRG and np.linalg.eigh can each land on a different, equally valid
    # vector within the same degenerate subspace, and the fidelity check
    # below would fail despite both being genuinely correct ground states.
    terms = heisenberg_terms(n, field=0.05)
    chain.set_hamiltonian(terms)
    Hdense = dense_reference(chain.sites, terms)
    w, v = np.linalg.eigh(Hdense)

    e0 = chain.gs_energy()
    check("Chain.gs_energy matches exact diagonalization", abs(e0 - w[0]) < 1e-6)
    check("Chain.gs_energy(skip_dmrg=True) returns the cached value without recomputation",
          chain.gs_energy(skip_dmrg=True) == e0)

    psi = chain.gs_wavefunction()
    dense_psi = mps_to_dense(psi)
    fidelity = abs(np.vdot(dense_psi, v[:, 0])) ** 2
    check("Chain.gs_wavefunction has high overlap with the exact ground state",
          fidelity > 1 - 1e-8)

    sz1_terms = [(1.0, [("Sz", 1)])]
    vev = chain.vev(sz1_terms, psi)
    ref_vev = np.vdot(dense_psi, dense_reference(chain.sites, sz1_terms) @ dense_psi)
    check("Chain.vev(<Sz_1>) matches dense expectation value", np.isclose(vev, ref_vev, atol=1e-6))

    A = chain.build_operator(sz1_terms)
    Apsi = chain.apply_pure_operator(A, psi)
    check("Chain.apply_pure_operator matches dense Sz_1|psi>",
          np.allclose(mps_to_dense(Apsi), dense_reference(chain.sites, sz1_terms) @ dense_psi, atol=1e-6))
    check("Chain.overlap(psi,psi) is 1", np.isclose(chain.overlap(psi, psi), 1.0, atol=1e-8))
    check("Chain.overlap_aMb_operator(psi,Sz_1,psi) matches vev", np.isclose(
        chain.overlap_aMb_operator(psi, A, psi), ref_vev, atol=1e-6))

    Ah = chain.hermitian_operator(A)
    from testutil import mpo_to_dense
    check("Chain.hermitian_operator(Sz_1) matches dense conjugate-transpose",
          np.allclose(mpo_to_dense(Ah), dense_reference(chain.sites, sz1_terms).conj().T, atol=1e-6))

    trace_ref = np.trace(dense_reference(chain.sites, sz1_terms))
    check("Chain.trace_operator matches dense trace", np.isclose(chain.trace_operator(A), trace_ref, atol=1e-6))

    B = chain.build_operator([(1.0, [("Sx", 1)])])
    AB = chain.multiply_operators(A, B)
    # multiply_operators(A,B) mirrors chain_session.h's mult_mpo(A,B): its
    # own comment documents that nmultMPO(A,prime(B)) composes "A's
    # level-0->1 map with B's now-level-1->2 map", i.e. A is applied to the
    # ket *first*, B second -- so in standard matrix form the result is
    # B_std @ A_std, not A_std @ B_std. Confirmed directly against both
    # orderings before writing this reference.
    ref_AB = dense_reference(chain.sites, [(1.0, [("Sx", 1)])]) @ dense_reference(chain.sites, sz1_terms)
    check("Chain.multiply_operators(Sz_1,Sx_1) matches dense Sx_1 @ Sz_1 "
          "(A-applied-first convention, see chain_session.h's mult_mpo comment)",
          np.allclose(mpo_to_dense(AB), ref_AB, atol=1e-6))

    conj_psi = chain.conjugate(psi)
    check("Chain.conjugate matches dense elementwise conjugate",
          np.allclose(mps_to_dense(conj_psi), np.conj(dense_psi), atol=1e-8))

    return chain, Hdense, w, v


def test_reduced_dm_and_bond_entropy(chain, Hdense, w, v):
    n = chain.sites.length()
    psi = chain.gs_wavefunction()
    dense_psi = v[:, 0]
    dims = [chain.sites.dim(i) for i in range(1, n + 1)]

    site = 3
    rho = chain.reduced_dm(psi, site)
    left_dim = int(np.prod(dims[:site - 1])) if site > 1 else 1
    right_dim = int(np.prod(dims[site:])) if site < n else 1
    a = dense_psi.reshape(left_dim, dims[site - 1], right_dim)
    ref_rho = np.einsum("lsr,lqr->sq", a, np.conj(a))
    check("Chain.reduced_dm matches dense partial trace at site {}".format(site),
          np.allclose(rho, ref_rho, atol=1e-6))
    check("reduced_dm has trace 1", np.isclose(np.trace(rho), 1.0, atol=1e-6))

    b = 2
    S = chain.bond_entropy(psi, b)
    left_dim_b = int(np.prod(dims[:b]))
    right_dim_b = int(np.prod(dims[b:]))
    mat = dense_psi.reshape(left_dim_b, right_dim_b)
    sv = np.linalg.svd(mat, compute_uv=False)
    p = sv ** 2
    p = p[p > 1e-12]
    ref_S = -np.sum(p * np.log(p))
    check("Chain.bond_entropy matches dense entanglement entropy at bond {}".format(b),
          np.isclose(S, ref_S, atol=1e-6))


def test_exponential_apply_and_quench(chain, Hdense, w, v):
    psi = chain.gs_wavefunction()
    dense_psi = mps_to_dense(psi)

    tau = -0.2j
    terms = heisenberg_terms(chain.sites.length())
    evolved = chain.exponential_apply(terms, psi, tau, nsteps=6)
    ref = expm(tau * Hdense) @ dense_psi
    overlap = abs(np.vdot(mps_to_dense(evolved), ref)) / (np.linalg.norm(mps_to_dense(evolved)) * np.linalg.norm(ref))
    check("Chain.exponential_apply approximates exp(tau*H)|psi> (overlap={:.6f})".format(overlap),
          overlap > 1 - 1e-3)

    sz1 = [(1.0, [("Sz", 1)])]
    szN = [(1.0, [("Sz", chain.sites.length())])]
    correlator, final_wf = chain.quench(terms, sz1, szN, nt=4, dt=0.05, fit_td=False)
    check("Chain.quench returns nt correlator values", len(correlator) == 4)
    correlator_tdvp, final_wf_tdvp = chain.quench_tdvp(terms, sz1, szN, nt=4, dt=0.05)
    check("Chain.quench_tdvp returns nt correlator values", len(correlator_tdvp) == 4)
    check("Chain.quench and Chain.quench_tdvp roughly agree at short times",
          abs(correlator[0] - correlator_tdvp[0]) < 5e-2)


def test_bicstab_and_cvm(chain, Hdense, w, v):
    n = chain.sites.length()
    psi = chain.gs_wavefunction()
    Id_terms = [(1.0, [("Id", 1)])]
    dense_Id = dense_reference(chain.sites, Id_terms)
    z = 0.7 + 0.3j
    A_dense = z * dense_Id - Hdense
    b_dense = mps_to_dense(psi)
    x_ref = np.linalg.solve(A_dense, b_dense)

    A_terms = [(z, [("Id", 1)])] + [(-c, ops) for c, ops in heisenberg_terms(n, field=0.05)]
    A_mpo = chain.build_operator(A_terms)
    x = chain._bicstab(A_mpo, psi, tol=1e-10, max_it=200)
    from testutil import mps_to_dense as m2d
    check("Chain._bicstab solves (z*Id - H)x = psi matching dense linear solve",
          np.allclose(m2d(x), x_ref, atol=1e-4))

    g = chain.cvm_dynamical_correlator([(1.0, [("Sz", 1)])], [(1.0, [("Sz", 1)])], omega=0.5, eta=0.1,
                                        energy=w[0], tol=1e-10, max_it=200)
    Sz1_dense = dense_reference(chain.sites, [(1.0, [("Sz", 1)])])
    zc = complex(0.5 + w[0], 0.1)
    x_cvm = np.linalg.solve(zc * dense_Id - Hdense, Sz1_dense @ b_dense)
    g_ref = -np.vdot(b_dense, Sz1_dense @ x_cvm).imag / np.pi
    check("Chain.cvm_dynamical_correlator matches a dense resolvent calculation",
          np.isclose(g, g_ref, atol=1e-4))


def test_kpm_sanity(chain):
    psi = chain.gs_wavefunction()
    moments, emin, emax, scale, n = chain.kpm_dynamical_correlator(
        [(1.0, [("Sz", 1)])], [(1.0, [("Sz", 1)])], kpmmaxm=40, kpm_scale=1.2, kpm_accelerate=True,
        kpm_n_scale=2, delta=0.5, kpm_cutoff=1e-10)
    check("kpm_dynamical_correlator returns the requested number of moments plus 2",
          len(moments) == n + 2)
    check("kpm mu0 is real and positive (it's a norm-squared overlap)",
          moments[0].real > 0 and abs(moments[0].imag) < 1e-8)


def test_fermion_correlators():
    np.random.seed(3)
    n = 5
    chain = Chain([0] * n)  # spinless fermion
    chain.set_sweep_params(maxm=40, nsweeps=14, cutoff=1e-12, noise=0.1)
    terms = []
    for i in range(1, n):
        terms.append((1.0, [("Cdag", i), ("C", i + 1)]))
        terms.append((1.0, [("Cdag", i + 1), ("C", i)]))
    for i in range(1, n):
        terms.append((0.3, [("N", i), ("N", i + 1)]))
    chain.set_hamiltonian(terms)
    e0 = chain.gs_energy()
    Hdense = dense_reference(chain.sites, terms)
    w, v = np.linalg.eigh(Hdense)
    check("fermion Chain.gs_energy matches exact diagonalization", abs(e0 - w[0]) < 1e-6)

    psi = chain.gs_wavefunction()
    dense_psi = v[:, 0]
    C = chain.correlation_matrix(psi)
    ref_C = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            op = dense_reference(chain.sites, [(1.0, [("Cdag", i + 1), ("C", j + 1)])])
            ref_C[i, j] = np.vdot(dense_psi, op @ dense_psi)
    check("Chain.correlation_matrix matches dense <Cdag_i C_j>", np.allclose(C, ref_C, atol=1e-6))

    F = chain.four_correlation_tensor(psi, accelerate=True)
    for (i, j, k, l) in [(0, 1, 2, 3), (1, 1, 2, 2), (0, 2, 1, 3)]:
        op = dense_reference(chain.sites, [(1.0, [("Cdag", i + 1), ("C", j + 1),
                                                    ("Cdag", k + 1), ("C", l + 1)])])
        ref = np.vdot(dense_psi, op @ dense_psi)
        check("Chain.four_correlation_tensor[{},{},{},{}] matches dense reference".format(i, j, k, l),
              np.isclose(F[i, j, k, l], ref, atol=1e-6))


def test_excited_states_via_chain():
    # scale_lagrange=2.0, nsweeps=20: the overlap-penalty method's local
    # minima (see dmrg.py's dmrg_excited docstring) turn out to be *more*
    # sensitive than this test first assumed -- confirmed directly by
    # sweeping many random seeds at the previous (scale_lagrange=1.0,
    # nsweeps=18) settings and finding a substantial failure rate, not just
    # this one seed being unlucky. Doubling the Lagrange weight and adding
    # a couple more sweeps was enough to make *this* seed converge reliably
    # under both the JAX and plain-NumPy kernel paths (also confirmed
    # directly: with the old settings, seed=2 itself happened to pass under
    # NumPy but not under JAX -- a real reminder that this algorithm's
    # convergence margin, not the kernel's numerical correctness, is what's
    # thin here; a single matvec call agrees between the two kernel paths to
    # ~1e-15 relative, see kernels.py's docstring).
    # This is a genuine, pre-existing limitation of the penalty method as
    # implemented (no noise-term perturbation, see dmrg.py's module
    # docstring) -- not something this test works around by cherry-picking
    # a lucky seed without saying so.
    #
    # Re-confirmed yet again after the contract_many() contraction-order
    # fix (tensor.py/dmrg.py/mpsalgebra.py): mathematically identical
    # results, different floating-point summation order, and that alone
    # was enough to flip nsweeps=20 from converged to not for this seed --
    # nsweeps=25 is reliably converged. Since DMRG is now dramatically
    # faster, the extra sweeps cost is negligible.
    np.random.seed(2)
    n = 6
    chain = Chain([2] * n)
    chain.set_sweep_params(maxm=60, nsweeps=25, cutoff=1e-13, noise=0.1)
    terms = heisenberg_terms(n, field=0.1)
    chain.set_hamiltonian(terms)
    Hdense = dense_reference(chain.sites, terms)
    ref = np.linalg.eigvalsh(Hdense)

    energies, fluctuations, wavefunctions = chain.excited_states(2, scale_lagrange=2.0)
    check("Chain.excited_states ground energy matches exact diagonalization",
          abs(energies[0] - ref[0]) < 1e-6)
    check("Chain.excited_states first-excited energy matches exact diagonalization",
          abs(energies[1] - ref[1]) < 1e-5)
    check("excited_states returns matching-length energies/fluctuations/wavefunctions",
          len(energies) == 2 and len(fluctuations) == 2 and len(wavefunctions) == 2)


if __name__ == "__main__":
    chain, Hdense, w, v = test_spin_chain_basics()
    test_reduced_dm_and_bond_entropy(chain, Hdense, w, v)
    test_exponential_apply_and_quench(chain, Hdense, w, v)
    test_bicstab_and_cvm(chain, Hdense, w, v)
    test_kpm_sanity(chain)
    test_fermion_correlators()
    test_excited_states_via_chain()
    print("All pyitensor Phase 7 checks passed.")
