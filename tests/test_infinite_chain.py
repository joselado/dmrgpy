"""Coverage for Infinite_Many_Body_Chain / Infinite_Spin_Chain (iDMRG,
pyitensor/idmrg.py, or mpscpp3/chain_session.h's Chain::idmrg_ground_state
for itensor_version=3) -- see infinitechain.py's own docstring for the
L/C/R-suffixed Hamiltonian convention this exercises. Most tests below run
only the default itensor_version="python" backend; the
test_itensor_version3_* tests near the end of this file are the only ones
exercising the ITensor v3 C++ backend, skipped automatically if mpscpp3
isn't compiled.

There is no ED oracle for a genuinely infinite system, so these tests
cross-check via: the exact Bethe-ansatz Heisenberg energy density
(1/4 - ln(2)), a model-agnostic self-consistency identity (<H_uc> must
equal n_uc*e0_density by translational invariance, independent of
whether the *value* of e0_density itself is even known in closed form),
and a finite-chain-size extrapolation cross-check via ED (mode="ED",
so this needs no compiled C++ backend).

n_uc>=3 is out of scope for now (Infinite_Many_Body_Chain's constructor
raises NotImplementedError -- see its own docstring and idmrg.py's
module docstring for why), so only n_uc in {1, 2} is exercised here.
"""
import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy import infinitechain
from dmrgpy import spinchain
from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor

EXACT_HEISENBERG_DENSITY = 0.25 - np.log(2)


def _converged_uniform_chain(n_uc, maxm=40, maxiter=200, etol=1e-9, itensor_version="python"):
    spins = ["1/2"] * n_uc
    ic = infinitechain.Infinite_Spin_Chain(spins, itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    terms = []
    for i in range(n_uc - 1):
        terms.append(ic.SxC[i] * ic.SxC[i + 1] + ic.SyC[i] * ic.SyC[i + 1]
                     + ic.SzC[i] * ic.SzC[i + 1])
    terms.append(ic.SxC[n_uc - 1] * ic.SxR[0] + ic.SyC[n_uc - 1] * ic.SyR[0]
                 + ic.SzC[n_uc - 1] * ic.SzR[0])
    h = terms[0]
    for t in terms[1:]:
        h = h + t
    ic.maxm = maxm
    ic.maxiter = maxiter
    ic.etol = etol
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic, h


def test_n_uc1_heisenberg_energy_density_matches_bethe_ansatz():
    """The uniform antiferromagnetic Heisenberg chain is gapless
    (critical): iDMRG's energy-density finite-difference convergence
    criterion (`.converged`) is a power-law, not exponential, in this
    regime at any *fixed* bond dimension, so it's not expected to trip
    within a practical `maxiter` here -- only the *value* is checked,
    with a tolerance loose enough for maxm=40's residual truncation
    error (confirmed independently: maxm=60 gets within ~3e-6)."""
    ic, _ = _converged_uniform_chain(1)
    assert ic.e0 == pytest.approx(EXACT_HEISENBERG_DENSITY, abs=1e-4)


def test_n_uc2_uniform_heisenberg_matches_bethe_ansatz():
    ic, _ = _converged_uniform_chain(2)
    assert ic.e0 == pytest.approx(EXACT_HEISENBERG_DENSITY, abs=1e-4)


def test_unit_cell_expectation_self_consistency():
    """<H_uc> (the same unit-cell MultiOperator passed to
    set_hamiltonian, evaluated via the converged transfer-matrix
    machinery) must equal n_uc*e0_density by translational invariance --
    a model-agnostic identity, true regardless of the specific model.

    Uses a moderately dimerized (gapped) chain rather than the uniform
    (gapless/critical) Heisenberg chain, and forces a fixed, generous
    number of macro-iterations via a near-zero etol rather than relying
    on the energy-density finite-difference stopping criterion: the
    *energy* converges quadratically in the wavefunction error and so
    reports "converged" well before IDMRGResult.U_list (built from the
    growing algorithm's raw per-iteration tensors, not yet a fully
    gauge-fixed periodic MPS) has settled enough for a *linear*-order
    quantity like this correlator identity to match tightly -- confirmed
    directly: the identical model/maxm with etol=1e-9's early stop left a
    ~1.8e-2 gap, while forcing the same run out to a fixed 150
    macro-iterations closed it to ~4e-4. Also confirmed directly on the
    exactly solvable fully-decoupled-dimer limit (J_weak=0) that the
    correlator formulas themselves are exact (both intra-cell and
    inter-cell/wraparound correlators matched -0.75 and 0 to machine
    precision there), ruling out a correlator-formula bug as the cause of
    the gap this test tolerates. `idmrg_ground_state`'s local 2-site solve
    is now warm-started from the previous macro-iteration's own local
    ground vector at the same unit-cell position (see
    `_local_two_site_solve`'s own docstring) specifically to damp this kind
    of run-to-run drift, but this *dimerized* (gapped) model was never the
    badly-broken case that fix targets -- confirmed directly, post-fix gaps
    here (3 trials: 0.004, 0.004, 0.008) land squarely inside the same
    ~3e-3 to ~2e-2 range observed before the fix, not dramatically tighter
    -- so the tolerance below is left as-is rather than over-tightened on
    thin evidence. The *uniform* (undimerized) n_uc=2 chain is where the
    fix's effect is large and reproducible, see
    `test_n_uc2_uniform_expectation_self_consistency` below.

    maxm=8 (well below this model's natural, cutoff-set bond dimension of
    ~11-12) is deliberately small: it makes truncation land exactly on
    maxm on *both* the intra- and inter-cell bonds, avoiding a separate,
    real edge case where two independent SVDs (one per micro-step of the
    unit cell) round a borderline singular value across the cutoff
    threshold differently for what is, by periodicity, the *same*
    wraparound bond -- confirmed to happen intermittently at maxm=40
    (natural dimension 11 vs 12, one singular value short), which
    idmrg.py's _dominant_right_fixed_point now raises a clear
    RuntimeError for instead of a cryptic reshape crash (see its own
    comment) -- a real, documented limitation of reusing the growing
    algorithm's raw per-iteration U_list as a stand-in unit cell, not
    something this test should also have to handle."""
    n_uc = 2
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    ic.gs_method = "idmrg"
    j_strong, j_weak = 1.0, 0.2
    h = (j_strong * (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1])
         + j_weak * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0]))
    ic.maxm = 8
    ic.maxiter = 150
    ic.etol = 1e-14  # forces the full 150 iterations, see docstring above
    ic.set_hamiltonian(h)
    ic.gs_energy()
    # <H_uc> = sum over every bond touched by the unit cell (intra-cell
    # 0-1, plus the inter-cell n_uc-1 -> 0-of-next-cell), via
    # correlator(name, p_i, name, r=1) = <name(p_i) name(p_i+1)>.
    total = 0.0
    for i in range(n_uc):
        for name in ("Sx", "Sy", "Sz"):
            total += ic.correlator(name, i, name, 1).real
    assert total == pytest.approx(n_uc * ic.e0, abs=5e-2)


def test_n_uc2_uniform_expectation_self_consistency():
    """Same <H_uc> self-consistency identity as
    test_unit_cell_expectation_self_consistency, but for the *uniform*
    (undimerized) n_uc=2 Heisenberg chain -- the regime where
    `idmrg_ground_state`'s warm-started local solve (see
    `_local_two_site_solve`'s own docstring) makes a large, reproducible
    difference. Before that fix, this exact configuration gave gaps of
    0.338-0.398 across independent runs (uniformly bad, not noise around a
    correct value -- confirmed by a best-of-6-seeds experiment that still
    landed 0.38-0.60 every time). After the fix, 7 independent trials all
    landed at <=0.0511 (several essentially exact, ~2e-6), a qualitative
    change, not just a quantitative improvement -- `abs=0.1` below sits
    comfortably above every trial observed post-fix while remaining well
    below the pre-fix range, so this test would fail on the old,
    always-fresh-random-restart behavior.

    Unlike the dimerized test above, `n_uc=1` is deliberately *not* also
    covered here: the same warm-start fix does not resolve n_uc=1's own
    analogous gap (confirmed directly: it stays in the 0.4-1.8 range, and
    IDMRGResult.state_overlap plateaus around 0.5-0.65 rather than trending
    toward 1 even after 80 macro-iterations) -- a distinct, still-open
    limitation isolated to n_uc=1's structurally unique micro-step (p_L
    always equals p_R there, unlike n_uc=2 where they never coincide; see
    idmrg.py's own `_fresh_physical_copy`)."""
    n_uc = 2
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    ic.gs_method = "idmrg"
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0])
    ic.maxm = 18
    ic.maxiter = 100
    ic.etol = 1e-14  # forces the full 100 iterations
    ic.set_hamiltonian(h)
    ic.gs_energy()
    total = 0.0
    for i in range(n_uc):
        for name in ("Sx", "Sy", "Sz"):
            total += ic.correlator(name, i, name, 1).real
    assert total == pytest.approx(n_uc * ic.e0, abs=0.1)


def test_n_uc2_dimerized_chain_matches_finite_size_extrapolation():
    """A dimerized (alternating strong/weak bond) n_uc=2 chain, cross-
    checked against a finite-size ED extrapolation: two different large
    open chains' total ground-state energies, differenced and divided by
    the site-count difference, cancel boundary effects to leading order
    and converge to the bulk energy density."""
    j_strong, j_weak = 1.0, 0.4

    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    ic.gs_method = "idmrg"
    h = (j_strong * (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1])
         + j_weak * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0]))
    ic.maxm = 40
    ic.maxiter = 200
    ic.etol = 1e-9
    ic.set_hamiltonian(h)
    density = ic.gs_energy()
    assert ic.converged

    def finite_energy(n_sites):
        sc = spinchain.Spin_Chain(["1/2"] * n_sites)
        hf = 0
        for i in range(n_sites - 1):
            j = j_strong if i % 2 == 0 else j_weak
            hf = hf + j * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                           + sc.Sz[i] * sc.Sz[i + 1])
        sc.set_hamiltonian(hf)
        return sc.gs_energy(mode="ED")

    n1, n2 = 12, 16
    e1, e2 = finite_energy(n1), finite_energy(n2)
    extrapolated_density = (e2 - e1) / (n2 - n1)
    assert density == pytest.approx(extrapolated_density, abs=5e-3)


def test_n_uc2_skip_site_coupling_via_R_suffix():
    """A coupling from site 0 of the central cell directly to site 0 of
    the next cell, skipping over site 1 -- exercises the R-suffix
    canonicalization path a strict nearest-neighbor-only design couldn't
    express. Only checks this runs and converges to *some* sane (real,
    finite, negative) energy density -- there's no simple closed form
    for this model, this is a smoke/regression test for the R-suffix
    validation and automaton-construction path, not a numerical oracle
    check."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    ic.gs_method = "idmrg"
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + 0.5 * (ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]))
    ic.maxm = 30
    ic.maxiter = 150
    ic.etol = 1e-8
    ic.set_hamiltonian(h)
    e0 = ic.gs_energy()
    assert np.isfinite(e0)
    assert e0 < 0


def test_n_uc_above_2_is_constructible_and_idmrg_still_rejects_it():
    """Construction no longer rejects a big cell: gs_method="vumps" (the
    default) handles any n_uc via the sequential multi-site algorithm. The
    restriction now lives with the algorithm that actually has it -- the
    growing algorithm's micro-step pairing, which is only adjacent for
    n_uc<=2 -- so gs_method="idmrg" must still refuse."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2", "1/2"])
    assert ic.n_uc == 3
    ic.gs_method = "idmrg"
    ic.set_hamiltonian(ic.SzC[0] * ic.SzC[1] + ic.SzC[1] * ic.SzC[2]
                        + ic.SzC[2] * ic.SzR[0])
    with pytest.raises(NotImplementedError):
        ic.gs_energy()


def test_itensor_version_other_than_python_or_3_not_implemented():
    """itensor_version=3 (the ITensor v3 C++ backend) is implemented -- see
    the itensor_version3_* tests below; 2 (no mpscpp2 iDMRG port) and
    "julia_live" (no mpsjulialive iDMRG port either) are not."""
    with pytest.raises(NotImplementedError):
        infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=2)
    with pytest.raises(NotImplementedError):
        infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="julia_live")


@pytest.mark.skipif(not cppext.available(3),
                     reason="requires the compiled mpscpp3 (ITensor v3) extension")
@pytest.mark.parametrize("n_uc", [1, 2])
def test_itensor_version3_energy_density_matches_bethe_ansatz(n_uc):
    """Lighter-weight than the itensor_version="python" versions of this
    check above (maxm=30, maxiter=60 rather than maxm=40, maxiter=200):
    the v3 C++ backend is ~2.5-3x slower than pyitensor at matched
    parameters (confirmed independently), so this uses a smaller
    configuration, cross-checked directly against pyitensor at the same
    parameters in test_itensor_version3_matches_python_backend below."""
    ic, _ = _converged_uniform_chain(n_uc, maxm=30, maxiter=60, etol=1e-9,
                                      itensor_version=3)
    assert ic.e0 == pytest.approx(EXACT_HEISENBERG_DENSITY, abs=1e-3)


@pytest.mark.skipif(not cppext.available(3),
                     reason="requires the compiled mpscpp3 (ITensor v3) extension")
def test_itensor_version3_matches_python_backend():
    """Cross-backend agreement at identical (maxm, maxiter, etol, niter):
    both solvers should land on the same truncated-bond-dimension energy
    density, not just agree in the maxm->infinity limit. Confirmed directly
    to agree to ~1e-8 (well inside the 1e-6 tolerance here) and to be stable
    run over run to ~1e-8 despite iDMRG's unseeded random starting MPS.

    Unlike the energy-only version this test used to be, it now also
    compares the converged *state* -- vev and a short correlator sweep --
    which is what actually guards the three things the C++ port was long
    missing relative to pyitensor/idmrg.py (McCulloch's wavefunction
    prediction, the theta-cell unit-cell extraction, and the per-site
    energy baseline; all three are ported now, see
    Chain::idmrg_ground_state's own comment). The energy density alone
    could not: the first two are state/gauge quality and leave it alone,
    and the third is exactly compensated on both sides.

    The state comparison uses a *dimerized* (gapped) chain rather than the
    uniform critical one the energy check above uses. On a gapless chain
    at finite maxm the growing algorithm settles onto a slightly,
    arbitrarily symmetry-broken state whose residual depends on its own
    unseeded random starting MPS -- measured at <Sz> ~ +4e-4 on one
    backend against ~ -3e-4 on the other here, both legitimate
    approximations to the exact 0, but not equal to each other. That is a
    property of the model, not of either implementation; comparing two
    backends' states needs a model whose state is unique. (The correlators
    would in fact have squeaked through, since the residual enters them
    quadratically -- which is exactly why the vev comparison is the one
    worth having.)"""
    ic_py, _ = _converged_uniform_chain(2, maxm=30, maxiter=30, etol=1e-9,
                                         itensor_version="python")
    ic_v3, _ = _converged_uniform_chain(2, maxm=30, maxiter=30, etol=1e-9,
                                         itensor_version=3)
    assert ic_v3.e0 == pytest.approx(ic_py.e0, abs=1e-6)

    dim_py = _dimerized_chain_for_state_check("python")
    dim_v3 = _dimerized_chain_for_state_check(3)
    assert dim_v3.e0 == pytest.approx(dim_py.e0, abs=1e-6)
    for p in (0, 1):
        assert dim_v3.vev("Sz", p) == pytest.approx(dim_py.vev("Sz", p), abs=1e-5)
    for r in (0, 1, 2, 3):
        assert (dim_v3.correlator("Sz", 0, "Sz", r)
                == pytest.approx(dim_py.correlator("Sz", 0, "Sz", r), abs=1e-5))


def _dimerized_chain_for_state_check(itensor_version, j_strong=1.0, j_weak=0.2):
    """A gapped, moderately dimerized n_uc=2 Heisenberg chain -- the same
    model test_unit_cell_expectation_self_consistency uses, including its
    etol=1e-14 "run the full maxiter" trick (see that test's docstring)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"],
                                            itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    h = (j_strong * (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1]
                     + ic.SzC[0] * ic.SzC[1])
         + j_weak * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0]
                     + ic.SzC[1] * ic.SzR[0]))
    ic.maxm, ic.maxiter, ic.etol = 8, 150, 1e-14
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


@pytest.mark.skipif(not cppext.available(3),
                     reason="requires the compiled mpscpp3 (ITensor v3) extension")
def test_itensor_version3_vev_and_correlator_rejects_bad_arguments():
    """vev/correlator ARE implemented on the v3 backend under
    gs_method="idmrg" now (see test_idmrg_correlator_v3.py for the value
    checks) -- what must still be rejected is out-of-range input."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "idmrg"
    h = ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]
    ic.set_hamiltonian(h)
    with pytest.raises(ValueError):
        ic.vev("Sz", 1)
    with pytest.raises(ValueError):
        ic.correlator("Sz", 0, "Sz", -1)


def test_l_and_r_in_same_term_rejected():
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    with pytest.raises(ValueError):
        ic.set_hamiltonian(ic.SxL[0] * ic.SxR[0])


# -- idmrg.apply_mpo: applying a periodic (bounded) MPO to the converged
# iMPS -- see idmrg.py's own "Applying a (bounded) MPO to the converged
# iMPS" section docstring for the algorithm and its scope restriction
# (bounded/local operators only, not idmrg.py's own extensive Hamiltonian
# automaton). These go through pyitensor.idmrg directly (not
# infinitechain.py, which has no public wrapper for this yet) --
# tests/test_metts_vev.py and others already establish that reaching into
# dmrgpy.pyitensor submodules directly from a test is a normal pattern
# here, not a workaround.

def _identity_mpo(sites_uc, n_uc):
    """chi_W=1 automaton: Id at every site, sharing sites_uc's own physical
    Indices (required, see grow_by_mpo's docstring)."""
    W = []
    for p in range(n_uc):
        d = sites_uc.dim(p + 1)
        s = sites_uc.si(p + 1)
        link_l = Index(1, tags="Link")
        link_r = Index(1, tags="Link")
        arr = np.eye(d, dtype=complex).reshape(1, d, d, 1)
        W.append(ITensor((link_l, s, s.prime(1), link_r), arr))
    # wraparound: every site's own fresh links are fine independently
    # (chi_W=1 throughout, so there is nothing to actually connect).
    return W


def _single_site_operator_mpo(sites_uc, n_uc, opname, p_active):
    """chi_W=1 automaton: `opname`'s matrix at sublattice p_active, Id
    everywhere else -- the simplest genuinely non-trivial bounded MPO."""
    W = []
    for p in range(n_uc):
        d = sites_uc.dim(p + 1)
        s = sites_uc.si(p + 1)
        mat = sites_uc.site_type(p + 1).matrix(opname) if p == p_active else np.eye(d, dtype=complex)
        link_l = Index(1, tags="Link")
        link_r = Index(1, tags="Link")
        W.append(ITensor((link_l, s, s.prime(1), link_r), mat.reshape(1, d, d, 1)))
    return W


@pytest.mark.parametrize("n_uc", [1, 2])
def test_apply_mpo_identity_is_noop(n_uc):
    """Applying the identity MPO must reproduce every existing observable
    to numerical precision -- isolates the canonicalization machinery
    (grow_by_mpo + _canonicalize_periodic) from any MPO-construction
    question, since chi_W=1 Id trivially can't grow the bond dimension."""
    ic, _ = _converged_uniform_chain(n_uc, maxm=20)
    W = _identity_mpo(ic._result.sites_uc, n_uc)
    new_result = idmrg.apply_mpo(ic._result, W, cutoff=1e-12, maxdim=None)
    assert new_result.eta == pytest.approx(1.0, abs=1e-8)
    for p in range(n_uc):
        orig = idmrg.onsite_expectation(ic._result, "Sz", p)
        new = idmrg.onsite_expectation(new_result, "Sz", p)
        assert new == pytest.approx(orig, abs=1e-8)
        for r in range(3):
            orig_c = idmrg.two_point_correlator(ic._result, "Sz", p, "Sz", r)
            new_c = idmrg.two_point_correlator(new_result, "Sz", p, "Sz", r)
            assert new_c == pytest.approx(orig_c, abs=1e-8)


def test_apply_mpo_pauli_x_flips_sz():
    """Pauli-X (2*Sx for spin-1/2) at every site is unitary and chi_W=1 --
    applying it must flip <Sz> -> -<Sz> exactly and leave <Sz(0)Sz(r)>
    unchanged exactly (two flipped operators multiply back to the
    original sign), a strong, closed-form check unrelated to any
    numerical oracle."""
    ic, _ = _converged_uniform_chain(1, maxm=20)
    sites_uc = ic._result.sites_uc
    d = sites_uc.dim(1)
    pauli_x = 2 * sites_uc.site_type(1).matrix("Sx")
    assert np.allclose(pauli_x @ pauli_x.conj().T, np.eye(d), atol=1e-10)
    link_l, link_r = Index(1, tags="Link"), Index(1, tags="Link")
    s = sites_uc.si(1)
    W = [ITensor((link_l, s, s.prime(1), link_r), pauli_x.reshape(1, d, d, 1))]

    new_result = idmrg.apply_mpo(ic._result, W, cutoff=1e-12, maxdim=None)
    assert new_result.eta == pytest.approx(1.0, abs=1e-8)

    orig_sz = idmrg.onsite_expectation(ic._result, "Sz", 0)
    new_sz = idmrg.onsite_expectation(new_result, "Sz", 0)
    assert new_sz == pytest.approx(-orig_sz, abs=1e-8)
    for r in range(1, 3):
        orig_c = idmrg.two_point_correlator(ic._result, "Sz", 0, "Sz", r)
        new_c = idmrg.two_point_correlator(new_result, "Sz", 0, "Sz", r)
        assert new_c == pytest.approx(orig_c, abs=1e-8)


def test_apply_mpo_two_site_gate_preserves_norm():
    """A genuinely bond-dimension>1 local gate (an SVD-split 2-site
    unitary rotation, tiled once per unit cell on the intra-cell bond
    only -- identity on the inter-cell bond, so this is a bounded/local
    operator, not idmrg.py's own extensive Hamiltonian automaton) --
    exercises the actual bond-growth path (grow_by_mpo produces a raw,
    non-canonical bond) unlike the chi_W=1 tests above, which barely
    touch _canonicalize_periodic's fixed-point machinery. Being unitary,
    it must preserve the norm (apply_mpo's own `eta` diagnostic ~1);
    _canonicalize_periodic's internal left-canonicality check (raises
    RuntimeError if it fails) is the other half of this test, implicitly
    exercised by apply_mpo not raising."""
    from scipy.linalg import expm

    n_uc = 2
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    ic.gs_method = "idmrg"
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + 0.4 * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0]))
    ic.maxm = 4
    ic.maxiter = 100
    ic.etol = 1e-12
    ic.set_hamiltonian(h)
    ic.gs_energy()

    sites_uc = ic._result.sites_uc
    d = sites_uc.dim(1)
    Sx = sites_uc.site_type(1).matrix("Sx")
    Sy = sites_uc.site_type(1).matrix("Sy")
    Sz = sites_uc.site_type(1).matrix("Sz")
    H2 = (np.kron(Sx, Sx) + np.kron(Sy, Sy) + np.kron(Sz, Sz)).real
    gate = expm(-1j * 0.37 * H2)  # d^2 x d^2 unitary 2-site rotation
    gate4 = np.transpose(gate.reshape(d, d, d, d), (2, 0, 3, 1))  # (s0,s0',s1,s1')
    U, S, Vh = np.linalg.svd(gate4.reshape(d * d, d * d), full_matrices=False)
    keep = int(np.sum(S > 1e-12))
    U, S, Vh = U[:, :keep], S[:keep], Vh[:keep, :]
    a_half = (U * S[None, :] ** 0.5).reshape(d, d, keep)
    b_half = (S[:, None] ** 0.5 * Vh).reshape(keep, d, d)

    s0, s1 = sites_uc.si(1), sites_uc.si(2)
    left_dummy, mid, right_dummy = (Index(1, tags="Link"), Index(keep, tags="Link"),
                                     Index(1, tags="Link"))
    W0 = ITensor((left_dummy, s0, s0.prime(1), mid), a_half.reshape(1, d, d, keep))
    W1 = ITensor((mid, s1, s1.prime(1), right_dummy), b_half.reshape(keep, d, d, 1))

    new_result = idmrg.apply_mpo(ic._result, [W0, W1], cutoff=1e-10, maxdim=None)
    assert new_result.eta == pytest.approx(1.0, abs=1e-6)


def test_grow_by_mpo_rejects_mismatched_length():
    ic, _ = _converged_uniform_chain(2, maxm=8)
    W = _identity_mpo(ic._result.sites_uc, 2)
    with pytest.raises(ValueError):
        idmrg.grow_by_mpo(W[:1], ic._result.cell_list, 2)


# -- idmrg.imps_overlap: per-site overlap/fidelity between two converged
# infinite MPS -- see idmrg.py's own docstring for the mixed-transfer-matrix
# construction and the physical meaning of the normalized (per-site
# fidelity) vs raw (normalize=False) return value.

@pytest.mark.parametrize("n_uc", [1, 2])
def test_imps_overlap_self_is_one(n_uc):
    """<psi|psi> per site must be exactly 1 (up to numerical precision),
    both normalized (the default) and raw (normalize=False, since a
    converged IDMRGResult.U_list is already left-canonical, so its own
    self-overlap transfer eigenvalue is already 1 without needing the
    normalization division)."""
    ic, _ = _converged_uniform_chain(n_uc, maxm=20)
    ov = idmrg.imps_overlap(ic._result, ic._result)
    assert ov == pytest.approx(1.0, abs=1e-8)
    ov_raw = idmrg.imps_overlap(ic._result, ic._result, normalize=False)
    assert ov_raw == pytest.approx(1.0, abs=1e-8)


@pytest.mark.parametrize("n_uc", [1, 2])
def test_imps_overlap_identity_mpo_preserves_state(n_uc):
    """Applying the identity MPO must reproduce the same physical state up
    to gauge -- |imps_overlap| between the original and the identity-MPO'd
    copy must be exactly 1, a gauge-independent cross-check complementing
    test_apply_mpo_identity_is_noop's own per-observable comparison."""
    ic, _ = _converged_uniform_chain(n_uc, maxm=20)
    W = _identity_mpo(ic._result.sites_uc, n_uc)
    new_result = idmrg.apply_mpo(ic._result, W, cutoff=1e-12, maxdim=None)
    ov = idmrg.imps_overlap(ic._result, new_result)
    assert abs(ov) == pytest.approx(1.0, abs=1e-6)


def test_imps_overlap_orthogonal_polarized_states_near_zero():
    """Two independently-converged n_uc=1 chains with a strong onsite field
    pinned to opposite Sz polarization converge to (near-exact) product
    states |up> and |down>, which are exactly orthogonal -- |imps_overlap|
    must be ~0, regardless of the fact that the two chains were never
    built or grown together. This exercises the same onsite-field
    machinery whose real bug (see _build_periodic_mpo's own docstring) was
    found and fixed while writing this test."""
    ic_up = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic_up.gs_method = "idmrg"
    ic_up.maxm, ic_up.maxiter, ic_up.etol = 4, 50, 1e-12
    ic_up.set_hamiltonian(-10.0 * ic_up.SzC[0])
    ic_up.gs_energy()

    ic_down = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic_down.gs_method = "idmrg"
    ic_down.maxm, ic_down.maxiter, ic_down.etol = 4, 50, 1e-12
    ic_down.set_hamiltonian(10.0 * ic_down.SzC[0])
    ic_down.gs_energy()

    ov = idmrg.imps_overlap(ic_up._result, ic_down._result)
    assert abs(ov) == pytest.approx(0.0, abs=1e-8)


def test_imps_overlap_rejects_n_uc_mismatch():
    ic1, _ = _converged_uniform_chain(1, maxm=8)
    ic2, _ = _converged_uniform_chain(2, maxm=8)
    with pytest.raises(ValueError):
        idmrg.imps_overlap(ic1._result, ic2._result)


def test_imps_overlap_rejects_dimension_mismatch():
    """Different local Hilbert-space dimensions per sublattice (spin-1/2
    vs spin-1) is a meaningless overlap -- must raise, not silently
    contract mismatched-size physical legs."""
    ic_half = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic_half.gs_method = "idmrg"
    ic_half.maxm, ic_half.maxiter, ic_half.etol = 4, 20, 1e-8
    ic_half.set_hamiltonian(ic_half.SxC[0] * ic_half.SxR[0]
                             + ic_half.SyC[0] * ic_half.SyR[0]
                             + ic_half.SzC[0] * ic_half.SzR[0])
    ic_half.gs_energy()

    ic_one = infinitechain.Infinite_Spin_Chain(["1"])
    ic_one.gs_method = "idmrg"
    ic_one.maxm, ic_one.maxiter, ic_one.etol = 4, 20, 1e-8
    ic_one.set_hamiltonian(ic_one.SxC[0] * ic_one.SxR[0]
                            + ic_one.SyC[0] * ic_one.SyR[0]
                            + ic_one.SzC[0] * ic_one.SzR[0])
    ic_one.gs_energy()

    with pytest.raises(ValueError):
        idmrg.imps_overlap(ic_half._result, ic_one._result)


# -- Regression tests for a real bug found while writing the imps_overlap
# tests above: _build_periodic_mpo wired a Hamiltonian's onsite terms onto
# the automaton's F,F self-loop instead of a direct S,F transition -- see
# that function's own docstring for the full derivation. This silently
# dropped onsite terms entirely for any bond-less Hamiltonian, and caused
# an exponential energy blow-up once at least one bond term was also
# present (both confirmed directly before the fix).

def test_onsite_field_matches_exact_solution():
    """A pure onsite (Zeeman-field) Hamiltonian with no bond term at all is
    exactly solvable (a decoupled product state at every site): the ground
    state aligns <Sz> with sign(B) (to minimize H=-B*Sz), so e0 must be
    -|B|/2 and <Sz> must be sign(B)*0.5 exactly, for either field
    direction. Before the fix, e0 stayed exactly 0 and <Sz> ~ 0 for every
    field strength tried, since W stayed block-diagonal between S and F
    with no bond term ever able to activate F -- confirmed directly."""
    for B in (5.0, -5.0, 2.0):
        ic = infinitechain.Infinite_Spin_Chain(["1/2"])
        ic.gs_method = "idmrg"
        ic.maxm, ic.maxiter, ic.etol = 4, 50, 1e-12
        ic.set_hamiltonian(-B * ic.SzC[0])
        e0 = ic.gs_energy()
        assert e0 == pytest.approx(-abs(B) / 2, abs=1e-6)
        sz = ic.vev("Sz", 0)
        assert sz.real == pytest.approx(0.5 if B > 0 else -0.5, abs=1e-6)


def test_onsite_field_with_bonds_does_not_diverge():
    """A small field added to an otherwise-normal (bond-coupled,
    dimerized) n_uc=2 chain must converge to a finite energy density --
    before the fix, once the bond term activated the automaton's F
    channel, every further absorbed site re-added the field content on top
    of the already-accumulated total (multiplicative, not additive):
    confirmed directly, energy density reached -1e23 at B=0.3 and -1e69 at
    B=1.0, with `.converged` staying False throughout. This dimerized
    model is gapped, so a field well below its gap should leave the energy
    density unchanged from B=0 (confirmed reproducible to ~1e-13 run over
    run at this exact configuration) -- a much stronger check than merely
    "finite"."""
    def density_at(B):
        ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
        ic.gs_method = "idmrg"
        h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
             + 0.5 * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0])
             - B * (ic.SzC[0] + ic.SzC[1]))
        ic.maxm, ic.maxiter, ic.etol = 20, 100, 1e-11
        ic.set_hamiltonian(h)
        e0 = ic.gs_energy()
        assert ic.converged
        return e0

    d0 = density_at(0.0)
    d1 = density_at(0.2)
    assert np.isfinite(d1) and abs(d1) < 10  # rules out the 1e23-scale blow-up
    assert d1 == pytest.approx(d0, abs=1e-6)


def test_two_point_correlator_same_site_composition_order():
    """two_point_correlator's r=0 branch must compose opname_i/opname_j in
    the same order idmrg.py's own onsite-term automaton builder
    (_classify_terms/_term_site_matrices, and through it, checked directly
    against HTerm.resolve()/MultiOperator.multiply_MO() -- see
    _term_site_matrices' own docstring) uses for an equivalent
    "opname_i[p]*opname_j[p]" term elsewhere in dmrgpy. An earlier version
    of this branch composed the opposite order (M_i @ M_j instead of
    M_j @ M_i, both in the stored (in,out) convention) -- invisible to
    every other test in this file, which all pass the *same* operator name
    twice (e.g. "Sz","Sz"), unable to distinguish the two orders for a
    commuting pair.

    Ground state: a strong field along Y gives an exact, known
    |Sy=+1/2> product state (same trick as
    test_onsite_field_matches_exact_solution) -- for this state <Sx*Sz>
    and <Sz*Sx> differ (their difference is the commutator expectation
    <[Sx,Sz]> = -i<Sy> != 0), so this is a genuine, order-sensitive check,
    with the expected value computed independently via
    idmrg._term_site_matrices (not via two_point_correlator itself)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "idmrg"
    ic.maxm, ic.maxiter, ic.etol = 4, 50, 1e-12
    ic.set_hamiltonian(-5.0 * ic.SyC[0])
    ic.gs_energy()
    result = ic._result

    op_term = [1.0, ["Sx", 0], ["Sz", 0]]
    _, _, mats, _ = idmrg._term_site_matrices(op_term, result.sites_uc, 1)
    M = mats[0]
    cell, n_cell = idmrg._correlator_cell(result)
    Es = idmrg._transfer_matrices(cell, n_cell)
    rho_after, _ = idmrg._all_right_fixed_points(Es, n_cell)
    A = idmrg._to_array_lpr(cell[0])
    Aop = np.einsum('io,lir->lor', M, A)
    E4 = np.einsum('lpr,LpR->lLrR', Aop, np.conj(A))
    expected = idmrg._expectation(Es, E4, Es[0], rho_after[0], 0, n_cell)

    got = idmrg.two_point_correlator(result, "Sx", 0, "Sz", 0)
    assert got == pytest.approx(expected, abs=1e-8)
    swapped = idmrg.two_point_correlator(result, "Sz", 0, "Sx", 0)
    assert abs(got - swapped) > 1e-3  # confirms this case is order-sensitive


# -- idmrg.imps_sum: direct sum of two converged infinite MPS -- see
# idmrg.py's own "Summing two converged iMPS" section docstring for the
# construction (periodic-chain analogue of mpsalgebra.sum, generalized
# from a finite chain's dimension-1 boundaries to a fully block-diagonal
# per-cut construction) and its physical scope: tiled to the
# thermodynamic limit, this only has a well-posed single-branch
# reduction when the two summands' own self-overlap eigenvalues (eta)
# differ -- summing two *ordinary* IDMRGResults (always eta=1 each, by
# left-canonical SVD construction) hits a genuine tie every time, which
# idmrg._dominant_right_fixed_point's own degeneracy check must catch
# and raise on, not silently resolve to one arbitrary branch.

def test_imps_sum_tied_norm_raises():
    """The common case: summing two separately-converged, ordinary
    IDMRGResults. Both are individually normalized to eta=1 exactly (SVD
    construction), so the combined transfer matrix's dominant eigenvalue
    is exactly degenerate -- idmrg._dominant_right_fixed_point's
    degeneracy check must raise RuntimeError here rather than silently
    collapsing to one of the two (physically distinct, oppositely
    Sz-polarized) branches."""
    ic_up = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic_up.gs_method = "idmrg"
    ic_up.maxm, ic_up.maxiter, ic_up.etol = 4, 50, 1e-12
    ic_up.set_hamiltonian(-10.0 * ic_up.SzC[0])
    ic_up.gs_energy()

    ic_down = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic_down.gs_method = "idmrg"
    ic_down.maxm, ic_down.maxiter, ic_down.etol = 4, 50, 1e-12
    ic_down.set_hamiltonian(10.0 * ic_down.SzC[0])
    ic_down.gs_energy()

    with pytest.raises(RuntimeError):
        idmrg.imps_sum(ic_up._result, ic_down._result)


@pytest.mark.parametrize("n_uc", [1, 2])
def test_imps_sum_dominant_branch_survives(n_uc):
    """A well-posed (non-degenerate) case: sum a converged state with a
    deliberately norm-rescaled copy of a *different* converged state
    (amplitude x0.9, i.e. self-overlap eta=0.81 rather than the ordinary
    1) -- an artificial construction (ordinary IDMRGResults never carry
    eta!=1 on their own), but the only way to exercise the "well-posed"
    branch of imps_sum's own scope at all, since two literal IDMRGResults
    always tie (see test_imps_sum_tied_norm_raises above). The
    smaller-norm branch must be exponentially suppressed: every static
    observable of the summed-and-truncated result must reproduce the
    larger-norm (unscaled) branch's own value exactly, and the surviving
    bond dimension must collapse back down to that branch's own (i.e. the
    rescaled branch is fully discarded, not merely down-weighted)."""
    ic_dom, _ = _converged_uniform_chain(n_uc, maxm=16, maxiter=60, etol=1e-11)

    spins = ["1/2"] * n_uc
    ic_other = infinitechain.Infinite_Spin_Chain(spins)
    ic_other.gs_method = "idmrg"
    # A different (XXZ-anisotropic, Delta=0.3) Hamiltonian so the two
    # states are genuinely distinct, not just re-derived copies of the
    # same ic_dom (isotropic Heisenberg) state -- works identically for
    # n_uc=1 and n_uc=2, unlike a dimerization (which needs n_uc>=2).
    terms = []
    for i in range(n_uc - 1):
        terms.append(ic_other.SxC[i] * ic_other.SxC[i + 1] + ic_other.SyC[i] * ic_other.SyC[i + 1]
                     + 0.3 * ic_other.SzC[i] * ic_other.SzC[i + 1])
    terms.append(ic_other.SxC[n_uc - 1] * ic_other.SxR[0] + ic_other.SyC[n_uc - 1] * ic_other.SyR[0]
                 + 0.3 * ic_other.SzC[n_uc - 1] * ic_other.SzR[0])
    h_other = terms[0]
    for t in terms[1:]:
        h_other = h_other + t
    ic_other.maxm, ic_other.maxiter, ic_other.etol = 16, 60, 1e-11
    ic_other.set_hamiltonian(h_other)
    ic_other.gs_energy()

    # Built from `cell_list`, not `U_list`: the raw per-micro-step U_list is
    # not a gauge-consistent periodic MPS (see idmrg.py's IDMRGResult
    # docstring), and every consumer here tiles the cell instead.
    scaled_other = idmrg.PeriodicMPS(
        ic_other._result.sites_uc, n_uc,
        [ITensor(T.inds, T.array * 0.9) for T in ic_other._result.cell_list],
        eta=0.81)

    result = idmrg.imps_sum(ic_dom._result, scaled_other, cutoff=1e-12, maxdim=None)
    assert result.eta == pytest.approx(1.0, abs=1e-6)
    for p in range(n_uc):
        orig = idmrg.onsite_expectation(ic_dom._result, "Sz", p)
        new = idmrg.onsite_expectation(result, "Sz", p)
        assert new == pytest.approx(orig, abs=1e-6)
        chi_dom = idmrg._to_array_lpr(ic_dom._result.cell_list[p]).shape[0]
        chi_new = idmrg._to_array_lpr(result.U_list[p]).shape[0]
        assert chi_new == chi_dom


def test_imps_sum_rejects_n_uc_mismatch():
    ic1, _ = _converged_uniform_chain(1, maxm=8)
    ic2, _ = _converged_uniform_chain(2, maxm=8)
    with pytest.raises(ValueError):
        idmrg.imps_sum(ic1._result, ic2._result)


def test_imps_sum_rejects_dimension_mismatch():
    """Mirrors imps_overlap's own dimension-mismatch check: different
    local Hilbert-space dimensions per sublattice (spin-1/2 vs spin-1) is
    a meaningless sum -- must raise, not silently place mismatched-size
    physical legs into the same combined tensor."""
    ic_half = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic_half.gs_method = "idmrg"
    ic_half.maxm, ic_half.maxiter, ic_half.etol = 4, 20, 1e-8
    ic_half.set_hamiltonian(ic_half.SxC[0] * ic_half.SxR[0]
                             + ic_half.SyC[0] * ic_half.SyR[0]
                             + ic_half.SzC[0] * ic_half.SzR[0])
    ic_half.gs_energy()

    ic_one = infinitechain.Infinite_Spin_Chain(["1"])
    ic_one.gs_method = "idmrg"
    ic_one.maxm, ic_one.maxiter, ic_one.etol = 4, 20, 1e-8
    ic_one.set_hamiltonian(ic_one.SxC[0] * ic_one.SxR[0]
                            + ic_one.SyC[0] * ic_one.SyR[0]
                            + ic_one.SzC[0] * ic_one.SzR[0])
    ic_one.gs_energy()

    with pytest.raises(ValueError):
        idmrg.imps_sum(ic_half._result, ic_one._result)


def test_periodic_direct_sum_bond_dimension():
    """idmrg._periodic_direct_sum's own raw (pre-canonicalization) output:
    an exact, block-diagonal construction -- the combined bond dimension
    at every cut (including the wraparound) must be exactly the sum of
    the two inputs' own bond dimensions there, with each input's own
    block placed at the expected array offset (not just the right shape
    by coincidence)."""
    ic_up, _ = _converged_uniform_chain(1, maxm=4, maxiter=50, etol=1e-12)
    ic_down = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic_down.gs_method = "idmrg"
    ic_down.maxm, ic_down.maxiter, ic_down.etol = 4, 50, 1e-12
    ic_down.set_hamiltonian(10.0 * ic_down.SzC[0])
    ic_down.gs_energy()

    raw = idmrg._periodic_direct_sum(ic_up._result.cell_list,
                                      ic_down._result.cell_list, 2)
    arr = raw[0].array
    chi_up = idmrg._to_array_lpr(ic_up._result.cell_list[0]).shape[0]
    chi_down = idmrg._to_array_lpr(ic_down._result.cell_list[0]).shape[0]
    assert arr.shape == (chi_up + chi_down, 2, chi_up + chi_down)
    a_up = idmrg._to_array_lpr(ic_up._result.cell_list[0])
    a_down = idmrg._to_array_lpr(ic_down._result.cell_list[0])
    assert np.allclose(arr[:chi_up, :, :chi_up], a_up)
    assert np.allclose(arr[chi_up:, :, chi_up:], a_down)
    assert np.allclose(arr[:chi_up, :, chi_up:], 0)
    assert np.allclose(arr[chi_up:, :, :chi_up], 0)


# -- kpm_finite: finite-window KPM dynamical correlator --
# (infinitechain.py's _window_hamiltonian + Infinite_Many_Body_Chain.
# kpm_finite, reusing kpmdmrg.get_dynamical_correlator /
# pyitensor/chain.py's KPM machinery verbatim on a temporary finite,
# open-boundary Many_Body_Chain -- see kpm_finite's own docstring for why
# this is a finite-window *approximation*, not an exact infinite-size
# method, and named accordingly).

def test_window_hamiltonian_matches_hand_built_finite_chain():
    """_window_hamiltonian's own construction (tiling h_intra/h_inter
    across n_window unit cells) must reproduce *exactly* the same
    Hamiltonian as hand-building the equivalent finite open chain term by
    term (the same pattern test_n_uc2_dimerized_chain_matches_finite_size_
    extrapolation already cross-checks iDMRG's own energy density
    against) -- checked via exact ED ground-state energy agreement, which
    only holds if every term (site indices, coefficients) matches
    exactly, not just approximately."""
    j_strong, j_weak = 1.0, 0.4
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    h = (j_strong * (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1])
         + j_weak * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0]))
    ic.set_hamiltonian(h)

    n_window = 4  # 8 sites
    n_sites = n_window * ic.n_uc
    h_window = infinitechain._window_hamiltonian(ic._h_intra, ic._h_inter, ic.n_uc, n_window)

    sc_hand = spinchain.Spin_Chain(["1/2"] * n_sites)
    hf = 0
    for i in range(n_sites - 1):
        j = j_strong if i % 2 == 0 else j_weak
        hf = hf + j * (sc_hand.Sx[i] * sc_hand.Sx[i + 1] + sc_hand.Sy[i] * sc_hand.Sy[i + 1]
                       + sc_hand.Sz[i] * sc_hand.Sz[i + 1])
    sc_hand.set_hamiltonian(hf)
    e_hand = sc_hand.gs_energy(mode="ED")

    sc_window = spinchain.Spin_Chain(["1/2"] * n_sites)
    sc_window.set_hamiltonian(h_window)
    e_window = sc_window.gs_energy(mode="ED")

    assert e_window == pytest.approx(e_hand, abs=1e-10)


def test_kpm_finite_runs_and_is_finite():
    """Smoke test: a uniform n_uc=1 Heisenberg chain's local (r=0, same
    operator) dynamical correlator on a modest finite window runs without
    error and returns a finite, non-trivial (not identically zero)
    profile -- mirrors test_n_uc2_skip_site_coupling_via_R_suffix's own
    smoke-test style (no simple closed form to check against for this
    combination of window size/DMRG/KPM parameters, just that the whole
    finite-window-then-KPM pipeline produces sane output)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    h = ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]
    ic.set_hamiltonian(h)

    es, ys = ic.kpm_finite(
        "Sz", 0, "Sz", 0, n_window=10,
        window_chain_kwargs=dict(maxm=20, nsweeps=8),
        delta=0.3, es=np.linspace(-1, 5, 60))

    assert es.shape == ys.shape == (60,)
    assert np.all(np.isfinite(ys))
    assert np.max(np.abs(ys)) > 1e-6  # not identically zero


def test_kpm_finite_before_set_hamiltonian_raises():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    with pytest.raises(RuntimeError):
        ic.kpm_finite("Sz", 0, "Sz", 0, n_window=5)


def test_kpm_finite_rejects_p_i_out_of_range():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0])
    with pytest.raises(ValueError):
        ic.kpm_finite("Sz", 5, "Sz", 0, n_window=5)


def test_kpm_finite_rejects_negative_r():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0])
    with pytest.raises(ValueError):
        ic.kpm_finite("Sz", 0, "Sz", -1, n_window=5)


def test_kpm_finite_rejects_window_too_small_for_r():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0])
    with pytest.raises(ValueError):
        ic.kpm_finite("Sz", 0, "Sz", 10, n_window=3)


# -- td_dynamical_correlator: the actual infinite-boundary-condition
# (IBC), real-time-TDVP dynamical correlator (arXiv:1804.09163 Sec. V.1)
# -- see pyitensor/idmrg_window.py's own module docstring and
# tests/test_idmrg_window.py for the underlying construction's own
# (extensively validated) coverage; the tests here only exercise the
# public infinitechain.py wrapper itself (argument validation, gs_energy()
# auto-call, itensor_version gating), not re-derive that validation.

def test_td_dynamical_correlator_runs_and_is_finite():
    """Smoke test: a gapped n_uc=1 transverse-field-like model's S(k,omega)
    runs without error and returns a finite, non-trivial profile -- mirrors
    test_kpm_finite_runs_and_is_finite's own smoke-test style."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "idmrg"
    h = 1.4 * ic.SzC[0] + ic.SxC[0] * ic.SxR[0]
    ic.maxm = 20
    ic.maxiter = 300
    ic.etol = 1e-12
    ic.niter = 150
    ic.set_hamiltonian(h)

    ks, es, Skw = ic.td_dynamical_correlator(
        "Sz", 0, "Sz", n_window=10, dt=0.05, nt=20,
        maxdim=40, cutoff=1e-10, niter=50, x_values=range(-4, 5),
        ks=np.linspace(-np.pi, np.pi, 9), delta=0.15, window=[-1, 5])

    assert ks.shape == (9,)
    assert Skw.shape == (9, es.shape[0])
    assert np.all(np.isfinite(Skw))
    assert np.max(np.abs(Skw)) > 1e-6  # not identically zero


def test_td_dynamical_correlator_calls_gs_energy_automatically():
    """Unlike kpm_finite (no dependency on a converged IDMRGResult),
    td_dynamical_correlator needs env_HL/env_HR -- it should call
    gs_energy() itself if the caller hasn't already, exactly like
    vev/correlator do."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "idmrg"
    ic.maxm = 20
    ic.maxiter = 300
    ic.etol = 1e-12
    ic.niter = 150
    ic.set_hamiltonian(1.4 * ic.SzC[0] + ic.SxC[0] * ic.SxR[0])
    assert ic._result is None
    ic.td_dynamical_correlator("Sz", 0, "Sz", n_window=8, dt=0.05, nt=4,
                                maxdim=30, cutoff=1e-10, niter=50,
                                x_values=range(-2, 3),
                                ks=np.linspace(-np.pi, np.pi, 5))
    assert ic._result is not None


def test_td_dynamical_correlator_rejects_p_i_out_of_range():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0])
    with pytest.raises(ValueError):
        ic.td_dynamical_correlator("Sz", 5, "Sz", n_window=5)


@pytest.mark.skipif(not cppext.available(3),
                     reason="ITensor v3 extension not compiled")
def test_td_dynamical_correlator_runs_on_v3_backend():
    """itensor_version=3 must return an S(k,omega) whose underlying
    S(x,t=0) reproduces the chain's own static correlator exactly.

    Between two earlier versions of this test, that identity is the whole
    point. It first asserted only that the call returned a finite,
    correctly-shaped S(k,omega) -- which it did, with wrong numbers
    inside, because this backend's window tiled the raw per-micro-step
    idmrg_U_ factors rather than the gauge-consistent unit cell, and
    neither shape nor finiteness can see a gauge error. It then asserted
    the call REFUSED, while that was being fixed. Now it asserts the
    thing that can actually fail."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "idmrg"
    ic.maxm, ic.maxiter, ic.etol, ic.niter = 8, 30, 1e-6, 30
    # A plain single-coupling XX term leaves the transfer matrix's
    # dominant eigenvalue (near-)degenerate at this loose convergence
    # (confirmed directly: a RuntimeError out of the transfer-matrix
    # fixed-point solve) -- the full XXX Heisenberg coupling (matching
    # examples/idmrg/heisenberg_infinite_python_VS_v3/main.py and
    # tests/test_idmrg_window_v3.py's own model) doesn't have this
    # problem.
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0]
                        + ic.SzC[0] * ic.SzR[0])
    ks = np.linspace(-np.pi, np.pi, 3)
    kk, ee, Sk = ic.td_dynamical_correlator(
        "Sz", 0, "Sz", n_window=8, dt=0.05, nt=3, maxdim=20,
        cutoff=1e-10, niter=30, x_values=[-2, 0, 2], ks=ks)
    assert np.all(np.isfinite(Sk))
    # The oracle, one level down: S(x,0) == correlator(x), exactly. Even
    # separations only -- with n_uc=1 the extracted cell is still two
    # sites long, so an odd separation measured from the window's centre
    # starts from the other cell position than correlator's own p_i=0
    # (see tests/test_idmrg_window_v3.py's own note).
    _ts, xarr, S = ic._session3.td_dynamical_correlator_window(
        8, "Sz", "Sz", 0.05, 1, [-2, 0, 2], 20, 1e-10, 30, False, 0)
    S = np.array(S).reshape(1, len(xarr))
    for ix, x in enumerate(xarr):
        assert S[0, ix] == pytest.approx(
            complex(ic.correlator("Sz", 0, "Sz", abs(int(x)))), abs=1e-8)


def test_td_dynamical_correlator_agrees_qualitatively_with_kpm_finite():
    """Cross-check the two independent dynamical-correlator approximations
    against each other: kpm_finite (open-boundary window + KPM/Chebyshev)
    and td_dynamical_correlator (IBC window + real-time TDVP) are
    different approximation schemes with different systematic errors, so
    an *exact* match isn't expected (the same way plotting both submodes
    for an ordinary *finite* chain, see
    examples/dynamical_correlator/dynamical_correlator_time_evolution/
    main.py, only ever compares them visually, never asserts numerical
    agreement) -- but both should agree on *where* the dominant spectral
    weight sits. This mirrors that same example's own convention of
    comparing KPM's real part against the TD submode's own magnitude
    (`np.abs`), not its real part -- confirmed directly that using
    `.real` for the TD side here gives a spurious sign/scale mismatch
    that `np.abs` does not."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "idmrg"
    ic.maxm = 20
    ic.maxiter = 60
    ic.etol = 1e-12
    ic.niter = 30
    ic.set_hamiltonian(1.4 * ic.SzC[0] + ic.SxC[0] * ic.SxR[0])
    ic.gs_energy()

    es = np.linspace(-1, 6, 100)
    es_kpm, y_kpm = ic.kpm_finite("Sz", 0, "Sz", 1, n_window=12,
                                    window_chain_kwargs=dict(maxm=16, nsweeps=5),
                                    delta=0.3, es=es)

    from dmrgpy.pyitensor import idmrg_window
    from dmrgpy.timedependent import _fourier_transform_correlator
    ts, xs, S = idmrg_window.dynamical_correlator_td(
        ic._result, n_window=10, opname_A="Sz", opname_B="Sz", dt=0.05,
        nt=25, cutoff=1e-10, maxdim=30, niter=20, x_values=[1])
    es_td, g_td = _fourier_transform_correlator(ts, S[:, 0], 0.05, es=es,
                                                  delta=0.3, window=[-1, 6])

    peak_kpm = es[np.argmax(np.abs(y_kpm))]
    # Check that the TD spectrum still carries substantial weight *at*
    # the KPM peak location, rather than requiring it to be the TD
    # spectrum's own single dominant (argmax) feature -- confirmed
    # directly (after idmrg_window.py's own `eshift` fix, see
    # tests/test_idmrg_window_free_fermion.py's module docstring, which
    # this run now correctly reflects) that a second, comparable-height
    # low-frequency feature in the TD spectrum occasionally edges out the
    # one matching KPM's own peak by argmax alone (about 1 run in 3, not
    # obviously tied to iDMRG's own state_overlap) -- both features are
    # genuine, low-frequency, physically-plausible spectral weight for
    # this gapped model, so requiring "substantial weight at the KPM
    # peak" is the more robust way to express "these two independent
    # approximations agree on where the dominant spectral weight sits"
    # than a bare argmax-vs-argmax comparison.
    ix_kpm = np.argmin(np.abs(es_td - peak_kpm))
    weight_at_kpm_peak = np.abs(g_td[ix_kpm])
    # 0.5, not a much looser bound: confirmed directly over 4 independent
    # runs that this fraction lands at 0.899-1.0 whenever the KPM peak
    # isn't the TD spectrum's own argmax -- 0.5 stays well clear of that
    # observed range while being tight enough to catch a real regression
    # (e.g. a reintroduced partial phase error genuinely moving the
    # dominant weight away from the KPM peak, not just a near-degenerate
    # second feature edging out the argmax).
    assert weight_at_kpm_peak > 0.5 * np.max(np.abs(g_td))


# -- excitation_energies/excitation_gap: the tangent-space/quasiparticle
# excitation ansatz (pyitensor/idmrg_excitations.py) -- see that module's
# own docstring for the algorithm. Requires gs_method="vumps" (not the
# default "idmrg", see infinitechain.py's own class docstring): any
# converged bond dimension D>=1 is supported, including a genuinely
# entangled (D>1) ground state -- see idmrg_excitations.py's own "History"
# section for the investigation that used to make D>1 an explicit,
# rejected scope limit here.

def _polarized_xx_chain(J=1.0, h=3.0, maxm=4, maxiter=50, etol=1e-12, gs_method="idmrg"):
    """A field-polarized n_uc=1 XX chain (h > J, the XX chain's own
    saturation field) -- the ground state is the exact fully-polarized
    product state, so both iDMRG and VUMPS converge to a bond dimension
    D=1 unit cell. Exact single-magnon dispersion (free fermions via
    Jordan-Wigner): E(k) = 2*h - J*cos(k). `gs_method` defaults to
    "idmrg" (what local_excitation_gap's own tests below need);
    excitation_energies/excitation_gap's own tests pass "vumps" explicitly
    (see infinitechain.py's own class docstring for why the two features
    need different gs_method values)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = gs_method
    h_op = -J * (ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0]) - 2 * h * ic.SzC[0]
    ic.maxm, ic.maxiter, ic.etol = maxm, maxiter, etol
    ic.set_hamiltonian(h_op)
    ic.gs_energy()
    return ic


def test_excitation_energies_matches_exact_xx_dispersion():
    """The exact ground state here is a product state, so it is represented
    exactly at bond dimension D=1 -- which is what `maxm=1` asks for, and
    what makes this test deterministic.

    It used to run at the helper's default `maxm=4` and was flaky, ~33% of
    runs missing the 1e-8 tolerance; that was blamed on `vumps.py`'s
    unseeded `_random_raw_tensor`. Seeding is why the failures move around
    run to run, but it is not the cause. Measured over 40 independent runs
    at maxm=4, the ground-state energy was always excellent (|e0 error| <=
    1.8e-9) while the excitation error ranged over six orders of magnitude,
    up to 1.6e-3 -- an amplification of ~1e6. Sweeping the bond dimension
    isolates it cleanly (12 runs each): maxm=1 and 2 never failed (worst
    error 1.8e-15 and 1.8e-10), maxm=4 failed 2/12, maxm=8 failed 10/12,
    with e0 staying <= 8e-10 throughout. Padding an exactly-D=1 state out to
    larger D leaves near-null directions in the bond space, and the
    excitation ansatz's (1-E)-type inverses are ill-conditioned there.

    So this asks for the bond dimension the state actually needs, which is
    both deterministic (30/30 runs within 1.8e-15) and a fair test of what
    it claims to test -- the magnon dispersion formula, not VUMPS's
    behaviour when over-parameterized. That behaviour is a real limitation
    and is pinned separately by
    test_excitation_accuracy_degrades_with_redundant_bond_dimension."""
    J, h = 1.0, 3.0
    ic = _polarized_xx_chain(J, h, maxm=1, gs_method="vumps")
    for k in np.linspace(0, 2 * np.pi, 9, endpoint=False):
        exact = 2 * h - J * np.cos(k)
        got = ic.excitation_energies(k, n=1)[0]
        assert got == pytest.approx(exact, abs=1e-8)


def test_excitation_gap_matches_exact_minimum():
    """min_k (2h - J*cos(k)) = 2h - J, attained at k=0. `maxm=1` for the
    same reason as test_excitation_energies_matches_exact_xx_dispersion --
    see its docstring."""
    J, h = 1.0, 3.0
    ic = _polarized_xx_chain(J, h, maxm=1, gs_method="vumps")
    assert ic.excitation_gap() == pytest.approx(2 * h - J, abs=1e-8)


def test_excitation_accuracy_degrades_with_redundant_bond_dimension():
    """Pins the real cause of what used to be a ~33% flake in the two
    exact-dispersion tests above: representing an exactly-D=1 state at a
    larger bond dimension costs excitation accuracy, badly, while leaving
    the ground-state energy essentially untouched.

    The mechanism is conditioning, not convergence -- padding the state out
    leaves near-null directions in the bond space and the excitation
    ansatz's (1-E)-type inverses amplify them. Measured over 12 runs at each
    bond dimension: worst excitation error 1.8e-15 (maxm=1), 1.8e-10
    (maxm=2), 1.7e-7 (maxm=4), 3.5e-3 (maxm=8), with |e0 error| <= 8e-10
    throughout.

    Asserted loosely and one-sidedly, because the run-to-run spread at large
    maxm is itself large (that is the phenomenon): only that maxm=1 is
    essentially exact and that maxm=8 is at least a thousand times worse.
    Uses a fixed comparison within a single run rather than pinned numbers,
    so it cannot go stale against a future conditioning fix -- if someone
    fixes the conditioning, this test fails and should be updated to say
    so."""
    J, h = 1.0, 3.0

    def worst_error(maxm):
        ic = _polarized_xx_chain(J, h, maxm=maxm, gs_method="vumps")
        return max(abs(ic.excitation_energies(k, n=1)[0] - (2 * h - J * np.cos(k)))
                   for k in np.linspace(0, 2 * np.pi, 5, endpoint=False))

    exact_D = worst_error(1)
    padded = max(worst_error(8) for _ in range(3))  # spread is large; take the worst
    assert exact_D < 1e-12, exact_D
    assert padded > 1000 * max(exact_D, 1e-15), (exact_D, padded)


def test_excitation_environment_is_cached_across_calls():
    """A second excitation_energies call (different k) must reuse
    self._excitation_env (identity, not just an equal rebuild) -- it is
    expensive to build (a null-space computation plus dense linear
    solves)."""
    ic = _polarized_xx_chain(maxm=1, gs_method="vumps")
    ic.excitation_energies(0.0)
    env = ic._excitation_env
    assert env is not None
    ic.excitation_energies(1.0)
    assert ic._excitation_env is env


def test_excitation_gap_before_set_hamiltonian_raises():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "vumps"
    with pytest.raises(RuntimeError):
        ic.excitation_gap()


def test_excitation_energies_matches_exact_tfim_dispersion_d2():
    """A genuinely entangled (D=2) ground state -- the transverse-field
    Ising chain (J=1, g=2.5, gapped paramagnetic phase), same convention
    as test_vumps.py's own `_tfim_chain` (H = -4*SxC*SxR - 2*g*SzC =
    -sigma^x sigma^x - g*sigma^z) -- reproduces the exact free-fermion
    single-magnon dispersion eps(k) = 2*sqrt(J^2+g^2-2*J*g*cos(k)) to
    within VUMPS's own D=2 convergence, cross-checked independently
    against MPSKit.jl's own D=2 result to 6 significant figures (see
    idmrg_excitations.py's own "History" section and
    docs/idmrg_excitation_mpskit_port_plan.md)."""
    J, g = 1.0, 2.5
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "vumps"
    ic.maxm = 2
    ic.etol = 1e-12
    ic.vumps_nrestarts = 6
    ic.set_hamiltonian(-4.0 * J * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0])
    ic.gs_energy()
    assert ic.converged
    for k in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
        exact = 2 * np.sqrt(J ** 2 + g ** 2 - 2 * J * g * np.cos(k))
        got = ic.excitation_energies(k, n=1)[0]
        assert got.real == pytest.approx(exact, abs=2e-3)


def test_excitation_energies_gs_method_idmrg_not_implemented():
    """excitation_energies/excitation_gap need VUMPSResult's own mixed-
    gauge {AL,AR,C,GL,GR} -- the default gs_method="idmrg" (the growing
    algorithm) has no equivalent, so it is rejected explicitly rather than
    silently misused."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "idmrg"
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] - 2.0 * ic.SzC[0])
    with pytest.raises(NotImplementedError):
        ic.excitation_gap()


def test_excitation_energies_itensor_version3_not_implemented():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "idmrg"
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] - 2.0 * ic.SzC[0])
    with pytest.raises(NotImplementedError):
        ic.excitation_gap()


def test_excitation_energies_rejects_reach_greater_than_one():
    """A deliberately constructed longer-range term
    (get_operator(..., group="R") with i>=n_uc) spans 2 supersites after
    grouping -- rejected by idmrg_excitations._check_reach_one, called
    from within vumps.vumps_ground_state itself (i.e. by the implicit
    gs_energy() call inside excitation_gap(), not by the excitation
    machinery directly)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "vumps"
    far = ic.get_operator("Sx", 1, group="R")  # site n_uc+1 = 2, reach=2
    h = ic.SxC[0] * far - 2.0 * ic.SzC[0]
    ic.maxm, ic.maxiter, ic.etol = 4, 50, 1e-12
    ic.set_hamiltonian(h)
    with pytest.raises(NotImplementedError):
        ic.excitation_gap()


# -- local_excitation_gap: the "local superblock gap" (pyitensor/idmrg.py's
# local_excitation_gap) -- a cheap, cruder alternative to excitation_gap
# above: re-diagonalizes the growing algorithm's own final 2-site effective
# Hamiltonian for its second-lowest eigenvalue instead of only its ground
# state. Unlike excitation_gap, it has no momentum label and does not
# require D=1 -- but is correspondingly less accurate (confirmed below: an
# exact ~10% overestimate on the one case an exact answer is available for).

def test_local_excitation_gap_is_approximate_for_polarized_xx_chain():
    """D=1 sanity check: unlike excitation_gap (exact here, see
    test_excitation_gap_matches_exact_minimum), local_excitation_gap is a
    genuinely cruder approximation -- lands in the right ballpark (~10%
    high, measured directly) but should NOT be expected to match the exact
    answer to the same precision excitation_gap does."""
    J, h = 1.0, 3.0
    ic = _polarized_xx_chain(J, h)
    exact_gap = 2 * h - J
    local_gap = ic.local_excitation_gap()
    assert local_gap == pytest.approx(exact_gap, rel=0.2)
    assert abs(local_gap - exact_gap) > 1e-6  # genuinely approximate, not exact


def test_local_excitation_gap_is_deterministic_across_repeated_calls():
    """local_excitation_gap uses a randomized Lanczos start internally (no
    warm start persisted the way the growing algorithm's own solve has,
    since this is a one-off post-hoc diagonalization) -- confirm the
    deflated eigenproblem still converges to the same value every call,
    not call-to-call noise."""
    ic = _polarized_xx_chain()
    g1 = ic.local_excitation_gap()
    g2 = ic.local_excitation_gap()
    assert g1 == pytest.approx(g2, abs=1e-8)


def test_local_excitation_gap_matches_finite_size_dimerized_gap():
    """A genuinely entangled (D>1) ground state -- exactly the case
    excitation_gap cannot handle yet -- cross-checked against a large
    finite open chain's own ED gap (get_gap(mode="ED")), reusing the same
    dimerized model test_n_uc2_dimerized_chain_matches_finite_size_
    extrapolation already validates its own (energy-density) number
    against. Measured directly: local_excitation_gap=0.7631 vs. a finite
    open chain's own ED gap of 0.7591 (n=16 sites) -- well within the
    generous absolute tolerance below, which allows for the finite-size
    drift of the ED reference itself (0.7704 at n=12, 0.7591 at n=16)."""
    j_strong, j_weak = 1.0, 0.4
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    ic.gs_method = "idmrg"
    h = (j_strong * (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1])
         + j_weak * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0]))
    ic.maxm, ic.maxiter, ic.etol = 40, 200, 1e-9
    ic.set_hamiltonian(h)
    ic.gs_energy()
    assert ic.converged
    local_gap = ic.local_excitation_gap()

    sc = spinchain.Spin_Chain(["1/2"] * 16)
    hf = 0
    for i in range(15):
        j = j_strong if i % 2 == 0 else j_weak
        hf = hf + j * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                       + sc.Sz[i] * sc.Sz[i + 1])
    sc.set_hamiltonian(hf)
    finite_gap = sc.get_gap(mode="ED")

    assert local_gap == pytest.approx(finite_gap, abs=0.05)


def test_local_excitation_gap_before_set_hamiltonian_raises():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "idmrg"
    with pytest.raises(RuntimeError):
        ic.local_excitation_gap()


def test_local_excitation_gap_itensor_version3_requires_gs_method_idmrg():
    """window=0 IS implemented on the v3 C++ backend now
    (Chain::idmrg_local_excitation_gap, see
    test_idmrg_correlator_v3.py's own value check) -- but, exactly like
    the pyitensor backend, only under gs_method="idmrg": the growing
    algorithm's final local superblock is what gets re-diagonalized, and
    VUMPS never produces one."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "vumps"
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] - 2.0 * ic.SzC[0])
    with pytest.raises(NotImplementedError):
        ic.local_excitation_gap()


def test_local_excitation_gap_raises_without_stored_superblock():
    """Defense-in-depth guard on idmrg.local_excitation_gap itself (not
    reachable through the public API, which always builds a real
    IDMRGResult via a completed idmrg_ground_state run)."""
    result = idmrg.IDMRGResult(None, 1, [None], 0.0, True, 1)
    with pytest.raises(RuntimeError):
        idmrg.local_excitation_gap(result)


def _tfim_chain(J=1.0, h=2.0, maxm=12, maxiter=200, etol=1e-10):
    """A gapped (paramagnetic-phase, h>J), n_uc=1 transverse-field Ising
    chain -- unlike _polarized_xx_chain, the ground state is genuinely
    entangled (D>1), giving local_excitation_gap_windowed something
    nontrivial to work with. H = -4J*Sx_i*Sx_{i+1} - 2h*Sz_i (the 4/2
    factors convert dmrgpy's spin-1/2 S operators into Pauli matrices,
    sigma=2S, matching the textbook TFIM convention H = -J*sigmax*sigmax
    - h*sigmaz)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "idmrg"
    h_op = -4 * J * ic.SxC[0] * ic.SxR[0] - 2 * h * ic.SzC[0]
    ic.maxm, ic.maxiter, ic.etol = maxm, maxiter, etol
    ic.set_hamiltonian(h_op)
    ic.gs_energy()
    return ic


# -- local_excitation_gap(window=...) / idmrg.local_excitation_gap_windowed:
# widens the local diagonalization block with real, free physical sites on
# each side instead of just re-diagonalizing the frozen 2-site block --
# window=0 (the default) reduces to exactly local_excitation_gap above;
# window>0 measurably tightens the gap at the cost of an exponentially
# larger local Hilbert space (d**(2*window) more), and is only supported
# for n_uc=1 (see local_excitation_gap_windowed's own docstring for why
# n_uc=2 isn't handled yet).

def test_local_excitation_gap_windowed_matches_unwindowed_at_window_zero():
    """window=0 must reduce to exactly the same effective Hamiltonian
    local_excitation_gap diagonalizes -- an internal consistency check on
    the windowed construction itself (built independently of it, re-solving
    the ground state from scratch via Lanczos rather than reusing the
    growing algorithm's own evec0), not just a regression pin."""
    ic = _polarized_xx_chain()
    g_plain = idmrg.local_excitation_gap(ic._result)
    g_window0 = idmrg.local_excitation_gap_windowed(
        ic._result, ic._h_intra.op, ic._h_inter.op, ic.site_types, ic.n_uc,
        window=0)
    assert g_window0 == pytest.approx(g_plain, abs=1e-6)


def test_local_excitation_gap_windowed_improves_accuracy_for_polarized_xx_chain():
    """D=1 sanity check, exact answer known (see
    test_local_excitation_gap_is_approximate_for_polarized_xx_chain for the
    window=0 baseline, ~10% high): widening the window with real free sites
    converges monotonically toward the exact gap -- measured directly,
    window=0/1/2/4/6 give rel. errors of 0.10/0.038/0.020/0.0081/0.0044."""
    J, h = 1.0, 3.0
    ic = _polarized_xx_chain(J, h)
    exact_gap = 2 * h - J
    g0 = ic.local_excitation_gap(window=0)
    g2 = ic.local_excitation_gap(window=2)
    assert abs(g2 - exact_gap) < abs(g0 - exact_gap)
    assert g2 == pytest.approx(exact_gap, rel=0.03)


def test_local_excitation_gap_windowed_matches_finite_size_tfim_gap():
    """D>1 sanity check: a gapped, genuinely entangled transverse-field
    Ising chain (n_uc=1), cross-checked against a finite open chain's own
    ED gap -- measured directly: window=2 gap=2.058, window=3 gap=2.047,
    vs. an 18-site open chain's own ED gap of 2.049 -- converging at least
    as fast as growing the finite chain itself does (n=10/12/14/16/18 ED
    gaps are 2.133/2.099/2.076/2.060/2.049)."""
    ic = _tfim_chain()
    assert ic.converged
    local_gap = ic.local_excitation_gap(window=2)

    sc = spinchain.Spin_Chain(["1/2"] * 16)
    hf = 0
    for i in range(15):
        hf = hf + (-4.0) * sc.Sx[i] * sc.Sx[i + 1]
    for i in range(16):
        hf = hf + (-4.0) * sc.Sz[i]
    sc.set_hamiltonian(hf)
    finite_gap = sc.get_gap(mode="ED")
    assert local_gap == pytest.approx(finite_gap, abs=0.05)


def test_local_excitation_gap_windowed_deterministic_across_repeated_calls():
    ic = _polarized_xx_chain()
    g1 = ic.local_excitation_gap(window=1)
    g2 = ic.local_excitation_gap(window=1)
    assert g1 == pytest.approx(g2, abs=1e-6)


def test_local_excitation_gap_windowed_rejects_n_uc2():
    """local_excitation_gap_windowed needs to know which sublattice
    position each extra inserted site takes, which isn't tracked for
    n_uc=2 -- see its own docstring."""
    j_strong, j_weak = 1.0, 0.4
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    ic.gs_method = "idmrg"
    h = (j_strong * (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1])
         + j_weak * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0]))
    ic.maxm, ic.maxiter, ic.etol = 10, 50, 1e-9
    ic.set_hamiltonian(h)
    with pytest.raises(NotImplementedError):
        ic.local_excitation_gap(window=1)


def test_local_excitation_gap_windowed_before_set_hamiltonian_raises():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "idmrg"
    with pytest.raises(RuntimeError):
        ic.local_excitation_gap(window=1)


def test_local_excitation_gap_windowed_itensor_version3_not_implemented():
    """window>0 is deliberately NOT ported to the v3 C++ backend --
    idmrg.local_excitation_gap_windowed is an explicit prototype rather
    than stable API (see its own docstring), so it stays
    itensor_version="python"-only even though window=0 is now supported
    on both backends."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "idmrg"
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] - 2.0 * ic.SzC[0])
    with pytest.raises(NotImplementedError):
        ic.local_excitation_gap(window=1)


def test_local_excitation_gap_windowed_raises_without_stored_superblock():
    """Defense-in-depth guard on idmrg.local_excitation_gap_windowed
    itself (not reachable through the public API)."""
    result = idmrg.IDMRGResult(None, 1, [None], 0.0, True, 1)
    with pytest.raises(RuntimeError):
        idmrg.local_excitation_gap_windowed(result, [], [], ["1/2"], 1, window=1)


def test_local_excitation_gap_windowed_rejects_n_uc2_directly():
    """The n_uc guard fires before even looking at result.local_superblock
    -- checked directly at the module level with a bare-bones result."""
    result = idmrg.IDMRGResult(None, 2, [None, None], 0.0, True, 1)
    with pytest.raises(NotImplementedError):
        idmrg.local_excitation_gap_windowed(result, [], [], ["1/2", "1/2"], 2, window=1)


# -- Jordan-Wigner threading in idmrg._term_site_matrices --
# Regression tests for a real bug: the function omitted the extra trailing F
# that autompo.HTerm.resolve() applies to a touched site whose own fermion
# parity differs from the parity carried in from its left (`need_F`). The
# omission is invisible for a Cdag-leading term (Cdag @ F == Cdag) but flips
# the sign of the first site's matrix for a C-leading one (F @ C == -C), so
# an ordinary Hermitian hopping `Cdag[i]*C[j] + C[i]*Cdag[j]` silently came
# out non-Hermitian. Not reachable through Infinite_Spin_Chain (spin sites
# are never fermionic), but Infinite_Many_Body_Chain accepts fermionic site
# codes and get_operator() accepts any operator name, so this is a latent
# public-API bug, not dead code. HTerm.resolve() is the oracle here: it is
# this codebase's own, independently tested transcription of ITensor's
# autompo.cc rewriteFermionic().

_FERMION_SITE_CODE = 0  # spinless fermion site (C/Cdag/F/N)


@pytest.mark.parametrize("ops", [
    [["C", 0], ["Cdag", 1]],        # C-leading, nearest neighbour
    [["C", 0], ["Cdag", 2]],        # C-leading, one site skipped (real JW string)
    [["Cdag", 0], ["C", 1]],        # Cdag-leading (the case that always passed)
    [["Cdag", 0], ["C", 2]],
    [["N", 0], ["N", 1]],           # bosonic-parity pair: no F anywhere
])
def test_term_site_matrices_matches_autompo_resolve(ops):
    """_term_site_matrices' per-site matrices (stored (in,out) convention)
    must equal HTerm.resolve()'s own (standard convention, hence the .T)
    for every touched site, for both operator orderings."""
    from dmrgpy.pyitensor.autompo import HTerm
    from dmrgpy.pyitensor.sites import SiteX

    n_uc = 2
    sites_uc = SiteX([_FERMION_SITE_CODE] * n_uc)
    sites_ref = SiteX([_FERMION_SITE_CODE] * 4)

    rel_sites, _coef, mats, _ferm = idmrg._term_site_matrices(
        [1.0] + ops, sites_uc, n_uc)

    ht = HTerm(1.0)
    for name, site in ops:
        ht.add(name, site + 1)  # resolve() is 1-based
    reference = ht.resolve(sites_ref)

    for k, site in enumerate(rel_sites):
        assert np.allclose(mats[k].T, reference[site]), (
            "site {} disagrees with HTerm.resolve()".format(site))


def test_term_site_matrices_c_leading_picks_up_the_extra_F():
    """The specific regression: F @ C == -C, so the first site's matrix for a
    C-leading term is the *negative* of the naive composition. Checked
    against the explicit F-multiplied form rather than only against
    resolve(), so this test still fails loudly if resolve() itself were ever
    changed in the same wrong direction."""
    from dmrgpy.pyitensor.sites import SiteX

    sites_uc = SiteX([_FERMION_SITE_CODE, _FERMION_SITE_CODE])
    st = sites_uc.site_type(1)
    _rel, _coef, mats, _ferm = idmrg._term_site_matrices(
        [1.0, ["C", 0], ["Cdag", 1]], sites_uc, 2)
    assert np.allclose(mats[0], st.matrix("F") @ st.matrix("C"))
    assert np.allclose(mats[0], -st.matrix("C"))
    # the closing site's parity matches the carry, so it gets no extra F
    assert np.allclose(mats[1], st.matrix("Cdag"))


def test_term_site_matrices_rejects_odd_total_fermion_parity():
    """A term with odd total fermion parity (here a bare single C) opens a
    Jordan-Wigner string that is never closed -- it would have to run to
    infinity, which _build_periodic_mpo's finite per-term pending channels
    cannot represent (they terminate every string at the term's own last
    touched site). Such a term is not parity-conserving and cannot appear in
    a physical Hamiltonian; it must raise rather than silently build an MPO
    with a truncated string."""
    from dmrgpy.pyitensor.sites import SiteX

    sites_uc = SiteX([_FERMION_SITE_CODE, _FERMION_SITE_CODE])
    with pytest.raises(ValueError):
        idmrg._term_site_matrices([1.0, ["C", 0]], sites_uc, 2)
    with pytest.raises(ValueError):
        idmrg._term_site_matrices([1.0, ["Cdag", 0], ["N", 1]], sites_uc, 2)

def test_term_site_matrices_applies_the_cross_site_reordering_sign():
    """`_term_site_matrices` sorts a term's factors by site itself, so it
    must supply the fermionic anticommutation sign that sorting implies --
    `autompo.HTerm.resolve()` never has to, because `HTerm.add()` already
    applied it while insertion-sorting.

    Without it `Cdag_1 C_0` and `C_0 Cdag_1` come out identical when they
    differ by exactly -1. Pinned against `HTerm.add`'s own coefficient, which
    is the authoritative implementation of the rule, over every 2-factor
    fermionic ordering plus random even-parity terms."""
    import random
    from dmrgpy.pyitensor.autompo import HTerm
    from dmrgpy.pyitensor.sites import SiteX

    sites_uc = SiteX([_FERMION_SITE_CODE, _FERMION_SITE_CODE])

    def coef_of(ops):
        _rel, coef, _mats, _ferm = idmrg._term_site_matrices([1.0] + ops, sites_uc, 2)
        ht = HTerm(1.0)
        for nm, st in ops:
            ht.add(nm, st + 1)
        return coef, ht.coef

    # the specific pair that motivated this
    a, b = coef_of([["Cdag", 1], ["C", 0]])
    assert a == b == -1
    a, b = coef_of([["C", 0], ["Cdag", 1]])
    assert a == b == 1

    random.seed(11)
    for _ in range(200):
        ops = []
        for _ in range(random.choice([1, 2])):
            ops.append([random.choice(["C", "Cdag"]), random.randint(0, 1)])
            ops.append([random.choice(["C", "Cdag"]), random.randint(0, 1)])
        if random.random() < 0.3:
            ops.insert(random.randint(0, len(ops)), ["N", random.randint(0, 1)])
        got, expected = coef_of(ops)
        assert abs(got - expected) < 1e-12, (ops, got, expected)


# ---------------------------------------------------------------------------
# End-to-end fermionic infinite chains (both backends).
#
# Everything above this point tests fermionic *term handling* in isolation
# (_term_site_matrices vs HTerm.resolve()); every test that actually runs
# gs_energy() is a spin chain. That gap is why itensor_version=3 shipped
# unable to run a fermionic infinite chain at all: infinitechain.py used to
# hand the C++ backend MultiOperator.to_terms(), i.e. the *finite*-chain
# Jordan-Wigner form whose strings start at site 1 of the chain, which both
# violates the documented "at most 2 distinct sites per term" contract (each
# explicit F factor lands on its own site) and hardcodes a string anchored
# at a site an infinite chain does not have. It happened to work whenever
# every fermionic term connected *adjacent* sites -- the global string then
# collapses to exactly the right endpoint-F composition with nothing left
# over -- which is the only case a spin-only test suite could ever have
# exercised. Both backends now take the raw, untransformed terms and thread
# the string locally themselves (pyitensor/idmrg.py's _term_site_matrices,
# mpscpp3/chain_session.h's idmrg_classify_terms + idmrg_build_row).
#
# The oracle is a free-fermion band integral: for a quadratic Hamiltonian
# the ground state fills every negative-energy single-particle level, so
# the exact energy density is an integral over the occupied bands of the
# Bloch matrix -- backend-independent, and exact to machine precision
# rather than a golden value pinned to some past commit.

_SPINFUL_SITE_CODE = 1  # native spinful (Electron/Hubbard) site, dim 4


def _free_fermion_energy_density(t1, t2, t3, nk=20001):
    """Exact ground-state energy per SITE of the periodic spinless chain
    with hoppings t1 (A_n-B_n), t2 (B_n-A_{n+1}) and t3 (A_n-A_{n+1}),
    two sites (A,B) per unit cell:

        H(k) = [[2 t3 cos k, t1 + t2 e^{-ik}], [t1 + t2 e^{ik}, 0]]

    t2 != t1 dimerizes the chain (opening a gap), and t3 is what makes this
    a genuine test of the Jordan-Wigner string: it couples the two A sites
    of adjacent cells, i.e. a term whose two endpoints have a site strictly
    between them, so its bond carries carry_ferm=True."""
    k = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    H = np.zeros((len(k), 2, 2), dtype=complex)
    H[:, 0, 0] = 2 * t3 * np.cos(k)
    H[:, 0, 1] = t1 + t2 * np.exp(-1j * k)
    H[:, 1, 0] = np.conj(H[:, 0, 1])
    w = np.linalg.eigvalsh(H)
    return np.where(w < 0, w, 0.0).sum(axis=1).mean() / 2


def _free_fermion_chain(t1, t2, t3, itensor_version, gs_method, maxm=8,
                         spinful=False):
    site = _SPINFUL_SITE_CODE if spinful else _FERMION_SITE_CODE
    ic = infinitechain.Infinite_Many_Body_Chain(
        [site, site], itensor_version=itensor_version)
    ic.gs_method = gs_method
    ic.maxm = maxm
    ic.maxiter = 40
    H = 0
    # A native spinful site carries both flavours, so the same hoppings are
    # written once per flavour and the two decouple exactly.
    for suffix in (["up", "dn"] if spinful else [""]):
        C = [ic.get_operator("C" + suffix, i, "C") for i in range(2)]
        Cd = [ic.get_operator("Cdag" + suffix, i, "C") for i in range(2)]
        CR = [ic.get_operator("C" + suffix, i, "R") for i in range(2)]
        CdR = [ic.get_operator("Cdag" + suffix, i, "R") for i in range(2)]
        H = H + t1 * (Cd[0] * C[1] + Cd[1] * C[0])
        H = H + t2 * (Cd[1] * CR[0] + CdR[0] * C[1])
        H = H + t3 * (Cd[0] * CR[0] + CdR[0] * C[0])
    ic.set_hamiltonian(H)
    return ic


@pytest.mark.parametrize("itensor_version,gs_method", [
    ("python", "idmrg"),
    ("python", "vumps"),
    pytest.param(3, "idmrg", marks=pytest.mark.skipif(
        not cppext.available(3), reason="mpscpp3 extension not compiled")),
    pytest.param(3, "vumps", marks=pytest.mark.skipif(
        not cppext.available(3), reason="mpscpp3 extension not compiled")),
])
def test_free_fermion_energy_density_matches_band_integral(itensor_version,
                                                            gs_method):
    """A gapped, dimerized spinless chain including a next-cell A-A hopping
    (the term that needs a Jordan-Wigner string on the site between its two
    endpoints) must reproduce the exact band-integral energy density on
    every backend. Without the string this lands at a visibly different
    number, not a slightly less converged one; before the fix
    itensor_version=3 aborted outright on this Hamiltonian."""
    t1, t2, t3 = 1.0, 0.4, 0.1
    exact = _free_fermion_energy_density(t1, t2, t3)
    ic = _free_fermion_chain(t1, t2, t3, itensor_version, gs_method)
    assert ic.gs_energy() == pytest.approx(exact, abs=1e-6)


@pytest.mark.parametrize("itensor_version", [
    "python",
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3), reason="mpscpp3 extension not compiled")),
])
def test_native_spinful_free_chain_is_twice_the_spinless_one(itensor_version):
    """The same free chain on native spinful (site_type=1) sites: two
    decoupled flavours, so the exact energy density is exactly 2x the
    spinless one. This is the only test of ElectronSite's own on-site spin
    convention (electron.h defines Cdn = Fup.Adn, an intra-site string no
    spinless site type has) -- do not fold it into the spinless test above.

    The tolerance is loose because the tolerance is not the point: a
    dim-4 site at maxm=8 is genuinely under-converged (each flavour gets
    the entanglement budget the spinless test spends on one), while a
    missing/incorrect Jordan-Wigner string moves this by O(1). The tight
    check is the backend-vs-backend one in
    test_native_spinful_backends_agree."""
    t1, t2, t3 = 1.0, 0.4, 0.1
    exact = 2 * _free_fermion_energy_density(t1, t2, t3)
    ic = _free_fermion_chain(t1, t2, t3, itensor_version, "idmrg", spinful=True)
    assert ic.gs_energy() == pytest.approx(exact, abs=2e-3)


@pytest.mark.skipif(not cppext.available(3),
                     reason="mpscpp3 extension not compiled")
def test_native_spinful_backends_agree():
    """itensor_version=3 and itensor_version="python" must agree on a
    native-spinful chain far more tightly than either agrees with the
    under-converged exact value -- they run the same algorithm at the same
    bond dimension, so any disagreement here is a porting bug in one of the
    two Jordan-Wigner implementations, not a convergence difference."""
    t1, t2, t3 = 1.0, 0.4, 0.1
    e_py = _free_fermion_chain(t1, t2, t3, "python", "idmrg", spinful=True).gs_energy()
    e_v3 = _free_fermion_chain(t1, t2, t3, 3, "idmrg", spinful=True).gs_energy()
    assert e_v3 == pytest.approx(e_py, abs=1e-6)


@pytest.mark.skipif(not cppext.available(3),
                     reason="mpscpp3 extension not compiled")
def test_interacting_spinful_cell_backends_agree():
    """A genuinely interacting native-spinful unit cell -- one itinerant
    ("c") and one correlated ("f") orbital, with a Hubbard U on f, a Kondo
    exchange, and a hybridization written as a three-operator product
    confined to two sites. That last shape is what originally motivated
    this: it is a perfectly legal 2-site term, but the finite-chain
    Jordan-Wigner form of its inter-cell partners spilled onto extra sites
    and the C++ classifier rejected the whole Hamiltonian.

    What this test does NOT do is check the two backends' converged energies
    against each other tightly, and the reason is measured rather than
    assumed. This Kondo-lattice-like cell converges slowly and both backends
    start from unseeded random states, so repeated runs at maxm=8 scatter by
    as much as 5e-2 -- while deliberately dropping the Jordan-Wigner string
    (forcing carry_ferm=False, i.e. reintroducing exactly the bug this work
    fixed) moves the energy by only 2.9e-2. The noise is larger than the
    signal, so a cross-backend energy tolerance here would either be too
    tight to pass reliably or too loose to catch the bug it appears to guard.

    The numerical guard lives in the free-fermion tests above instead, where
    an exact band integral pins every backend to ~1e-6 and the same string
    machinery is exercised. What is left for this test is the thing those
    cannot cover and the thing that actually broke: that a three-operator
    product confined to two sites, with interacting partners across the cell
    boundary, is *accepted and runs* on both backends. Before the fix it
    aborted the v3 run outright, so "both produce a finite energy in the same
    ballpark" is a real regression check, and it is stated at the strength it
    actually has."""
    U1, tc, J = 0.5, 0.2, 0.2

    def build(itensor_version):
        ic = infinitechain.Infinite_Many_Body_Chain(
            [_SPINFUL_SITE_CODE] * 2, itensor_version=itensor_version)
        ic.gs_method = "idmrg"
        ic.maxm = 8
        ic.maxiter = 30
        g = ic.get_operator
        Cup = [g("Cup", i, "C") for i in range(2)]
        Cdup = [g("Cdagup", i, "C") for i in range(2)]
        Cdn = [g("Cdn", i, "C") for i in range(2)]
        Cddn = [g("Cdagdn", i, "C") for i in range(2)]
        Nup = [g("Nup", i, "C") for i in range(2)]
        Ndn = [g("Ndn", i, "C") for i in range(2)]
        Sx = [0.5 * Cdup[i] * Cdn[i] + 0.5 * Cddn[i] * Cup[i] for i in range(2)]
        Sy = [-0.5j * Cdup[i] * Cdn[i] + 0.5j * Cddn[i] * Cup[i] for i in range(2)]
        Sz = [0.5 * Nup[i] - 0.5 * Ndn[i] for i in range(2)]
        H = 0
        for suffix in ["up", "dn"]:
            Cd = [g("Cdag" + suffix, i, "C") for i in range(2)]
            C = [g("C" + suffix, i, "C") for i in range(2)]
            CdR = [g("Cdag" + suffix, i, "R") for i in range(2)]
            CR = [g("C" + suffix, i, "R") for i in range(2)]
            H = H + 1.0 * (Cd[0] * CR[0] + CdR[0] * C[0])   # c-orbital hopping
            H = H + 0.3 * (Cd[1] * CR[1] + CdR[1] * C[1])   # f-orbital hopping
        c, f = 0, 1
        Tup = -1j * (Cdup[f] * Cup[c] - Cdup[c] * Cup[f])
        H = H + U1 / 2 * (Nup[f] - 0.5) * (Ndn[f] - 0.5)
        H = H + tc / 2 * (Ndn[f] - 0.5) * Tup               # 3 factors, 2 sites
        H = H + J * (Sx[f] * Sx[c] + Sy[f] * Sy[c] + Sz[f] * Sz[c])
        ic.set_hamiltonian(H)
        return ic

    e_v3 = build(3).gs_energy()
    e_py = build("python").gs_energy()
    assert np.isfinite(e_v3) and np.isfinite(e_py)
    # 0.1 is above the ~5e-2 run-to-run scatter measured for this model and
    # far below the O(1) nonsense a genuinely mis-built Hamiltonian gives;
    # see the docstring for why it is deliberately not tighter.
    assert e_v3 == pytest.approx(e_py, abs=0.1)


def test_odd_fermion_parity_term_is_rejected_before_the_backend_sees_it():
    """A term with odd total fermion parity opens a Jordan-Wigner string
    that can never close on an infinite chain. Both backends reject it, but
    the v3 one does so with an ITensor Error() that aborts the process, so
    set_hamiltonian has to catch it first -- on the raw terms, before any
    backend is chosen."""
    ic = infinitechain.Infinite_Many_Body_Chain([_FERMION_SITE_CODE] * 2)
    C = ic.get_operator("C", 0, "C")
    with pytest.raises(ValueError, match="odd total fermion parity"):
        ic.set_hamiltonian(C + ic.get_operator("Cdag", 0, "C"))


def test_raw_and_jordan_wigner_terms_agree_for_spin_operators():
    """infinitechain.py now hands the v3 backend
    to_terms(jordan_wigner_transform=False). For a spin Hamiltonian the
    Jordan-Wigner transform is the identity, so the two forms must be
    literally the same term list -- which is what makes the whole existing
    spin-chain v3 test suite above a regression net for that plumbing
    change."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    Sx = [ic.get_operator("Sx", i, "C") for i in range(2)]
    SxR = [ic.get_operator("Sx", i, "R") for i in range(2)]
    Sz = [ic.get_operator("Sz", i, "C") for i in range(2)]
    ic.set_hamiltonian(Sx[0] * Sx[1] + Sx[1] * SxR[0] + 0.3 * Sz[0])
    for mo in (ic._h_intra, ic._h_inter):
        assert mo.to_terms() == mo.to_terms(jordan_wigner_transform=False)


# ---------------------------------------------------------------------------
# Dense vs iterative linear algebra.
#
# The transfer-matrix eigenproblem and the environment solves both used to
# be done densely on chi^2 x chi^2 matrices, which is O(chi^6) twice over
# and is what made VUMPS impractical at any useful bond dimension (a
# profile of one chi=12 run spent 14.7s of ~20s inside np.linalg.eig
# alone). Both now switch to a Krylov method above a size threshold, with
# the dense route kept as an exact fallback. The thresholds are set low
# enough to be crossed in practice, and these tests pin them to 0/infinity
# to check the two routes actually agree.


def _random_transfer_tensors(chi, d, n_uc, seed):
    from dmrgpy.pyitensor import idmrg_excitations as idmrg_exc

    rng = np.random.default_rng(seed)
    Es = []
    for _ in range(n_uc):
        A = rng.normal(size=(chi, d, chi)) + 1j * rng.normal(size=(chi, d, chi))
        Es.append(idmrg_exc._op_transfer_matrix(A, A))
    return Es


@pytest.mark.parametrize("n_uc", [1, 2])
def test_iterative_dominant_fixed_point_matches_dense(n_uc, monkeypatch):
    """The ARPACK route must return the same dominant eigenpair as the
    dense np.linalg.eig route, for both the left and the right problem."""
    Es = _random_transfer_tensors(5, 3, n_uc, seed=4)
    for fixed_point in (idmrg._dominant_right_fixed_point,
                         idmrg._dominant_left_fixed_point):
        monkeypatch.setattr(idmrg, "_DENSE_EIG_MAX", 10 ** 9)
        rho_dense, eta_dense = fixed_point(Es)
        monkeypatch.setattr(idmrg, "_DENSE_EIG_MAX", 0)
        rho_iter, eta_iter = fixed_point(Es)
        assert eta_iter == pytest.approx(eta_dense, rel=1e-10)
        assert np.allclose(rho_iter, rho_dense, atol=1e-8)


def test_excitation_iterative_eigensolver_matches_dense(monkeypatch):
    """`excitation_energies`' Lanczos route must reproduce the dense
    assemble-and-`eigh` route it replaces above `_DENSE_EIG_MAX`.

    The point of the iterative route is that it needs only tens of
    `_h_eff_action` applications instead of `dim = D*D*(d_g-1)` of them, so
    the two paths are only ever compared here -- in production one or the
    other runs, picked by size. Uses a D=2 TFIM (dim=4, normally dense) with
    the threshold monkeypatched both ways, the same style as
    test_iterative_linear_solve_matches_dense below."""
    from dmrgpy.pyitensor import idmrg_excitations as idmrg_exc

    J, g = 1.0, 2.5
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "vumps"
    ic.maxm = 2
    ic.etol = 1e-12
    ic.vumps_nrestarts = 6
    ic.set_hamiltonian(-4.0 * J * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0])
    ic.gs_energy()
    env = ic._get_excitation_environment()

    for k in [0.0, 0.7, 2.6]:
        monkeypatch.setattr(idmrg_exc, "_DENSE_EIG_MAX", 10 ** 9)
        dense = idmrg_exc.excitation_energies(env, k, n=2)
        monkeypatch.setattr(idmrg_exc, "_DENSE_EIG_MAX", 0)
        iterative = idmrg_exc.excitation_energies(env, k, n=2)
        assert np.allclose(iterative, dense, atol=1e-9), (k, dense, iterative)


def test_excitation_dense_path_survives_dim_one():
    """A D=1 spin-1/2 chain has dim = D*D*(d_g-1) = 1, where ARPACK cannot
    be used at all (it needs nev < dim). The dense path is what keeps this
    working, so it is not an optimization that could be dropped once the
    iterative one exists -- pinned here because the threshold that normally
    protects it is a tunable constant."""
    from dmrgpy.pyitensor import idmrg_excitations as idmrg_exc

    J, h = 1.0, 3.0
    ic = _polarized_xx_chain(J, h, maxm=1, gs_method="vumps")
    env = ic._get_excitation_environment()
    assert env.D * env.D * (env.d_g - 1) == 1
    got = idmrg_exc.excitation_energies(env, 0.7, n=1)
    assert got[0] == pytest.approx(2 * h - J * np.cos(0.7), abs=1e-8)


def test_excitation_energies_can_return_eigenvectors():
    """`return_vectors=True` hands back the tangent-space parameters X
    alongside the energies -- what a spectral weight is built from. Checks
    the shape contract and that each X really is an eigenvector of the same
    H_eff(k) the energy came from."""
    from dmrgpy.pyitensor import idmrg_excitations as idmrg_exc

    J, g = 1.0, 2.5
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "vumps"
    ic.maxm = 3
    ic.etol = 1e-12
    ic.vumps_nrestarts = 6
    ic.set_hamiltonian(-4.0 * J * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0])
    ic.gs_energy()
    env = ic._get_excitation_environment()

    k = 0.9
    w, xs = idmrg_exc.excitation_energies(env, k, n=2, return_vectors=True)
    assert len(xs) == len(w) == 2
    Dx = env.D * (env.d_g - 1)
    for X, lam in zip(xs, w):
        assert X.shape == (Dx, env.D)
        # H_eff(k)[X] = (lam + lam_AC) * X -- `excitation_energies` reports
        # the energy *above* the ground state, i.e. with lam_AC subtracted.
        residual = idmrg_exc._h_eff_action(k, X, env) - (lam + env.lam_AC) * X
        assert np.linalg.norm(residual) < 1e-8, np.linalg.norm(residual)


def test_excitation_resolvents_are_cached_per_momentum():
    """The GBL/GBR channel resolvents depend only on the momentum and the
    environment, never on the excitation tensor B -- so one per momentum is
    enough, and rebuilding them inside every `_h_eff_action` call (which is
    what used to happen) repeated a dense build plus an O(D^6)
    factorization for nothing."""
    from dmrgpy.pyitensor import idmrg_excitations as idmrg_exc

    ic = _polarized_xx_chain(maxm=1, gs_method="vumps")
    env = ic._get_excitation_environment()
    assert env.resolvent_cache == {}

    built = []
    real_resolvent = idmrg_exc._channel_resolvent

    def counting_resolvent(*args, **kwargs):
        built.append(1)
        return real_resolvent(*args, **kwargs)

    idmrg_exc._channel_resolvent = counting_resolvent
    try:
        idmrg_exc.excitation_energies(env, 0.4, n=1)
        after_first = len(built)
        idmrg_exc.excitation_energies(env, 0.4, n=1)
        assert len(built) == after_first  # same k: fully cached
        idmrg_exc.excitation_energies(env, 1.1, n=1)
        assert len(built) == after_first + 2  # a new k: one per direction
    finally:
        idmrg_exc._channel_resolvent = real_resolvent
    assert set(env.resolvent_cache) == {0.4, 1.1}


def test_iterative_linear_solve_matches_dense(monkeypatch):
    """_solve_linear_map's GMRES route must reproduce the dense LU route."""
    from dmrgpy.pyitensor import idmrg_excitations as idmrg_exc

    D = 5
    rng = np.random.default_rng(9)
    M = (rng.normal(size=(D * D, D * D)) + 1j * rng.normal(size=(D * D, D * D))
         + 3 * np.eye(D * D))
    rhs = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))

    def action(X):
        return (M @ X.reshape(-1)).reshape(D, D)

    monkeypatch.setattr(idmrg_exc, "_DENSE_SOLVE_MAX", 10 ** 9)
    dense = idmrg_exc._solve_linear_map(D, action, rhs)
    monkeypatch.setattr(idmrg_exc, "_DENSE_SOLVE_MAX", 0)
    iterative = idmrg_exc._solve_linear_map(D, action, rhs)
    assert np.allclose(iterative, dense, atol=1e-10)


def test_d_ramp_doubles_and_ends_at_the_target():
    """The VUMPS driver solves at every bond dimension the ramp lists, so
    its length is a direct multiplier on the whole run's cost."""
    from dmrgpy.pyitensor.vumps import _d_ramp

    for D in range(1, 65):
        ramp = _d_ramp(D)
        assert ramp[-1] == D
        assert ramp[0] == 1
        assert all(b > a for a, b in zip(ramp, ramp[1:]))
        assert len(ramp) <= D
    assert _d_ramp(30) == [1, 2, 4, 8, 16, 30]


# ---------------------------------------------------------------------------
# Fermionic two-point correlators (the Jordan-Wigner string on the sites
# BETWEEN the two operators).
#
# These used to be refused outright: every correlator path inserted the bare
# operator into an otherwise-identity transfer matrix, which is the right
# thing for a spin/boson operator and silently wrong for a fermionic one --
# a stringless <Cdag_i C_j> is not merely imprecise, it is a different
# quantity, and one whose error grows with separation (so it corrupts a
# decay rate, which is usually what such a correlator is measured for).
#
# All three paths now thread the string, and all three take their endpoint
# matrices and their "is a string open at all" flag from the SAME helper
# that builds the Hamiltonian's own 2-site terms (idmrg._term_site_matrices
# / the C++ vumps_correlator_endpoints port of it), so a correlator and a
# Hamiltonian term written with the same operator names cannot drift apart.
#
# Oracle: a free-fermion one-body density matrix, exact by construction.
#
# Coverage note: every test here uses spinless (site_type=0) sites, because
# a native spinful (site_type=1, dim 4) chain is too slow to belong in a
# routine suite -- one v3/VUMPS correlator check at maxm=8 takes ~255s. It
# was verified by hand instead, on the same dimerized model, and the
# ElectronSite-specific part of the convention (electron.h defines
# Cdn = Fup.Adn, an intra-site string no spinless site type has) is
# exercised there: <Cdagup_0 Cup_r> for r=0..3 came out
# [0.511161, -0.478454, -0.029260, 0.089760] against an exact
# [0.510503, -0.477696, -0.027524, 0.092718], i.e. max error 3.0e-03 at
# maxm=8 and 4.5e-04 at maxm=10 -- ordinary under-convergence for dim-4
# sites at that bond dimension, with the sign structure exact, whereas a
# mis-threaded string moves these by O(0.1). The Hamiltonian side of the
# same ElectronSite convention IS covered automatically, by
# test_native_spinful_free_chain_is_twice_the_spinless_one above.


def _exact_one_body_density_matrix(t1, t2, t3, L=200):
    """P[i,j] = <c^dag_i c_j> for the periodic dimerized spinless ring of L
    cells (2 sites per cell, site index 2*cell+sublattice), all negative
    single-particle levels filled. Gapped, so the bulk of a large ring is
    the infinite-chain answer to far better than the tolerances below."""
    N = 2 * L
    H = np.zeros((N, N))
    for n in range(L):
        a, b, a2 = 2 * n, 2 * n + 1, (2 * n + 2) % N
        H[a, b] += t1; H[b, a] += t1
        H[b, a2] += t2; H[a2, b] += t2
        H[a, a2] += t3; H[a2, a] += t3
    w, v = np.linalg.eigh(H)
    occ = v[:, w < 0]
    return (occ.conj() @ occ.T).real


@pytest.mark.parametrize("itensor_version,gs_method,maxm,tol", [
    ("python", "idmrg", 20, 1e-6),
    ("python", "vumps", 12, 1e-5),
    pytest.param(3, "vumps", 12, 1e-5, marks=pytest.mark.skipif(
        not cppext.available(3), reason="mpscpp3 extension not compiled")),
])
def test_fermionic_correlator_matches_free_fermion_exact(itensor_version,
                                                          gs_method, maxm, tol):
    """<Cdag_0 C_r> across a range of separations, against the exact
    free-fermion one-body density matrix. r=2 and beyond have at least one
    site strictly between the two operators, so they are the ones that need
    the string; without it these come out as different numbers entirely
    (dropping the string was measured to move the r=0..8 values by up to
    ~0.5), not as slightly-less-converged ones."""
    t1, t2, t3 = 1.0, 0.4, 0.1
    P = _exact_one_body_density_matrix(t1, t2, t3)
    i0 = 2 * (P.shape[0] // 4)          # a bulk A site
    ic = _free_fermion_chain(t1, t2, t3, itensor_version, gs_method, maxm=maxm)
    for r in range(7):
        assert complex(ic.correlator("Cdag", 0, "C", r)).real == \
            pytest.approx(P[i0, i0 + r], abs=tol), "separation r={}".format(r)


@pytest.mark.parametrize("itensor_version,gs_method", [
    ("python", "idmrg"),
    ("python", "vumps"),
    pytest.param(3, "vumps", marks=pytest.mark.skipif(
        not cppext.available(3), reason="mpscpp3 extension not compiled")),
])
def test_fermionic_correlator_n_uc_1_crosses_whole_cells(itensor_version,
                                                          gs_method):
    """With n_uc=1 every separation r>1 puts one or more COMPLETE unit cells
    between the two operators, which is a distinct code path from a string
    that stays inside one cell (in the VUMPS backends it is the branch that
    applies F to every sub-site of each fully-crossed cell's transfer
    tensor). The uniform half-filled chain has <Cdag_0 C_r> = 0 for every
    even r>0 by particle-hole symmetry -- a sharp, tolerance-free signature
    that survives even though this model is gapless and so converges only
    slowly in maxm. Without the string those entries are NOT small."""
    ic = infinitechain.Infinite_Many_Body_Chain(
        [_FERMION_SITE_CODE], itensor_version=itensor_version)
    ic.gs_method = gs_method
    ic.maxm = 16
    ic.maxiter = 60
    # This chain is gapless, which is VUMPS's documented weak spot (see
    # pyitensor/vumps.py's "Convergence robustness" section). The extra
    # restarts are about reaching a converged state to measure, not about
    # the correlator being delicate.
    #
    # They used to be covering for something else as well. Being
    # half-filled, this chain is critical with 2k_F = pi, so its transfer
    # matrix carries a peripheral eigenvalue at phase pi whose magnitude
    # approaches the leading one -- and idmrg's dominant-fixed-point guard
    # compared MAGNITUDES, so it rejected that perfectly well-posed state
    # as "(near-)degenerate". Measured at D=16: ~42% of individual VUMPS
    # attempts and ~16% of whole solves died there, i.e. this test failed
    # intermittently even at nrestarts=10. Fixed in
    # idmrg._check_dominant_eigenvalue_nondegenerate (a repeated
    # EIGENVALUE is degenerate; +rho and -rho are not) and pinned by
    # tests/test_transfer_degeneracy_guard.py, so the restarts below are
    # back to covering only what they say they cover.
    ic.vumps_nrestarts = 10
    Cd = ic.get_operator("Cdag", 0, "C")
    C = ic.get_operator("C", 0, "C")
    CdR = ic.get_operator("Cdag", 0, "R")
    CR = ic.get_operator("C", 0, "R")
    ic.set_hamiltonian(Cd * CR + CdR * C)
    ic.gs_energy()
    for r in (2, 4):
        assert abs(complex(ic.correlator("Cdag", 0, "C", r)).real) < 5e-3
    # ...while the odd separations are large and alternate in sign
    assert complex(ic.correlator("Cdag", 0, "C", 1)).real < -0.2
    assert complex(ic.correlator("Cdag", 0, "C", 3)).real > 0.05


def test_stringless_fermionic_correlator_would_be_a_different_number():
    """Guards the guard: confirms the string is actually load-bearing here,
    so the tests above would fail if it were dropped rather than passing for
    some unrelated reason. Recomputes r=2..6 with the string suppressed and
    requires the result to differ far beyond the tolerances used above."""
    from dmrgpy.pyitensor import idmrg as _idmrg

    t1, t2, t3 = 1.0, 0.4, 0.1
    ic = _free_fermion_chain(t1, t2, t3, "python", "idmrg", maxm=20)
    ic.gs_energy()
    strung = [complex(ic.correlator("Cdag", 0, "C", r)).real for r in range(2, 7)]

    original = _idmrg._term_site_matrices

    def no_string(op_term, sites_uc, n_uc):
        rel, coef, mats, ferm = original(op_term, sites_uc, n_uc)
        return rel, coef, mats, [False] * len(ferm)

    _idmrg._term_site_matrices = no_string
    try:
        bare = [complex(ic.correlator("Cdag", 0, "C", r)).real for r in range(2, 7)]
    finally:
        _idmrg._term_site_matrices = original

    assert max(abs(a - b) for a, b in zip(strung, bare)) > 1e-3


def test_odd_total_parity_correlator_is_rejected():
    """A pair whose total fermion parity is odd (one fermionic operator and
    one parity-even one) opens a string that never closes. It vanishes
    identically in any parity-conserving state, and the v3 backend rejects
    it with an ITensor Error() that would abort the process, so it is caught
    in infinitechain.py before any backend sees it."""
    ic = _free_fermion_chain(1.0, 0.4, 0.1, "python", "idmrg", maxm=6)
    with pytest.raises(ValueError, match="odd total fermion parity"):
        ic.correlator("Cdag", 0, "N", 1)
