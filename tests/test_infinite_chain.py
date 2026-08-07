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
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0])
    ic.maxm = 30
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
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + 0.5 * (ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]))
    ic.maxm = 30
    ic.maxiter = 150
    ic.etol = 1e-8
    ic.set_hamiltonian(h)
    e0 = ic.gs_energy()
    assert np.isfinite(e0)
    assert e0 < 0


def test_n_uc_above_2_not_implemented():
    with pytest.raises(NotImplementedError):
        infinitechain.Infinite_Spin_Chain(["1/2", "1/2", "1/2"])


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
    mpscpp3/chain_session.h's Chain::idmrg_ground_state is a line-by-line
    port of pyitensor/idmrg.py's own idmrg_ground_state, so the two should
    land on the same truncated-bond-dimension energy density, not just
    agree in the maxm->infinity limit. Confirmed directly to agree to
    ~1e-8 (well inside the 1e-6 tolerance here) and to be stable run over
    run to ~1e-8 despite iDMRG's unseeded random starting MPS."""
    ic_py, _ = _converged_uniform_chain(2, maxm=30, maxiter=60, etol=1e-9,
                                         itensor_version="python")
    ic_v3, _ = _converged_uniform_chain(2, maxm=30, maxiter=60, etol=1e-9,
                                         itensor_version=3)
    assert ic_v3.e0 == pytest.approx(ic_py.e0, abs=1e-6)


@pytest.mark.skipif(not cppext.available(3),
                     reason="requires the compiled mpscpp3 (ITensor v3) extension")
def test_itensor_version3_vev_and_correlator_not_implemented():
    """The v3 C++ backend has no static-correlator machinery yet (see
    Infinite_Many_Body_Chain.gs_energy's own comment) -- vev/correlator
    must raise a clear NotImplementedError rather than silently misusing
    a stale/absent self._result."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    h = ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]
    ic.set_hamiltonian(h)
    with pytest.raises(NotImplementedError):
        ic.vev("Sz", 0)
    with pytest.raises(NotImplementedError):
        ic.correlator("Sz", 0, "Sz", 1)


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
        idmrg.grow_by_mpo(W[:1], ic._result.U_list, 2)


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
    ic_up.maxm, ic_up.maxiter, ic_up.etol = 4, 50, 1e-12
    ic_up.set_hamiltonian(-10.0 * ic_up.SzC[0])
    ic_up.gs_energy()

    ic_down = infinitechain.Infinite_Spin_Chain(["1/2"])
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
    ic_half.maxm, ic_half.maxiter, ic_half.etol = 4, 20, 1e-8
    ic_half.set_hamiltonian(ic_half.SxC[0] * ic_half.SxR[0]
                             + ic_half.SyC[0] * ic_half.SyR[0]
                             + ic_half.SzC[0] * ic_half.SzR[0])
    ic_half.gs_energy()

    ic_one = infinitechain.Infinite_Spin_Chain(["1"])
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
    ic.maxm, ic.maxiter, ic.etol = 4, 50, 1e-12
    ic.set_hamiltonian(-5.0 * ic.SyC[0])
    ic.gs_energy()
    result = ic._result

    op_term = [1.0, ["Sx", 0], ["Sz", 0]]
    _, _, mats, _ = idmrg._term_site_matrices(op_term, result.sites_uc, 1)
    M = mats[0]
    Es = idmrg._transfer_matrices(result.U_list, 1)
    rho_after, _ = idmrg._all_right_fixed_points(Es, 1)
    A = idmrg._to_array_lpr(result.U_list[0])
    Aop = np.einsum('io,lir->lor', M, A)
    E4 = np.einsum('lpr,LpR->lLrR', Aop, np.conj(A))
    expected = complex(np.trace(idmrg._apply_transfer(E4, rho_after[0])))

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
    ic_up.maxm, ic_up.maxiter, ic_up.etol = 4, 50, 1e-12
    ic_up.set_hamiltonian(-10.0 * ic_up.SzC[0])
    ic_up.gs_energy()

    ic_down = infinitechain.Infinite_Spin_Chain(["1/2"])
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
    ic_dom, _ = _converged_uniform_chain(n_uc, maxm=20, maxiter=100, etol=1e-11)

    spins = ["1/2"] * n_uc
    ic_other = infinitechain.Infinite_Spin_Chain(spins)
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
    ic_other.maxm, ic_other.maxiter, ic_other.etol = 20, 100, 1e-11
    ic_other.set_hamiltonian(h_other)
    ic_other.gs_energy()

    scaled_other = idmrg.PeriodicMPS(
        ic_other._result.sites_uc, n_uc,
        [ITensor(T.inds, T.array * 0.9) for T in ic_other._result.U_list],
        eta=0.81)

    result = idmrg.imps_sum(ic_dom._result, scaled_other, cutoff=1e-12, maxdim=None)
    assert result.eta == pytest.approx(1.0, abs=1e-6)
    for p in range(n_uc):
        orig = idmrg.onsite_expectation(ic_dom._result, "Sz", p)
        new = idmrg.onsite_expectation(result, "Sz", p)
        assert new == pytest.approx(orig, abs=1e-6)
        chi_dom = idmrg._to_array_lpr(ic_dom._result.U_list[p]).shape[0]
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
    ic_half.maxm, ic_half.maxiter, ic_half.etol = 4, 20, 1e-8
    ic_half.set_hamiltonian(ic_half.SxC[0] * ic_half.SxR[0]
                             + ic_half.SyC[0] * ic_half.SyR[0]
                             + ic_half.SzC[0] * ic_half.SzR[0])
    ic_half.gs_energy()

    ic_one = infinitechain.Infinite_Spin_Chain(["1"])
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
    ic_down.maxm, ic_down.maxiter, ic_down.etol = 4, 50, 1e-12
    ic_down.set_hamiltonian(10.0 * ic_down.SzC[0])
    ic_down.gs_energy()

    raw = idmrg._periodic_direct_sum(ic_up._result.U_list, ic_down._result.U_list, 1)
    arr = raw[0].array
    chi_up = idmrg._to_array_lpr(ic_up._result.U_list[0]).shape[0]
    chi_down = idmrg._to_array_lpr(ic_down._result.U_list[0]).shape[0]
    assert arr.shape == (chi_up + chi_down, 2, chi_up + chi_down)
    a_up = idmrg._to_array_lpr(ic_up._result.U_list[0])
    a_down = idmrg._to_array_lpr(ic_down._result.U_list[0])
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
def test_td_dynamical_correlator_works_on_v3_backend():
    """itensor_version=3 is a supported backend for td_dynamical_correlator
    (Chain::td_dynamical_correlator_window, mpscpp3/chain_session.h) --
    superseded what used to be
    test_td_dynamical_correlator_rejects_non_python_backend before that
    port landed; see tests/test_idmrg_window_v3.py for the dedicated,
    physics-accuracy-focused v3 coverage (cross-checked against the
    "python" backend) this smoke test doesn't attempt to duplicate."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.maxm, ic.maxiter, ic.etol, ic.niter = 8, 30, 1e-6, 30
    # A plain single-coupling XX term leaves the transfer matrix's
    # dominant eigenvalue (near-)degenerate at this loose convergence
    # (confirmed directly: RuntimeError from
    # idmrg_dominant_right_fixed_point's own power iteration) -- the full
    # XXX Heisenberg coupling (matching
    # examples/idmrg/heisenberg_infinite_python_VS_v3/main.py and
    # tests/test_idmrg_window_v3.py's own model) doesn't have this
    # problem.
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0]
                        + ic.SzC[0] * ic.SzR[0])
    ks, es, Skw = ic.td_dynamical_correlator(
        "Sz", 0, "Sz", n_window=4, dt=0.05, nt=3, maxdim=20, cutoff=1e-10,
        niter=30, x_values=[-1, 0, 1], ks=np.linspace(-np.pi, np.pi, 3))
    assert Skw.shape == (len(ks), len(es))
    assert np.all(np.isfinite(Skw))
    assert ic._session3 is not None


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
    ic.maxm = 20
    ic.maxiter = 300
    ic.etol = 1e-12
    ic.niter = 150
    ic.set_hamiltonian(1.4 * ic.SzC[0] + ic.SxC[0] * ic.SxR[0])
    ic.gs_energy()

    es = np.linspace(-1, 6, 100)
    es_kpm, y_kpm = ic.kpm_finite("Sz", 0, "Sz", 1, n_window=20,
                                    window_chain_kwargs=dict(maxm=30, nsweeps=10),
                                    delta=0.3, es=es)

    from dmrgpy.pyitensor import idmrg_window
    from dmrgpy.timedependent import _fourier_transform_correlator
    ts, xs, S = idmrg_window.dynamical_correlator_td(
        ic._result, n_window=16, opname_A="Sz", opname_B="Sz", dt=0.05,
        nt=60, cutoff=1e-10, maxdim=60, niter=50, x_values=[1])
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
    J, h = 1.0, 3.0
    ic = _polarized_xx_chain(J, h, gs_method="vumps")
    for k in np.linspace(0, 2 * np.pi, 9, endpoint=False):
        exact = 2 * h - J * np.cos(k)
        got = ic.excitation_energies(k, n=1)[0]
        assert got == pytest.approx(exact, abs=1e-8)


def test_excitation_gap_matches_exact_minimum():
    """min_k (2h - J*cos(k)) = 2h - J, attained at k=0."""
    J, h = 1.0, 3.0
    ic = _polarized_xx_chain(J, h, gs_method="vumps")
    assert ic.excitation_gap() == pytest.approx(2 * h - J, abs=1e-8)


def test_excitation_environment_is_cached_across_calls():
    """A second excitation_energies call (different k) must reuse
    self._excitation_env (identity, not just an equal rebuild) -- it is
    expensive to build (a null-space computation plus dense linear
    solves)."""
    ic = _polarized_xx_chain(gs_method="vumps")
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
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] - 2.0 * ic.SzC[0])
    with pytest.raises(NotImplementedError):
        ic.excitation_gap()


def test_excitation_energies_itensor_version3_not_implemented():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
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
    with pytest.raises(RuntimeError):
        ic.local_excitation_gap()


def test_local_excitation_gap_itensor_version3_not_implemented():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
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
    h = (j_strong * (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1])
         + j_weak * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0]))
    ic.maxm, ic.maxiter, ic.etol = 10, 50, 1e-9
    ic.set_hamiltonian(h)
    with pytest.raises(NotImplementedError):
        ic.local_excitation_gap(window=1)


def test_local_excitation_gap_windowed_before_set_hamiltonian_raises():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    with pytest.raises(RuntimeError):
        ic.local_excitation_gap(window=1)


def test_local_excitation_gap_windowed_itensor_version3_not_implemented():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
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
