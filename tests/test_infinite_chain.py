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
    the gap this test tolerates. DMRG never seeds Lanczos from a product
    state here (see CLAUDE.md), so the *size* of that gap varies
    noticeably run to run at fixed (maxm, maxiter) -- observed anywhere
    from ~3e-3 to ~2e-2 across different random seeds at this exact
    configuration -- hence the generous tolerance below; it still
    confirms the identity holds in the right ballpark and with the right
    sign, which a genuine correlator bug would not survive.

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
