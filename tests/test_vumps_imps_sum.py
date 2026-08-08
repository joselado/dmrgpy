"""Coverage for pyitensor/vumps.py's imps_sum/UniformMPS/_complete_mixed_gauge
-- the VUMPS-mixed-gauge analogue of pyitensor.idmrg.imps_sum, see vumps.py's
own "Summing two converged VUMPS iMPS" section docstring for the construction
and physical scope (identical to idmrg.imps_sum's own scope, see
tests/test_infinite_chain.py's own imps_sum tests for the idmrg-side
counterparts these mirror).

Cross-checks used: (1) the same degenerate-dominant-eigenvalue RuntimeError
idmrg.imps_sum's own tests rely on, here triggered by summing two ordinary
(both eta=1 on both the left AND right transfer eigenvalue, by mixed-gauge
construction) VUMPSResults; (2) the well-posed norm-mismatch case, exercised
at D=1 (field-polarized product states, exact closed form) and at genuinely
entangled D>1 (TFIM, gapped) where the surviving branch's onsite_expectation
AND two_point_correlator (built from the summed-and-reconstructed AC/AR) are
checked directly against the untouched dominant branch's own values -- unlike
the D=1 case, this also exercises `_complete_mixed_gauge`'s eigh-based C
construction on a genuinely rank->1 (non-scalar) fixed point, not just a
trivial 1x1 case.
"""
import pytest

from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps


def _field_chain(B, maxm=1):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = maxm
    ic.set_hamiltonian(-B * ic.SzC[0])
    return ic


def _tfim_chain(g, maxm):
    """H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i] (same convention as
    test_vumps_correlator.py's own _tfim_chain) -- gapped away from
    criticality (g!=1), so VUMPS converges reliably at modest D."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = maxm
    ic.etol = 1e-11
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    return ic


def _rescaled_copy(result, factor):
    """A UniformMPS identical to `result`'s own physical content but with
    every mixed-gauge tensor rescaled by `factor` (so the ket amplitude
    scales by `factor`, and the self-overlap eta by `factor**2`) -- the
    VUMPS analogue of test_infinite_chain.py's own scaled_other
    (idmrg.PeriodicMPS with .array*0.9), used to manufacture a genuine
    per-site norm mismatch since ordinary VUMPSResults never carry one on
    their own (AL always exactly left-canonical, AR always exactly
    right-canonical, by mixed-gauge construction -- eta=1 on both sides).
    AC scales by factor**2 (== factor(AL) * factor(C)), matching AL@C's own
    bilinearity, so the returned UniformMPS's {AL,AR,C,AC} stay mutually
    consistent (AL@C == C@AR == AC) despite AL/AR no longer being exactly
    iso/right-canonical (deliberately, to represent an eta=factor**2 branch,
    exactly mirroring idmrg.imps_sum's own worked example)."""
    return vumps.UniformMPS(
        result.sites_uc, result.n_uc, result.D, result.d_g,
        result.AL * factor, result.AR * factor,
        result.C * factor, result.AC * factor ** 2,
        eta=factor ** 2)


# -- pyitensor.vumps.imps_sum: direct sum of two converged VUMPS iMPS -- see
# vumps.py's own "Summing two converged VUMPS iMPS" section docstring for the
# construction (block-diagonal direct sum of AL, re-canonicalized/truncated
# via idmrg._canonicalize_periodic, then completed to the full mixed gauge
# via _complete_mixed_gauge) and its physical scope: tiled to the
# thermodynamic limit, this only has a well-posed single-branch reduction
# when the two summands' own self-overlap eigenvalues (eta) differ --
# summing two *ordinary* VUMPSResults (always eta=1 on both sides, by
# mixed-gauge construction) hits a genuine tie every time, which
# idmrg._dominant_right_fixed_point's own degeneracy check (reused
# unmodified inside imps_sum's own _canonicalize_periodic call) must catch
# and raise on, not silently resolve to one arbitrary branch.

def test_imps_sum_tied_norm_raises():
    """The common case: summing two separately-converged, ordinary
    VUMPSResults. Both are individually normalized to eta=1 exactly on
    both the left and right transfer eigenvalue (mixed-gauge construction),
    so the combined transfer matrix's dominant eigenvalue is exactly
    degenerate -- must raise RuntimeError rather than silently collapsing
    to one of the two (physically distinct, oppositely Sz-polarized)
    branches."""
    ic_up = _field_chain(10.0)
    ic_up.gs_energy()
    ic_down = _field_chain(-10.0)
    ic_down.gs_energy()

    with pytest.raises(RuntimeError):
        vumps.imps_sum(ic_up._vumps_result, ic_down._vumps_result)


def test_imps_sum_dominant_branch_survives_d1():
    """A well-posed (non-degenerate) case at D=1: sum a converged
    field-polarized product state with a deliberately norm-rescaled copy
    of a *different* converged product state (amplitude x0.9, i.e.
    self-overlap eta=0.81 rather than the ordinary 1). The smaller-norm
    branch must be exponentially suppressed: onsite_expectation of the
    summed-and-truncated result must reproduce the dominant branch's own
    exact closed-form value (0.5), and the surviving bond dimension must
    collapse back down to D=1 (the rescaled branch fully discarded, not
    merely down-weighted)."""
    ic_dom = _field_chain(1.0)
    ic_dom.gs_energy()
    ic_other = _field_chain(-0.6)
    ic_other.gs_energy()

    scaled_other = _rescaled_copy(ic_other._vumps_result, 0.9)
    result = vumps.imps_sum(ic_dom._vumps_result, scaled_other, cutoff=1e-12, maxdim=None)

    assert result.eta == pytest.approx(1.0, abs=1e-6)
    assert result.D == ic_dom._vumps_result.D == 1
    assert vumps.onsite_expectation(result, "Sz", 0) == pytest.approx(0.5, abs=1e-6)


@pytest.mark.parametrize("D", [2, 3])
def test_imps_sum_dominant_branch_survives_entangled(D):
    """The same well-posed construction at a genuinely entangled D>1 (TFIM,
    gapped, away from criticality) -- exercises _complete_mixed_gauge's
    eigh-based C construction on a non-trivial (rank>1) fixed point, not
    just the D=1 scalar case above. Checks BOTH onsite_expectation and a
    multi-r two_point_correlator scan against the untouched dominant
    branch's own values, computed independently (same AC/AR-based formula,
    but on the freshly-reconstructed mixed gauge of the summed-and-
    truncated result, not a copy of the original tensors)."""
    ic_dom = _tfim_chain(1.5, maxm=D)
    ic_dom.gs_energy()
    ic_other = _tfim_chain(2.0, maxm=D)
    ic_other.gs_energy()

    scaled_other = _rescaled_copy(ic_other._vumps_result, 0.9)
    result = vumps.imps_sum(ic_dom._vumps_result, scaled_other, cutoff=1e-12, maxdim=None)

    assert result.eta == pytest.approx(1.0, abs=1e-6)

    sz_dom = vumps.onsite_expectation(ic_dom._vumps_result, "Sz", 0)
    sz_new = vumps.onsite_expectation(result, "Sz", 0)
    assert sz_new == pytest.approx(sz_dom, abs=1e-6)

    for r in (1, 2, 3):
        corr_dom = vumps.two_point_correlator(ic_dom._vumps_result, "Sz", 0, "Sz", r)
        corr_new = vumps.two_point_correlator(result, "Sz", 0, "Sz", r)
        assert corr_new == pytest.approx(corr_dom, abs=1e-6)


def test_imps_sum_rejects_n_uc_mismatch():
    ic1 = _field_chain(1.0)
    ic1.gs_energy()

    ic2 = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
    ic2.gs_method = "vumps"
    ic2.maxm = 1
    h = ic2.SxC[0] * ic2.SxC[1] + ic2.SyC[0] * ic2.SyC[1] + ic2.SzC[0] * ic2.SzC[1]
    ic2.set_hamiltonian(h)
    ic2.gs_energy()

    with pytest.raises(ValueError):
        vumps.imps_sum(ic1._vumps_result, ic2._vumps_result)


def test_imps_sum_rejects_dimension_mismatch():
    """Mirrors idmrg.imps_sum's own dimension-mismatch check: different
    local Hilbert-space dimensions (spin-1/2 vs spin-1) give different
    grouped physical dimensions d_g -- a meaningless sum, must raise
    ValueError rather than silently placing mismatched-size physical legs
    into the same combined tensor."""
    ic_half = _field_chain(1.0)
    ic_half.gs_energy()

    ic_one = infinitechain.Infinite_Spin_Chain(["1"], itensor_version="python")
    ic_one.gs_method = "vumps"
    ic_one.maxm = 1
    ic_one.set_hamiltonian(-1.0 * ic_one.SzC[0])
    ic_one.gs_energy()

    with pytest.raises(ValueError):
        vumps.imps_sum(ic_half._vumps_result, ic_one._vumps_result)
