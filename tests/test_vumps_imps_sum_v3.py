"""Coverage for the ITensor v3 C++ port of VUMPS's imps_sum -- the direct
sum of two converged mixed-gauge VUMPS iMPS (`mpscpp3/chain_session.h`'s
`Chain::vumps_imps_sum`, plus `Chain::vumps_load_uniform_state`/
`Chain::vumps_get_snapshot`) -- mirrors `tests/test_vumps_imps_sum.py`'s
own `itensor_version="python"` coverage of `pyitensor/vumps.py`'s
imps_sum: the degenerate-dominant-eigenvalue RuntimeError for two ordinary
(both eta=1) VUMPS ground states, and the well-posed norm-mismatch case at
D=1 (exact) and genuinely entangled D>1 (TFIM), where the surviving
branch's onsite_expectation/two_point_correlator are checked directly.

There is no `Infinite_Many_Body_Chain`-level wrapper for imps_sum on
EITHER backend (matches pyitensor's own scope) -- v3's own imps_sum is
reached directly via `ic._session3.vumps_imps_sum(...)`.

Skipped automatically if mpscpp3 isn't compiled, exactly like
test_vumps_v3.py's own tests.
"""
import pytest

from dmrgpy import cppext
from dmrgpy import infinitechain

pytestmark = pytest.mark.skipif(
    not cppext.available(3), reason="requires the compiled mpscpp3 (ITensor v3) extension")


def _field_chain(B, maxm=1):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "vumps"
    ic.maxm = maxm
    ic.set_hamiltonian(-B * ic.SzC[0])
    ic.gs_energy()
    return ic


def _tfim_chain(g, maxm):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "vumps"
    ic.maxm = maxm
    ic.etol = 1e-11
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


def _rescaled_AL_b(ic, factor):
    """`ic`'s own converged AL, rescaled by `factor` -- the v3 analogue of
    test_vumps_imps_sum.py's own `_rescaled_copy` (Chain::vumps_imps_sum
    only ever reads AL_b off the second operand, exactly mirroring
    pyitensor/vumps.py's own imps_sum, which only reads .AL off
    result_b -- see that function's own construction: `raw[Da:,:,Da:] =
    result_b.AL`)."""
    D_b, d_g, AL, _AR, _C = ic._session3.vumps_get_snapshot()
    return D_b, (factor * AL).flatten().tolist()


def test_imps_sum_tied_norm_raises():
    """Two separately-converged, ordinary VUMPS ground states: both are
    individually normalized to eta=1 exactly (mixed-gauge construction),
    so the combined transfer matrix's dominant eigenvalue is exactly
    degenerate -- must raise (Chain::vx_canonicalize_n1's own
    vx_dominant_*_fixed_point degeneracy guard), not silently collapse to
    one of the two (physically distinct, oppositely Sz-polarized)
    branches."""
    ic_up = _field_chain(10.0)
    ic_down = _field_chain(-10.0)
    D_b, d_g, AL_b, _AR, _C = ic_down._session3.vumps_get_snapshot()

    with pytest.raises(Exception):
        ic_up._session3.vumps_imps_sum(D_b, AL_b.flatten().tolist(), 1e-12, 0)


def test_imps_sum_dominant_branch_survives_d1():
    """A well-posed (non-degenerate) case at D=1: sum a converged
    field-polarized product state with a deliberately norm-rescaled copy
    (amplitude x0.9, i.e. self-overlap eta=0.81) of a DIFFERENT converged
    product state. The smaller-norm branch must be exponentially
    suppressed: onsite_expectation of the summed-and-truncated result must
    reproduce the dominant branch's own exact closed-form value (0.5), and
    the surviving bond dimension must collapse back down to D=1."""
    ic_dom = _field_chain(1.0)
    ic_other = _field_chain(-0.6)
    D_b, AL_b = _rescaled_AL_b(ic_other, 0.9)

    D, d_g, AL, AR, C, AC, eta = ic_dom._session3.vumps_imps_sum(D_b, AL_b, 1e-12, 0)

    assert eta == pytest.approx(1.0, abs=1e-6)
    assert D == 1

    ic_dom._session3.vumps_load_uniform_state(
        D, d_g, AL.flatten().tolist(), AR.flatten().tolist(), C.flatten().tolist())
    assert ic_dom.vev("Sz", 0) == pytest.approx(0.5, abs=1e-6)


@pytest.mark.parametrize("D", [2, 3])
def test_imps_sum_dominant_branch_survives_entangled(D):
    """The same well-posed construction at a genuinely entangled D>1
    (TFIM, gapped) -- exercises vumps_complete_mixed_gauge's eigh-based C
    construction on a non-trivial (rank>1) fixed point, not just the D=1
    scalar case above. Checks BOTH onsite_expectation and a multi-r
    two_point_correlator scan against the untouched dominant branch's own
    values."""
    ic_dom = _tfim_chain(1.5, maxm=D)
    ic_other = _tfim_chain(2.0, maxm=D)

    sz_dom = ic_dom.vev("Sz", 0)
    corr_dom = [ic_dom.correlator("Sz", 0, "Sz", r) for r in (1, 2, 3)]

    D_b, AL_b = _rescaled_AL_b(ic_other, 0.9)
    D_new, d_g, AL, AR, C, AC, eta = ic_dom._session3.vumps_imps_sum(D_b, AL_b, 1e-12, 0)

    assert eta == pytest.approx(1.0, abs=1e-6)

    ic_dom._session3.vumps_load_uniform_state(
        D_new, d_g, AL.flatten().tolist(), AR.flatten().tolist(), C.flatten().tolist())
    assert ic_dom.vev("Sz", 0) == pytest.approx(sz_dom, abs=1e-6)
    for r, c_dom in zip((1, 2, 3), corr_dom):
        assert ic_dom.correlator("Sz", 0, "Sz", r) == pytest.approx(c_dom, abs=1e-6)


def test_imps_sum_rejects_dimension_mismatch():
    """Different local Hilbert-space dimensions (spin-1/2 vs spin-1) give
    different grouped physical dimensions d_g -- a meaningless sum, must
    raise rather than silently placing mismatched-size physical legs into
    the same combined tensor. Caught C++-side (AL_b's own flat size no
    longer matches D_b*d_g*D_b for this Chain's own d_g) -- asserted as
    `Exception`, same convention as test_vumps_apply_mpo_v3.py's own
    dimension-mismatch test."""
    ic_half = _field_chain(1.0)

    ic_one = infinitechain.Infinite_Spin_Chain(["1"], itensor_version=3)
    ic_one.gs_method = "vumps"
    ic_one.maxm = 1
    ic_one.set_hamiltonian(-1.0 * ic_one.SzC[0])
    ic_one.gs_energy()

    D_b, d_g_b, AL_b, _AR, _C = ic_one._session3.vumps_get_snapshot()
    assert d_g_b == 3  # spin-1, vs ic_half's own d_g=2

    with pytest.raises(Exception):
        ic_half._session3.vumps_imps_sum(D_b, AL_b.flatten().tolist(), 1e-12, 0)


def test_imps_sum_requires_converged_snapshot():
    """Called on a fresh Chain that never ran vumps_ground_state/
    vumps_load_uniform_state -- must raise, not silently read garbage."""
    backend = cppext.get_backend(3)
    chain = backend.Chain([2])
    with pytest.raises(Exception):
        chain.vumps_imps_sum(1, [1.0 + 0j, 0j, 0j, 1.0 + 0j], 1e-12, 0)
