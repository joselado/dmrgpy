"""Coverage for the ITensor v3 C++ port of VUMPS's apply_mpo -- the
mixed-gauge MPO application (`mpscpp3/chain_session.h`'s
`Chain::vumps_apply_mpo`, plus `Chain::vumps_load_uniform_state`/
`Chain::vumps_get_snapshot`) -- mirrors `tests/test_vumps_apply_mpo.py`'s
own `itensor_version="python"` coverage of `pyitensor/vumps.py`'s
apply_mpo (same D=1 exact closed-form case, same D>1 unitary invariants on
TFIM, same n_uc=2 bond-growth case), plus direct cross-backend agreement
against `itensor_version="python"` -- the C++ port is a from-scratch
dense-array translation sharing no code with pyitensor/vumps.py beyond the
algorithm itself, so this is a genuine cross-check, not a shared-bug blind
spot (same reasoning as test_vumps_correlator_v3.py's own docstring).

There is no `Infinite_Many_Body_Chain`-level wrapper for apply_mpo on
EITHER backend (matches pyitensor's own scope, see docs/user_guide.md's
own "Applying an operator/gate to a converged VUMPS iMPS" section) -- v3's
own apply_mpo is reached directly via `ic._session3.vumps_apply_mpo(...)`,
working against `ic._session3`'s own converged VUMPS snapshot.

Skipped automatically if mpscpp3 isn't compiled, exactly like
test_vumps_v3.py's own tests.
"""
import numpy as np
import pytest
from scipy.linalg import expm

from dmrgpy import cppext
from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps

pytestmark = pytest.mark.skipif(
    not cppext.available(3), reason="requires the compiled mpscpp3 (ITensor v3) extension")

_PAULI_X = np.array([[0, 1], [1, 0]], dtype=complex)
_SX = 0.5 * _PAULI_X
_SY = 0.5 * np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = 0.5 * np.array([[1, 0], [0, -1]], dtype=complex)


def _field_chain(B, maxm=1, itensor_version=3):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "vumps"
    ic.maxm = maxm
    ic.set_hamiltonian(-B * ic.SzC[0])
    ic.gs_energy()
    return ic


def _tfim_chain(g, maxm, itensor_version=3):
    """H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i], same convention as
    test_vumps_correlator_v3.py's own _tfim_chain."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "vumps"
    ic.maxm = maxm
    ic.etol = 1e-11
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


def _apply_pauli_x_v3(ic):
    """Applies a chi_W=1 Pauli-X (unitary, single-sublattice) to `ic`'s own
    converged VUMPS snapshot and LOADS the result back into
    `ic._session3` -- so ic.vev/ic.correlator (Python-level dispatch,
    unaffected by any of this since ic._session3_has_vumps stays True) see
    the flipped state afterward, mirroring how the pyitensor test suite
    passes the returned UniformMPS directly to vumps.onsite_expectation/
    two_point_correlator."""
    W = [_PAULI_X.flatten().tolist()]
    D, d_g, AL, AR, C, AC, eta = ic._session3.vumps_apply_mpo(W, [1], [1], 1e-12, 0)
    ic._session3.vumps_load_uniform_state(
        D, d_g, AL.flatten().tolist(), AR.flatten().tolist(), C.flatten().tolist())
    return eta


def test_apply_mpo_d1_pauli_x_flip_exact():
    """D=1 field-polarized product state: applying Pauli-X flips <Sz>
    exactly and preserves the trivially-factorized correlator, with eta
    staying exactly 1 (unitary operator, no truncation)."""
    ic = _field_chain(5.0)
    assert ic.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)

    eta = _apply_pauli_x_v3(ic)
    assert eta == pytest.approx(1.0, abs=1e-8)
    assert ic.vev("Sz", 0) == pytest.approx(-0.5, abs=1e-8)
    for r in (1, 2, 3):
        assert ic.correlator("Sz", 0, "Sz", r) == pytest.approx(0.25, abs=1e-8)


def test_apply_mpo_matches_pyitensor_apply_mpo():
    """The SAME chi_W=1 Pauli-X flip, applied independently on v3
    (Chain::vumps_apply_mpo) and pyitensor (vumps.apply_mpo) at the same
    exact D=1 field-polarized point, must agree -- no convergence-
    tolerance ambiguity to hide a real discrepancy behind."""
    ic_v3 = _field_chain(5.0, itensor_version=3)
    ic_py = _field_chain(5.0, itensor_version="python")

    _apply_pauli_x_v3(ic_v3)
    from dmrgpy.pyitensor.index import Index
    from dmrgpy.pyitensor.tensor import ITensor
    sites_uc = ic_py._vumps_result.sites_uc
    d = sites_uc.dim(1)
    link_l, link_r = Index(1, tags="Link"), Index(1, tags="Link")
    s = sites_uc.si(1)
    W_py = [ITensor((link_l, s, s.prime(1), link_r), _PAULI_X.reshape(1, d, d, 1))]
    flipped_py = vumps.apply_mpo(ic_py._vumps_result, W_py, cutoff=1e-12, maxdim=None)

    sz_v3 = ic_v3.vev("Sz", 0)
    sz_py = vumps.onsite_expectation(flipped_py, "Sz", 0)
    assert sz_v3 == pytest.approx(sz_py, abs=1e-8)


@pytest.mark.parametrize("D", [2, 3])
def test_apply_mpo_unitary_invariants_entangled(D):
    """Genuinely entangled D>1 (TFIM): a unitary chi_W=1 Pauli-X flip must
    flip <Sz>, leave <Sz(0)Sz(r)> unchanged, and preserve the norm
    diagnostic eta."""
    ic = _tfim_chain(1.5, maxm=D)
    sz0 = ic.vev("Sz", 0)
    c0 = [ic.correlator("Sz", 0, "Sz", r) for r in (1, 2, 3)]

    eta = _apply_pauli_x_v3(ic)
    assert eta == pytest.approx(1.0, abs=1e-6)
    assert ic.vev("Sz", 0) == pytest.approx(-sz0, abs=1e-6)
    for r, c0r in zip((1, 2, 3), c0):
        assert ic.correlator("Sz", 0, "Sz", r) == pytest.approx(c0r, abs=1e-6)


def test_apply_mpo_bond_growth_two_site_gate_n_uc2():
    """A genuinely bond-dimension>1 local gate (an SVD-split 2-site
    unitary rotation, tiled once per an n_uc=2 unit cell) -- exercises the
    actual bond-growth path (vumps_group_automaton's n_uc=2 branch, plus
    vx_grow_by_mpo_n1's Kronecker-merge); being unitary must preserve the
    norm (apply_mpo's own eta diagnostic)."""
    ic2 = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=3)
    ic2.gs_method = "vumps"
    h2 = (ic2.SxC[0] * ic2.SxC[1] + ic2.SyC[0] * ic2.SyC[1] + ic2.SzC[0] * ic2.SzC[1]
          + 0.4 * (ic2.SxC[1] * ic2.SxR[0] + ic2.SyC[1] * ic2.SyR[0] + ic2.SzC[1] * ic2.SzR[0]))
    ic2.set_hamiltonian(h2)
    ic2.maxm, ic2.maxiter, ic2.etol = 3, 100, 1e-10
    ic2.gs_energy()

    d = 2
    H2 = (np.kron(_SX, _SX) + np.kron(_SY, _SY) + np.kron(_SZ, _SZ)).real
    gate = expm(-1j * 0.37 * H2)
    gate4 = np.transpose(gate.reshape(d, d, d, d), (2, 0, 3, 1))
    U, S, Vh = np.linalg.svd(gate4.reshape(d * d, d * d), full_matrices=False)
    keep = int(np.sum(S > 1e-12))
    U, S, Vh = U[:, :keep], S[:keep], Vh[:keep, :]
    a_half = (U * S[None, :] ** 0.5).reshape(d, d, keep)
    b_half = (S[:, None] ** 0.5 * Vh).reshape(keep, d, d)

    W0 = a_half.reshape(1, d, d, keep).flatten().tolist()
    W1 = b_half.reshape(keep, d, d, 1).flatten().tolist()

    D, d_g, AL, AR, C, AC, eta = ic2._session3.vumps_apply_mpo(
        [W0, W1], [1, keep], [keep, 1], 1e-10, 0)

    assert eta == pytest.approx(1.0, abs=1e-4)
    assert d_g == 4


def test_apply_mpo_rejects_dimension_mismatch():
    """W_bulk's own grouped physical dimension must match the Chain's own
    vumps snapshot d_g -- feeding an operator built for a DIFFERENT
    (larger) local Hilbert space must raise rather than silently
    corrupting the contraction. Caught C++-side (Chain::vumps_apply_mpo's
    own ITError, surfaced by pybind11 as a plain RuntimeError) -- asserted
    as `Exception`, mirroring test_vumps_correlator_v3.py's own convention
    for this backend's input-validation paths."""
    ic = _field_chain(5.0)
    pauli_1 = np.eye(3, dtype=complex)  # a spin-1 (d=3) operator, wrong d_g=2
    with pytest.raises(Exception):
        ic._session3.vumps_apply_mpo([pauli_1.flatten().tolist()], [1], [1], 1e-12, 0)


def test_apply_mpo_requires_converged_snapshot():
    """Called on a fresh Chain that never ran vumps_ground_state/
    vumps_load_uniform_state -- must raise, not silently read garbage
    (mirrors test_idmrg_window_v3.py's own
    test_td_dynamical_correlator_v3_requires_gs_energy_first pattern)."""
    backend = cppext.get_backend(3)
    chain = backend.Chain([2])  # site_types: one spin-1/2
    with pytest.raises(Exception):
        chain.vumps_apply_mpo([_PAULI_X.flatten().tolist()], [1], [1], 1e-12, 0)
