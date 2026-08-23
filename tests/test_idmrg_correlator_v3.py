"""Coverage for the ITensor v3 C++ port of iDMRG's static-correlator and
local-gap machinery (`mpscpp3/chain_session.h`'s
`Chain::idmrg_onsite_expectation`/`Chain::idmrg_two_point_correlator`/
`Chain::idmrg_local_excitation_gap`, wired into
`Infinite_Many_Body_Chain.vev`/`correlator`/`local_excitation_gap` via
`itensor_version=3`, `gs_method="idmrg"`).

Mirrors what `tests/test_infinite_chain.py` already covers for
`itensor_version="python"` (`pyitensor/idmrg.py`'s own `onsite_expectation`/
`two_point_correlator`/`local_excitation_gap`), plus direct cross-backend
agreement. The C++ side is an independent dense-array translation of the
same algorithm -- it shares no code with `pyitensor/idmrg.py` -- so a
direct value comparison between the two is a genuine cross-check rather
than a shared-bug blind spot, exactly as `tests/test_vumps_correlator_v3.py`
argues for its own VUMPS counterpart.

The two exact-value anchors here are model-level, not implementation-level:
a field-polarized product state (whose every static observable is known in
closed form) and the translational-invariance identity `<H_uc> = n_uc*e0`,
which holds for any translationally invariant state independently of
whether `e0` itself is known in closed form. That identity is what actually
guards the *gauge* of the extracted unit cell: tiling the growing
algorithm's raw per-micro-step factors instead of the theta cell
(`Chain::idmrg_theta_cell`) violates it by orders of magnitude while
leaving the energy itself correct -- see `pyitensor/idmrg.py`'s own
`_theta_cell`/`IDMRGResult` docstrings for the measured numbers.

Skipped automatically if mpscpp3 isn't compiled, exactly like
`test_vumps_correlator_v3.py`'s own tests.
"""
import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy import infinitechain

pytestmark = pytest.mark.skipif(
    not cppext.available(3), reason="requires the compiled mpscpp3 (ITensor v3) extension")


def _polarized_xx_chain(itensor_version, J=1.0, h=3.0, maxm=4, maxiter=50,
                         etol=1e-12):
    """A field-polarized n_uc=1 XX chain (h above the XX chain's own
    saturation field) -- the ground state is the exact fully-polarized
    product state, so every static observable is known in closed form:
    <Sz>=1/2, <Sz(0)Sz(r)>=1/4 for every r, <Sx(0)Sx(r)>=0 for r>=1.
    Same model as test_infinite_chain.py's own `_polarized_xx_chain`."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    h_op = -J * (ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0]) - 2 * h * ic.SzC[0]
    ic.maxm, ic.maxiter, ic.etol = maxm, maxiter, etol
    ic.set_hamiltonian(h_op)
    ic.gs_energy()
    return ic


def _heisenberg_chain(itensor_version, n_uc, maxm=30, maxiter=40, etol=1e-9):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc,
                                            itensor_version=itensor_version)
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
    ic.maxm, ic.maxiter, ic.etol = maxm, maxiter, etol
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


def _dimerized_chain(itensor_version, j_strong=1.0, j_weak=0.2, maxm=8,
                     maxiter=150, etol=1e-14):
    """The same gapped, moderately dimerized n_uc=2 chain
    test_infinite_chain.py's own `test_unit_cell_expectation_self_consistency`
    uses (including its etol=1e-14 "run the full maxiter" trick, see that
    test's own docstring)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"],
                                            itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    h = (j_strong * (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1]
                     + ic.SzC[0] * ic.SzC[1])
         + j_weak * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0]
                     + ic.SzC[1] * ic.SzR[0]))
    ic.maxm, ic.maxiter, ic.etol = maxm, maxiter, etol
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


def _ssh_fermion_chain(itensor_version, t1=1.0, t2=0.4, maxm=20, maxiter=60,
                        etol=1e-11):
    """A dimerized (SSH-like) spinless-fermion chain, n_uc=2 -- the
    smallest model whose *correlators* need a real Jordan-Wigner string
    (a <Cdag(0) C(2)> skips over the site in between). Site code 0 is a
    spinless fermion site, see examples/idmrg/fermionic_infinite_chain."""
    ic = infinitechain.Infinite_Many_Body_Chain([0, 0],
                                                 itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    cdag = lambda g, i: ic.get_operator("Cdag", i, group=g)
    c = lambda g, i: ic.get_operator("C", i, group=g)
    h = (t1 * (cdag("C", 0) * c("C", 1) + cdag("C", 1) * c("C", 0))
         + t2 * (cdag("C", 1) * c("R", 0) + cdag("R", 0) * c("C", 1)))
    ic.maxm, ic.maxiter, ic.etol = maxm, maxiter, etol
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


# == Exact closed-form checks: field-polarized product state ================

def test_polarized_onsite_matches_exact():
    ic = _polarized_xx_chain(3)
    assert ic.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)


@pytest.mark.parametrize("r", [0, 1, 2, 5])
def test_polarized_zz_correlator_factorizes(r):
    """<Sz(0)Sz(r)> = 1/4 for every r in a fully polarized state,
    including r=0 (where Sz^2 = 1/4 for spin-1/2)."""
    ic = _polarized_xx_chain(3)
    assert ic.correlator("Sz", 0, "Sz", r) == pytest.approx(0.25, abs=1e-7)


@pytest.mark.parametrize("r", [1, 2, 3])
def test_polarized_transverse_correlator_vanishes(r):
    ic = _polarized_xx_chain(3)
    assert ic.correlator("Sx", 0, "Sx", r) == pytest.approx(0.0, abs=1e-7)


def test_polarized_same_site_mixed_operator_order():
    """r=0 composes as M_j @ M_i (the last-written factor is the left
    matrix factor -- see Chain::idmrg_correlator_endpoints, shared with
    the VUMPS path). For a fully Sz-polarized state,
    <Sx(0)Sy(0)> = <(i/2)Sz> = i/4, whereas the opposite order would give
    -i/4 -- so this pins the ordering convention, which a same-name pair
    like ("Sz","Sz") cannot."""
    ic = _polarized_xx_chain(3)
    assert ic.correlator("Sx", 0, "Sy", 0) == pytest.approx(0.25j, abs=1e-7)


# == Gauge check: <H_uc> == n_uc * e0 =======================================

@pytest.mark.parametrize("n_uc", [1, 2])
def test_unit_cell_expectation_self_consistency_v3(n_uc):
    """The identity this whole port exists to satisfy: <H_uc> must equal
    n_uc*e0 by translational invariance. It is a *linear*-order property of
    the converged state, so -- unlike the energy, which is quadratic in the
    wavefunction error -- it only holds once the extracted unit cell is
    genuinely gauge-consistent."""
    ic = _dimerized_chain(3) if n_uc == 2 else _heisenberg_chain(3, 1, maxm=20,
                                                                  maxiter=80)
    couplings = ([1.0, 0.2] if n_uc == 2 else [1.0])
    total = 0.0
    for i in range(n_uc):
        for name in ("Sx", "Sy", "Sz"):
            total += couplings[i] * ic.correlator(name, i, name, 1).real
    assert total == pytest.approx(n_uc * ic.e0, abs=5e-2)


# == Direct cross-check against itensor_version="python" ====================

def _tfim_chain(itensor_version, J=1.0, h=2.0, maxm=16, maxiter=120, etol=1e-12):
    """A gapped (paramagnetic-phase, h>J) n_uc=1 transverse-field Ising
    chain, H = -4J*Sx_i*Sx_{i+1} - 2h*Sz_i -- same convention/model as
    test_infinite_chain.py's own `_tfim_chain`. Gapped and
    non-degenerate, so both backends converge to the *same* state rather
    than to two different members of a degenerate manifold."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    h_op = -4.0 * J * (ic.SxC[0] * ic.SxR[0]) - 2.0 * h * ic.SzC[0]
    ic.maxm, ic.maxiter, ic.etol = maxm, maxiter, etol
    ic.set_hamiltonian(h_op)
    ic.gs_energy()
    return ic


@pytest.mark.parametrize("n_uc", [1, 2])
def test_correlator_matches_python_backend_on_gapped_chain(n_uc):
    """Same model, same algorithm, two independent implementations.

    Deliberately gapped models (TFIM in its paramagnetic phase for n_uc=1,
    a dimerized Heisenberg chain for n_uc=2) rather than the *critical*
    uniform Heisenberg chain: on a gapless chain at finite maxm the
    growing algorithm settles onto a slightly, arbitrarily
    symmetry-broken state whose residual depends on its own (unseeded)
    random starting MPS -- measured at <Sz> ~ +4e-4 on one backend against
    ~ -3e-4 on the other at maxm=30/maxiter=40, both a legitimate
    approximation to the exact 0, but not equal to each other. That is a
    property of the model, not of either implementation; comparing two
    backends' *states* needs a model whose state is unique."""
    ic_v3 = _tfim_chain(3) if n_uc == 1 else _dimerized_chain(3)
    ic_py = _tfim_chain("python") if n_uc == 1 else _dimerized_chain("python")
    assert ic_v3.e0 == pytest.approx(ic_py.e0, abs=1e-6)
    for p in range(n_uc):
        assert ic_v3.vev("Sz", p) == pytest.approx(ic_py.vev("Sz", p), abs=1e-5)
    for r in (0, 1, 2, 3):
        assert (ic_v3.correlator("Sz", 0, "Sz", r)
                == pytest.approx(ic_py.correlator("Sz", 0, "Sz", r), abs=1e-5))


def test_fermionic_correlator_matches_python_backend():
    """<Cdag(0) C(r)> on the SSH chain, for r=1 (adjacent, no string) and
    r=2/r=3 (a genuine Jordan-Wigner string across the sites in between)
    -- the string threading is the part of the correlator path that has no
    analogue in a purely bosonic/spin model."""
    ic_v3 = _ssh_fermion_chain(3)
    ic_py = _ssh_fermion_chain("python")
    assert ic_v3.e0 == pytest.approx(ic_py.e0, abs=1e-6)
    assert ic_v3.vev("N", 0) == pytest.approx(ic_py.vev("N", 0), abs=1e-6)
    for r in (1, 2, 3):
        assert (ic_v3.correlator("Cdag", 0, "C", r)
                == pytest.approx(ic_py.correlator("Cdag", 0, "C", r), abs=1e-6))


def test_fermionic_half_filling_v3():
    """The SSH chain above is particle-hole symmetric at zero chemical
    potential, so it sits at half filling: <N> = 1/2 on both sublattices,
    exactly, independent of (t1, t2)."""
    ic = _ssh_fermion_chain(3)
    for p in (0, 1):
        assert ic.vev("N", p) == pytest.approx(0.5, abs=1e-6)


def test_odd_parity_operator_pair_rejected():
    """A pair with odd total fermion parity has a Jordan-Wigner string
    that can never close -- rejected Python-side before either backend
    sees it (Chain::idmrg_correlator_endpoints would also reject it)."""
    ic = _ssh_fermion_chain(3)
    with pytest.raises(ValueError):
        ic.correlator("Cdag", 0, "N", 2)


# == local_excitation_gap ===================================================

def test_local_excitation_gap_matches_python_backend():
    """The local superblock gap is a property of the growing algorithm's
    own final 2-site effective Hamiltonian, so the two backends only agree
    to the extent that they converge to the same one -- checked on the
    gapped, D=1 polarized chain, where that Hamiltonian is essentially
    exactly reproducible."""
    ic_v3 = _polarized_xx_chain(3)
    ic_py = _polarized_xx_chain("python")
    assert ic_v3.local_excitation_gap() == pytest.approx(
        ic_py.local_excitation_gap(), rel=1e-3)


def test_local_excitation_gap_is_positive_and_near_exact_xx_gap():
    """Exact single-magnon dispersion of the polarized XX chain is
    E(k) = 2h - J*cos(k), minimized at 2h-J = 5.0 for the defaults here.
    local_excitation_gap is a deliberately cruder estimator (no momentum
    label, frozen environments -- see its own docstring), hence the loose
    tolerance; this checks it lands in the right place at all."""
    ic = _polarized_xx_chain(3)
    gap = ic.local_excitation_gap()
    assert gap > 0.0
    assert gap == pytest.approx(5.0, rel=0.15)


# == Error paths ============================================================

def test_vev_rejects_p_out_of_range():
    ic = _polarized_xx_chain(3)
    with pytest.raises(ValueError):
        ic.vev("Sz", 1)


def test_correlator_rejects_negative_r():
    ic = _polarized_xx_chain(3)
    with pytest.raises(ValueError):
        ic.correlator("Sz", 0, "Sz", -1)


def test_correlator_before_set_hamiltonian_raises():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "idmrg"
    with pytest.raises(RuntimeError):
        ic.correlator("Sz", 0, "Sz", 1)


def test_switching_gs_method_after_vumps_run_still_works():
    """Mirror of test_vumps_correlator_v3.py's own
    test_switching_gs_method_after_idmrg_run_still_works: a self._session3
    left over from a gs_method="vumps" run must not be reused by the
    gs_method="idmrg" correlator path (its Chain has no iDMRG unit-cell
    snapshot at all), it must transparently rerun gs_energy()."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.maxm = 1
    ic.set_hamiltonian(-1.0 * (ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0])
                        - 6.0 * ic.SzC[0])
    ic.gs_method = "vumps"
    ic.gs_energy()
    ic.gs_method = "idmrg"
    assert ic.vev("Sz", 0) == pytest.approx(0.5, abs=1e-6)
