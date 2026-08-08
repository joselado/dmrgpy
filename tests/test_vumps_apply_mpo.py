"""Coverage for pyitensor/vumps.py's apply_mpo -- the VUMPS-mixed-gauge
analogue of pyitensor.idmrg.apply_mpo, see vumps.py's own "Applying a
(bounded) MPO to a converged VUMPS iMPS" section docstring for the
construction (group W_bulk via _group_automaton, grow AL via
idmrg.grow_by_mpo on the trivial n_uc=1 "periodic chain" VUMPS already
works at, re-canonicalize/truncate via idmrg._canonicalize_periodic, then
complete the mixed gauge via _complete_mixed_gauge) and the W_bulk
convention (identical to idmrg.apply_mpo's own -- a list of n_uc ITensors,
one per unit-cell sublattice position, each rank-4 (Left, in, out, Right)).

Cross-checks used: (1) an exact D=1 closed-form case (field-polarized
product state, Pauli-X flip); (2) the SAME W_bulk list fed to both
idmrg.apply_mpo and vumps.apply_mpo on the same (exact D=1) model, checked
to agree; (3) unitary-operator invariants at genuinely entangled D>1 (TFIM)
-- <Sz> flips sign, <Sz(0)Sz(r)> is unchanged, and the norm diagnostic eta
stays at 1; (4) a genuinely bond-growing chi_W>1 two-site gate tiled once
per a n_uc=2 unit cell, checked for the same unitarity invariant (eta~1);
(5) the dimension-mismatch guard this module's own apply_mpo adds on top of
idmrg.apply_mpo's shared machinery.
"""
import numpy as np
import pytest
from scipy.linalg import expm

from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor import vumps
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor


def _field_chain_vumps(B, maxm=1):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = maxm
    ic.set_hamiltonian(-B * ic.SzC[0])
    ic.gs_energy()
    return ic


def _field_chain_idmrg(B, maxm=1):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "idmrg"
    ic.maxm = maxm
    ic.set_hamiltonian(-B * ic.SzC[0])
    ic.gs_energy()
    return ic


def _tfim_chain(g, maxm):
    """H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i] (same convention as
    test_vumps_imps_sum.py's own _tfim_chain) -- gapped away from
    criticality, so VUMPS converges reliably at modest D."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = maxm
    ic.etol = 1e-11
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


def _pauli_x_mpo(sites_uc):
    """chi_W=1 W_bulk = [Pauli-X] (a single-sublattice unitary, applied at
    every site) -- exactly idmrg's own apply_mpo example
    (examples/idmrg/apply_mpo_to_infinite_chain/main.py)."""
    d = sites_uc.dim(1)
    pauli_x = 2 * sites_uc.site_type(1).matrix("Sx")  # unitary for spin-1/2
    link_l, link_r = Index(1, tags="Link"), Index(1, tags="Link")
    s = sites_uc.si(1)
    return [ITensor((link_l, s, s.prime(1), link_r), pauli_x.reshape(1, d, d, 1))]


def test_apply_mpo_d1_pauli_x_flip_exact():
    """D=1 field-polarized product state: applying Pauli-X flips <Sz>
    exactly and preserves the trivially-factorized correlator (a product
    state has <Sz(0)Sz(r)>=<Sz>^2 for every r, unaffected by a purely local
    unitary), with eta staying exactly 1 (unitary operator, no truncation)."""
    ic = _field_chain_vumps(5.0)
    W = _pauli_x_mpo(ic._vumps_result.sites_uc)

    flipped = vumps.apply_mpo(ic._vumps_result, W, cutoff=1e-12, maxdim=None)

    sz0 = vumps.onsite_expectation(ic._vumps_result, "Sz", 0)
    sz0_flipped = vumps.onsite_expectation(flipped, "Sz", 0)
    assert sz0 == pytest.approx(0.5, abs=1e-8)
    assert sz0_flipped == pytest.approx(-0.5, abs=1e-8)
    assert flipped.eta == pytest.approx(1.0, abs=1e-8)
    assert flipped.D == 1

    for r in (1, 2, 3):
        c0 = vumps.two_point_correlator(ic._vumps_result, "Sz", 0, "Sz", r)
        c1 = vumps.two_point_correlator(flipped, "Sz", 0, "Sz", r)
        assert c1 == pytest.approx(c0, abs=1e-8)


def test_apply_mpo_matches_idmrg_apply_mpo():
    """The SAME W_bulk list, fed to both idmrg.apply_mpo and
    vumps.apply_mpo on the same exact D=1 model, must agree -- both are
    independent implementations of "apply this MPO to the converged iMPS",
    and D=1 field-polarization is exact on both backends, so there is no
    convergence-tolerance ambiguity to hide a real discrepancy behind."""
    ic_v = _field_chain_vumps(5.0)
    ic_i = _field_chain_idmrg(5.0)
    W = _pauli_x_mpo(ic_v._vumps_result.sites_uc)

    flipped_v = vumps.apply_mpo(ic_v._vumps_result, W, cutoff=1e-12, maxdim=None)
    flipped_i = idmrg.apply_mpo(ic_i._result, W, cutoff=1e-12, maxdim=None)

    sz_v = vumps.onsite_expectation(flipped_v, "Sz", 0)
    sz_i = idmrg.onsite_expectation(flipped_i, "Sz", 0)
    assert sz_v == pytest.approx(sz_i, abs=1e-8)
    assert flipped_v.eta == pytest.approx(flipped_i.eta, abs=1e-8)


@pytest.mark.parametrize("D", [2, 3])
def test_apply_mpo_unitary_invariants_entangled(D):
    """Genuinely entangled D>1 (TFIM): a unitary chi_W=1 Pauli-X flip must
    flip <Sz>, leave <Sz(0)Sz(r)> unchanged, and preserve the norm
    diagnostic eta -- exactly the invariants
    examples/idmrg/apply_mpo_to_infinite_chain/main.py already checks for
    idmrg.apply_mpo, here on vumps.apply_mpo's own mixed-gauge output."""
    ic = _tfim_chain(1.5, maxm=D)
    W = _pauli_x_mpo(ic._vumps_result.sites_uc)

    flipped = vumps.apply_mpo(ic._vumps_result, W, cutoff=1e-12, maxdim=None)

    assert flipped.eta == pytest.approx(1.0, abs=1e-6)
    sz0 = vumps.onsite_expectation(ic._vumps_result, "Sz", 0)
    sz0_flipped = vumps.onsite_expectation(flipped, "Sz", 0)
    assert sz0_flipped == pytest.approx(-sz0, abs=1e-6)

    for r in (1, 2, 3):
        c0 = vumps.two_point_correlator(ic._vumps_result, "Sz", 0, "Sz", r)
        c1 = vumps.two_point_correlator(flipped, "Sz", 0, "Sz", r)
        assert c1 == pytest.approx(c0, abs=1e-6)


def test_apply_mpo_bond_growth_two_site_gate_n_uc2():
    """A genuinely bond-dimension>1 local gate (an SVD-split 2-site unitary
    rotation, tiled once per an n_uc=2 unit cell) -- exercises the actual
    bond-growth path (_group_automaton's n_uc=2 einsum branch, plus
    idmrg.grow_by_mpo's Kronecker-merge), and being unitary must preserve
    the norm (apply_mpo's own eta diagnostic), mirroring
    examples/idmrg/apply_mpo_to_infinite_chain/main.py's own chi_W>1 case."""
    ic2 = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
    ic2.gs_method = "vumps"
    h2 = (ic2.SxC[0] * ic2.SxC[1] + ic2.SyC[0] * ic2.SyC[1] + ic2.SzC[0] * ic2.SzC[1]
          + 0.4 * (ic2.SxC[1] * ic2.SxR[0] + ic2.SyC[1] * ic2.SyR[0] + ic2.SzC[1] * ic2.SzR[0]))
    ic2.set_hamiltonian(h2)
    ic2.maxm, ic2.maxiter, ic2.etol = 3, 100, 1e-10
    ic2.gs_energy()

    sites_uc2 = ic2._vumps_result.sites_uc
    d = sites_uc2.dim(1)
    Sx, Sy, Sz = (sites_uc2.site_type(1).matrix(n) for n in ("Sx", "Sy", "Sz"))
    H2 = (np.kron(Sx, Sx) + np.kron(Sy, Sy) + np.kron(Sz, Sz)).real
    gate = expm(-1j * 0.37 * H2)
    gate4 = np.transpose(gate.reshape(d, d, d, d), (2, 0, 3, 1))
    U, S, Vh = np.linalg.svd(gate4.reshape(d * d, d * d), full_matrices=False)
    keep = int(np.sum(S > 1e-12))
    U, S, Vh = U[:, :keep], S[:keep], Vh[:keep, :]
    a_half = (U * S[None, :] ** 0.5).reshape(d, d, keep)
    b_half = (S[:, None] ** 0.5 * Vh).reshape(keep, d, d)

    s0, s1 = sites_uc2.si(1), sites_uc2.si(2)
    left_dummy, mid, right_dummy = (Index(1, tags="Link"), Index(keep, tags="Link"),
                                     Index(1, tags="Link"))
    W0 = ITensor((left_dummy, s0, s0.prime(1), mid), a_half.reshape(1, d, d, keep))
    W1 = ITensor((mid, s1, s1.prime(1), right_dummy), b_half.reshape(keep, d, d, 1))

    gated = vumps.apply_mpo(ic2._vumps_result, [W0, W1], cutoff=1e-10, maxdim=None)

    assert gated.eta == pytest.approx(1.0, abs=1e-4)
    assert gated.n_uc == 2
    assert gated.d_g == ic2._vumps_result.d_g


def test_apply_mpo_rejects_dimension_mismatch():
    """W_bulk's own grouped physical dimension must match result.d_g --
    feeding a Pauli-X built for a DIFFERENT (larger) local Hilbert space
    must raise ValueError rather than silently corrupting the contraction."""
    ic = _field_chain_vumps(5.0)

    ic_one = infinitechain.Infinite_Spin_Chain(["1"], itensor_version="python")
    ic_one.gs_method = "vumps"
    ic_one.maxm = 1
    ic_one.set_hamiltonian(-1.0 * ic_one.SzC[0])
    ic_one.gs_energy()
    W_wrong = _pauli_x_mpo(ic_one._vumps_result.sites_uc)

    with pytest.raises(ValueError):
        vumps.apply_mpo(ic._vumps_result, W_wrong, cutoff=1e-12, maxdim=None)
