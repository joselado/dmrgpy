"""Leaving a conserved sector while keeping the state computed inside it,
on the pure-Python backend: promote_to_dense() / promote_mps() with
itensor_version="python".

Same feature, same motivation and the same ED-restricted references as
test_sector_promotion.py checks on mpscpp3 -- the operations a sector
*forbids* (applying a bare C, an Sx correlator, a charge-changing
dynamical correlator vertex) become legal on a state that was nevertheless
obtained inside one.

What "promotion" converts is backend-specific: mpscpp3 turns a genuinely
block-sparse MPS into a dense one (removeQNs), while pyitensor's state was
dense all along and only its site-index *labels* change (see
pyitensor/sector.py). Both are exact and involve no re-solve, and the two
tests below that check a promoted state against ED on a sector-violating
observable are what actually verifies that -- an implementation that got
the local basis ordering wrong would read garbage there.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain
from dmrgpy.multioperator import MO2matrix

DMRG_TOL = 1e-6


def tV_chain(n, v=2.0):
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + v * fc.N[i] * fc.N[i + 1]
    fc.maxm, fc.nsweeps = 40, 8
    return fc, h, sum(fc.N)


def ed_sector_ground_state(chain, hamiltonian, charge, target):
    """(energy, statevector) of `hamiltonian` restricted to the sector
    where the diagonal operator `charge` equals `target`, with the vector
    embedded back into the full Hilbert space so any operator -- including
    ones that leave the sector -- can be sandwiched with it."""
    ed = chain.get_ED_obj()
    H = np.array(MO2matrix(hamiltonian, ed).todense())
    N = np.array(MO2matrix(charge, ed).todense())
    assert np.allclose(N, np.diag(np.diag(N))), "charge operator is not diagonal"
    sel = np.abs(np.diag(N).real - target) < 1e-9
    assert sel.any(), "empty sector in the ED reference"
    es, vs = np.linalg.eigh(H[np.ix_(sel, sel)])
    v = np.zeros(H.shape[0], dtype=complex)
    v[sel] = vs[:, 0]
    return float(es[0]), v


def ed_vev(chain, op, v):
    M = np.array(MO2matrix(op, chain.get_ED_obj()).todense())
    return complex(np.conjugate(v) @ (M @ v))


def solved_in_sector(n=8, nf=3, v=2.0):
    fc, h, number = tV_chain(n, v=v)
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=nf)
    e = fc.gs_energy()
    e_ed, v_ed = ed_sector_ground_state(fc, h, number, nf)
    assert e == pytest.approx(e_ed, abs=DMRG_TOL)
    return fc, h, number, v_ed


def test_promotion_keeps_the_sector_state():
    """The state survives the promotion unchanged: same energy, same
    sector, same site densities as ED restricted to that sector."""
    n, nf = 8, 3
    fc, h, number, v_ed = solved_in_sector(n=n, nf=nf)
    e_before = fc.gs_energy()
    fc.promote_to_dense()
    assert fc.conserved_sector is None
    assert fc.vev(h).real == pytest.approx(e_before, abs=DMRG_TOL)
    assert fc.vev(number).real == pytest.approx(nf, abs=1e-8)
    for i in range(n):
        assert fc.vev(fc.N[i]).real == pytest.approx(
            ed_vev(fc, fc.N[i], v_ed).real, abs=DMRG_TOL)


def test_annihilation_operator_is_rejected_in_the_sector_and_legal_after():
    """The actual motivation: c_i applied to a fixed-N ground state."""
    n, nf, site = 8, 3, 2
    fc, h, number, v_ed = solved_in_sector(n=n, nf=nf)
    with pytest.raises(ValueError):  # sector_terms() names the offending term
        fc.applyoperator(fc.C[site], fc.wf0)
    fc.promote_to_dense()
    wf = fc.applyoperator(fc.C[site], fc.wf0)
    # |c_i|gs>|^2 = <gs|n_i|gs> exactly -- an identity that needs no reference
    assert wf.norm() ** 2 == pytest.approx(fc.vev(fc.N[site]).real, abs=DMRG_TOL)
    assert fc.vev(number, wf=wf).real == pytest.approx(nf - 1, abs=1e-6)


def test_one_body_density_matrix_from_applied_operators_matches_ed():
    """<c_i gs|c_j gs> = <gs|c^dag_i c_j|gs>, computed by actually applying
    the (sector-violating) operators to the promoted state."""
    n, nf = 6, 3
    fc, h, number, v_ed = solved_in_sector(n=n, nf=nf)
    fc.promote_to_dense()
    applied = [fc.applyoperator(fc.C[i], fc.wf0) for i in range(n)]
    for i in range(n):
        for j in range(n):
            got = fc.overlap(applied[i], applied[j])
            ref = ed_vev(fc, fc.Cdag[i] * fc.C[j], v_ed)
            assert got.real == pytest.approx(ref.real, abs=DMRG_TOL)
            assert got.imag == pytest.approx(ref.imag, abs=DMRG_TOL)


def test_promoted_state_matches_ed_on_a_sector_violating_correlator():
    """A spin chain solved at Sz=0, then measured with Sx -- an operator
    that carries no definite Sz charge at all and so cannot even be built
    while the sector is set."""
    n, two_sz = 8, 0
    sc = spinchain.Spin_Chain([2] * n)
    sc.setup_python()
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sz = sum(sc.Sz)
    sc.maxm, sc.nsweeps = 40, 8
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=two_sz)
    sc.gs_energy()
    _, v_ed = ed_sector_ground_state(sc, h, sz, two_sz / 2.0)
    with pytest.raises(ValueError):
        sc.vev(sc.Sx[0] * sc.Sx[1])
    sc.promote_to_dense()
    for j in range(1, n):
        assert sc.vev(sc.Sx[0] * sc.Sx[j]).real == pytest.approx(
            ed_vev(sc, sc.Sx[0] * sc.Sx[j], v_ed).real, abs=DMRG_TOL)


def test_promote_mps_converts_a_handle_taken_before_promotion():
    """promote_to_dense() only reaches the chain's own ground state; an MPS
    Python already holds needs promote_mps()."""
    n, nf = 8, 3
    fc, h, number, v_ed = solved_in_sector(n=n, nf=nf)
    held = fc.get_gs().copy()
    fc.promote_to_dense()
    promoted = fc.promote_mps(held)
    assert fc.vev(number, wf=promoted).real == pytest.approx(nf, abs=1e-8)
    assert abs(fc.overlap(promoted, fc.wf0)) == pytest.approx(1.0, abs=DMRG_TOL)
    # idempotent: promoting an already-dense state changes nothing
    again = fc.promote_mps(promoted)
    assert abs(fc.overlap(again, promoted)) == pytest.approx(1.0, abs=DMRG_TOL)


def test_an_unpromoted_state_is_refused_rather_than_silently_contracted():
    """Handing back a state built under a different sector setting must
    raise. With dense storage the mismatched site indices would not fail
    on their own -- they would simply not contract, silently producing an
    outer product -- so the site set is checked explicitly."""
    n, nf = 6, 3
    fc, h, number, v_ed = solved_in_sector(n=n, nf=nf)
    held = fc.get_gs().copy()
    fc.promote_to_dense()
    with pytest.raises(ValueError, match="different site set"):
        fc.vev(number, wf=held)


def test_promotion_without_a_sector_is_a_noop():
    fc, h, number = tV_chain(6)
    fc.set_hamiltonian(h)
    e = fc.gs_energy()
    fc.promote_to_dense()
    assert fc.conserved_sector is None
    assert fc.gs_energy() == pytest.approx(e, abs=DMRG_TOL)


def test_a_sector_can_be_set_again_after_promoting():
    """Promotion clears the sector rather than disabling the machinery."""
    n = 8
    fc, h, number = tV_chain(n)
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=3)
    e3 = fc.gs_energy()
    fc.promote_to_dense()
    fc.set_conserved_sector(Nf=5)
    e5 = fc.gs_energy()
    assert fc.vev(number).real == pytest.approx(5, abs=1e-8)
    e_ed, _ = ed_sector_ground_state(fc, h, number, 5)
    assert e5 == pytest.approx(e_ed, abs=DMRG_TOL)
    assert e5 != pytest.approx(e3, abs=1e-3)


def test_states_promoted_from_different_sectors_are_comparable():
    """Two sectors solved one after the other on the same chain, each
    promoted, end up on the *same* dense site indices -- which is what a
    photoemission weight Z_i = |<gs_{N-1}|c_i|gs_N>|^2 needs. This is why
    promotion rebases onto the chain's own site set, kept from
    construction, rather than a freshly minted one."""
    n, nf = 6, 3
    fc, h, number = tV_chain(n)
    fc.set_hamiltonian(h)

    fc.set_conserved_sector(Nf=nf)
    fc.gs_energy()
    _, v_N = ed_sector_ground_state(fc, h, number, nf)
    fc.promote_to_dense()
    wf_N = fc.wf0.copy()

    fc.set_conserved_sector(Nf=nf - 1)
    fc.gs_energy()
    _, v_Nm = ed_sector_ground_state(fc, h, number, nf - 1)
    fc.promote_to_dense()
    wf_Nm = fc.wf0.copy()

    assert fc.vev(number, wf=wf_N).real == pytest.approx(nf, abs=1e-8)
    assert fc.vev(number, wf=wf_Nm).real == pytest.approx(nf - 1, abs=1e-8)

    ed = fc.get_ED_obj()
    for i in range(n):
        z = abs(fc.overlap(wf_Nm, fc.applyoperator(fc.C[i], wf_N))) ** 2
        C = np.array(MO2matrix(fc.C[i], ed).todense())
        z_ed = abs(np.conjugate(v_Nm) @ (C @ v_N)) ** 2
        assert z == pytest.approx(z_ed, abs=DMRG_TOL)
