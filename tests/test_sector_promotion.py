"""Leaving a conserved sector while keeping the state computed inside it:
Many_Body_Chain.promote_to_dense() / promote_mps().

The point of the feature is the operations a sector *forbids*. While
set_conserved_sector() is active every operator built on the chain must
conserve the sector's quantum numbers -- a bare C changes Nf, Sx changes
Sz -- and that is enforced rather than merely discouraged, because
ITensor's AutoMPO aborts the process over a flux-violating term instead
of raising anything catchable. Promotion converts the block-sparse state
to its exact dense equivalent, after which all of those are legal again.

Every check here is against ED restricted to the sector *exactly* (a
submatrix in the ED product basis), for the same reason
test_sector_conservation.py spells out: filtering full-spectrum
eigenvectors by <N> instead is not equivalent, because sector-degenerate
eigenvalues come back as arbitrary superpositions.

These tests run on itensor_version=3; itensor_version="python"
implements the same promotion API (an index relabeling rather than a
block scatter, see pyitensor/sector.py) and is covered by
test_sector_promotion_python.py.
"""
import numpy as np
import pytest

from dmrgpy import bosonchain, fermionchain, spinchain
from dmrgpy.multioperator import MO2matrix

DMRG_TOL = 1e-6


def tV_chain(n, v=2.0):
    """Spinless t-V chain plus its total-number operator."""
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + v * fc.N[i] * fc.N[i + 1]
    return fc, h, sum(fc.N)


def heisenberg_chain(n):
    sc = spinchain.Spin_Chain([2] * n)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    return sc, h, sum(sc.Sz)


def ed_sector_ground_state(chain, hamiltonian, charge, target):
    """(energy, statevector) of `hamiltonian` restricted to the sector where
    the diagonal operator `charge` equals `target`, with the vector embedded
    back into the full Hilbert space so any operator -- including ones that
    leave the sector -- can be sandwiched with it."""
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
    """<v|op|v> for a full-space ED vector `v`."""
    M = np.array(MO2matrix(op, chain.get_ED_obj()).todense())
    return complex(np.conjugate(v) @ (M @ v))


def solved_in_sector(n=8, nf=3, v=2.0):
    """A t-V chain solved at fixed particle number, plus its ED reference."""
    fc, h, number = tV_chain(n, v=v)
    fc.set_hamiltonian(h)
    fc.maxm = 40
    fc.nsweeps = 8
    fc.set_conserved_sector(Nf=nf)
    e = fc.gs_energy()
    e_ed, v_ed = ed_sector_ground_state(fc, h, number, nf)
    assert e == pytest.approx(e_ed, abs=DMRG_TOL)
    return fc, h, number, v_ed


def test_promotion_keeps_the_sector_state():
    """The state survives the promotion unchanged: same energy, same sector,
    same site densities as ED restricted to that sector."""
    n, nf = 8, 3
    fc, h, number, v_ed = solved_in_sector(n=n, nf=nf)
    e_before = fc.gs_energy()
    fc.promote_to_dense()
    assert fc.conserved_sector is None
    # gs_energy() is cached, so check the energy through the state itself
    assert fc.vev(h).real == pytest.approx(e_before, abs=DMRG_TOL)
    assert fc.vev(number).real == pytest.approx(nf, abs=1e-8)
    for i in range(n):
        assert fc.vev(fc.N[i]).real == pytest.approx(
            ed_vev(fc, fc.N[i], v_ed).real, abs=DMRG_TOL)


def test_annihilation_operator_is_rejected_in_the_sector_and_legal_after():
    """The actual motivation: c_i applied to a fixed-N ground state."""
    n, nf, site = 8, 3, 2
    fc, h, number, v_ed = solved_in_sector(n=n, nf=nf)
    with pytest.raises(Exception):  # sector_terms() names the offending term
        fc.applyoperator(fc.C[site], fc.wf0)
    fc.promote_to_dense()
    wf = fc.applyoperator(fc.C[site], fc.wf0)
    # |c_i|gs>|^2 = <gs|n_i|gs> exactly -- an identity that needs no reference
    assert wf.norm() ** 2 == pytest.approx(
        fc.vev(fc.N[site]).real, abs=DMRG_TOL)
    # and the result really has one particle fewer
    assert fc.vev(number, wf=wf).real == pytest.approx(nf - 1, abs=1e-6)


def test_one_body_density_matrix_from_applied_operators_matches_ed():
    """<c_i gs|c_j gs> = <gs|c^dag_i c_j|gs>, computed by actually applying
    the (sector-violating) operators to the promoted state."""
    n, nf = 8, 3
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
    sc, h, sz = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.maxm = 40
    sc.nsweeps = 8
    sc.set_conserved_sector(Sz=two_sz)
    sc.gs_energy()
    e_ed, v_ed = ed_sector_ground_state(sc, h, sz, two_sz / 2.0)
    with pytest.raises(Exception):
        sc.vev(sc.Sx[0] * sc.Sx[1])
    sc.promote_to_dense()
    for j in range(1, n):
        got = sc.vev(sc.Sx[0] * sc.Sx[j]).real
        assert got == pytest.approx(ed_vev(sc, sc.Sx[0] * sc.Sx[j], v_ed).real,
                                    abs=DMRG_TOL)


def test_promote_mps_converts_a_handle_taken_before_promotion():
    """promote_to_dense() only reaches the chain's own ground state; an MPS
    Python already holds needs promote_mps()."""
    n, nf = 8, 3
    fc, h, number, v_ed = solved_in_sector(n=n, nf=nf)
    held = fc.get_gs().copy()
    fc.promote_to_dense()
    promoted = fc.promote_mps(held)
    assert fc.vev(number, wf=promoted).real == pytest.approx(nf, abs=1e-8)
    assert fc.overlap(promoted, fc.wf0) == pytest.approx(1.0, abs=DMRG_TOL)
    # idempotent: promoting an already-dense state changes nothing
    again = fc.promote_mps(promoted)
    assert fc.overlap(again, promoted) == pytest.approx(1.0, abs=DMRG_TOL)


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
    fc.maxm = 40
    fc.nsweeps = 8
    fc.set_conserved_sector(Nf=3)
    e3 = fc.gs_energy()
    fc.promote_to_dense()
    fc.set_conserved_sector(Nf=5)
    e5 = fc.gs_energy()
    assert fc.vev(number).real == pytest.approx(5, abs=1e-8)
    _, v5 = ed_sector_ground_state(fc, h, number, 5)
    assert e5 == pytest.approx(ed_vev(fc, h, v5).real, abs=DMRG_TOL)
    assert e5 != pytest.approx(e3, abs=1e-3)


def test_states_promoted_from_different_sectors_are_comparable():
    """Two sectors solved one after the other on the same chain, each
    promoted, end up on the *same* dense site indices -- so they can be
    overlapped, which is what a photoemission weight
    Z_i = |<gs_{N-1}|c_i|gs_N>|^2 needs. This is why promotion rebases onto
    the chain's own site set (kept from construction) rather than onto a
    freshly minted one, which would give every promotion new index ids and
    make two promoted states silently incomparable."""
    n, nf = 8, 4
    fc, h, number = tV_chain(n)
    fc.set_hamiltonian(h)
    fc.maxm = 40
    fc.nsweeps = 8

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


def test_hubbard_promotion_matches_ed_on_a_spin_flip_observable():
    """Native Hubbard sites (four local states, Nf and Sz both conserved)
    exercise the one assumption the conversion makes: that the QN-carrying
    and dense site indices enumerate their local basis states in the same
    order. A spin-flip operator c^dag_up c_dn breaks Sz, so it cannot even
    be built inside the sector -- and would read garbage if the local basis
    ordering did not survive the promotion."""
    n = 4
    hc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = 0
    for i in range(n - 1):
        for c, cdag in [(hc.Cup, hc.Cdagup), (hc.Cdn, hc.Cdagdn)]:
            h = h + cdag[i] * c[i + 1] + cdag[i + 1] * c[i]
    h = h + sum(4.0 * hc.Nup[i] * hc.Ndn[i] for i in range(n))
    number = sum(hc.Ntot)
    sz = sum(hc.Nup[i] - hc.Ndn[i] for i in range(n))

    hc.set_hamiltonian(h)
    hc.maxm = 40
    hc.nsweeps = 8
    hc.set_conserved_sector(Nf=n, Sz=0)
    e = hc.gs_energy()

    ed = hc.get_ED_obj()
    H = np.array(MO2matrix(h, ed).todense())
    N = np.array(MO2matrix(number, ed).todense())
    S = np.array(MO2matrix(sz, ed).todense())
    sel = (np.abs(np.diag(N).real - n) < 1e-9) & (np.abs(np.diag(S).real) < 1e-9)
    es, vs = np.linalg.eigh(H[np.ix_(sel, sel)])
    v = np.zeros(H.shape[0], dtype=complex)
    v[sel] = vs[:, 0]
    assert e == pytest.approx(float(es[0]), abs=DMRG_TOL)

    flip = hc.Cdagup[0] * hc.Cdn[1]
    with pytest.raises(Exception):
        hc.vev(flip)
    hc.promote_to_dense()
    for i in range(n):
        for j in range(n):
            op = hc.Cdagup[i] * hc.Cdn[j]
            got = hc.vev(op)
            ref = ed_vev(hc, op, v)
            assert got.real == pytest.approx(ref.real, abs=DMRG_TOL)
            assert got.imag == pytest.approx(ref.imag, abs=DMRG_TOL)
    # occupations too, one per spin species, as an ordering check on the
    # four-state local basis
    for i in range(n):
        for op in (hc.Nup[i], hc.Ndn[i]):
            assert hc.vev(op).real == pytest.approx(ed_vev(hc, op, v).real,
                                                    abs=DMRG_TOL)


def test_boson_number_sector_promotion_matches_ed():
    """Boson sites are dmrgpy's own site type (extra/bosonfour.h), not an
    upstream ITensor one, so their QN-vs-dense local basis ordering is
    worth checking separately. A bare `A` changes Nb and is refused inside
    the sector."""
    n, nb = 4, 3
    bc = bosonchain.Bosonic_Chain(n, maxnb=[3] * n)
    h = 0
    for i in range(n - 1):
        h = h + bc.Adag[i] * bc.A[i + 1] + bc.Adag[i + 1] * bc.A[i]
    h = h + sum(1.5 * bc.N[i] * bc.N[i] for i in range(n))
    number = sum(bc.N)

    bc.set_hamiltonian(h)
    bc.maxm = 40
    bc.nsweeps = 8
    bc.set_conserved_sector(Nb=nb)
    e = bc.gs_energy()
    e_ed, v_ed = ed_sector_ground_state(bc, h, number, nb)
    assert e == pytest.approx(e_ed, abs=DMRG_TOL)

    with pytest.raises(Exception):
        bc.applyoperator(bc.A[0], bc.wf0)
    bc.promote_to_dense()
    for i in range(n):
        assert bc.vev(bc.N[i]).real == pytest.approx(
            ed_vev(bc, bc.N[i], v_ed).real, abs=DMRG_TOL)
    wf = bc.applyoperator(bc.A[0], bc.wf0)
    assert wf.norm() ** 2 == pytest.approx(bc.vev(bc.N[0]).real, abs=DMRG_TOL)
    assert bc.vev(number, wf=wf).real == pytest.approx(nb - 1, abs=1e-6)


def test_dynamical_correlator_of_a_charge_changing_operator_after_promotion():
    """The user-guide claim that a dynamical correlator with a
    sector-violating vertex becomes available after promoting. The
    reference is the *same* calculation on a chain that was dense from the
    start, at a filling whose sector ground state is the global one -- not
    ED, whose exact Lorentzian broadening differs from KPM's by far more
    than any promotion error."""
    n, nf = 6, 3  # half filling: the global ground state of a hopping chain

    def hopping_chain():
        fc = fermionchain.Fermionic_Chain(n)
        h = 0
        for i in range(n - 1):
            h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        fc.set_hamiltonian(h)
        fc.maxm = 40
        fc.nsweeps = 8
        return fc

    sector = hopping_chain()
    sector.set_conserved_sector(Nf=nf)
    e_sector = sector.gs_energy()
    with pytest.raises(Exception):  # the vertex leaves the sector
        sector.get_dynamical_correlator(name=(sector.C[0], sector.Cdag[0]),
                                        mode="DMRG", delta=0.2)
    sector.promote_to_dense()

    dense = hopping_chain()
    assert dense.gs_energy() == pytest.approx(e_sector, abs=DMRG_TOL)

    w = np.linspace(-3.5, 3.5, 301)

    def spectrum(chain):
        x, y = chain.get_dynamical_correlator(
            name=(chain.C[0], chain.Cdag[0]), mode="DMRG", delta=0.2)
        return np.interp(w, x, y.real)

    got, ref = spectrum(sector), spectrum(dense)
    assert np.max(np.abs(got - ref)) < 1e-6 * max(1.0, np.max(np.abs(ref)))
