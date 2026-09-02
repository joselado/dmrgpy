"""Conserved-sector (quantum-number) targeting on the ED backend:
Many_Body_Chain.set_conserved_sector() with mode="ED".

The user-facing feature is the one test_sector_conservation.py checks on
mpscpp3 and test_sector_conservation_python.py on the pure-Python DMRG
backend -- same call, same quantum-number names, same units (Sz in
integer 2*Sz units). What differs is the mechanism, and it is the
simplest of the three: a conserved charge is diagonal in the ED product
basis, so a sector is a *set of basis states* and confining a calculation
to it is taking the corresponding submatrix of every assembled operator
(edtk/edchain.py). The full Hilbert space is still built on the way
there, so a sector buys a smaller eigenproblem, not a smaller
construction.

The reference here is deliberately *not* that same restriction called
twice. `ed_sector_energy` below assembles the Hamiltonian on the full
Hilbert space through the module-level multioperator.MO2matrix (which
reaches the ED object's single-site building blocks and never its
sector-aware entry points) and restricts it by hand, exactly as the other
two files' references do -- so an error in edtk/edchain.py's masking or
in its charge bookkeeping shows up as a disagreement rather than
cancelling out on both sides.
"""
import numpy as np
import pytest

from dmrgpy import bosonchain, fermionchain, parafermionchain, spinchain
from dmrgpy.multioperator import MO2matrix

TOL = 1e-8


def ed_sector_energy(chain, hamiltonian, charge, target):
    """Lowest eigenvalue of `hamiltonian` restricted to the sector where
    the (diagonal) operator `charge` equals `target`, by dense ED on the
    full Hilbert space."""
    ed = chain.get_ED_obj()
    H = np.array(MO2matrix(hamiltonian, ed).todense())
    N = np.array(MO2matrix(charge, ed).todense())
    assert np.allclose(N, np.diag(np.diag(N))), "charge operator is not diagonal"
    sel = np.abs(np.diag(N).real - target) < 1e-9
    assert sel.any(), "empty sector in the ED reference"
    return float(np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0])


def tV_chain(n, v=2.0):
    """Spinless t-V chain, forced onto the ED solver."""
    fc = fermionchain.Fermionic_Chain(n)
    fc.mode = "ED"
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + v * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    return fc, h, sum(fc.N)


def heisenberg_chain(n, sites=2):
    """Heisenberg chain written the way dmrgpy's examples write it --
    Sx.Sx+Sy.Sy+Sz.Sz -- so that no individual term conserves Sz and only
    their sum does, which is what the conservation check has to see."""
    sc = spinchain.Spin_Chain([sites] * n)
    sc.mode = "ED"
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
              + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    return sc, h, sum(sc.Sz)


@pytest.mark.parametrize("nf", [0, 1, 3, 4, 6, 8])
def test_fermion_number_sector_matches_ed(nf):
    """Every particle-number sector of an 8-site t-V chain. Most of these
    are not the global ground state, which is the point."""
    n = 8
    fc, h, number = tV_chain(n)
    fc.set_conserved_sector(Nf=nf)
    e = fc.gs_energy(mode="ED")
    assert e == pytest.approx(ed_sector_energy(fc, h, number, nf), abs=TOL)
    # the state really is in that sector, not merely close in energy
    assert fc.vev(number, mode="ED").real == pytest.approx(nf, abs=TOL)


def test_lowest_sector_reproduces_the_unconstrained_ground_state():
    """min over sectors == the ordinary, sector-free ground state, and
    switching sector targeting back off restores the same answer."""
    n = 8
    fc, h, _ = tV_chain(n)
    e_dense = fc.gs_energy(mode="ED")
    energies = []
    for nf in range(n + 1):
        fc.set_conserved_sector(Nf=nf)
        energies.append(fc.gs_energy(mode="ED"))
    assert min(energies) == pytest.approx(e_dense, abs=TOL)
    fc.set_conserved_sector()
    assert fc.gs_energy(mode="ED") == pytest.approx(e_dense, abs=TOL)


@pytest.mark.parametrize("two_sz", [-6, -2, 0, 2, 8])
def test_total_sz_sector_matches_ed(two_sz):
    """Total-Sz sectors of a spin-1/2 Heisenberg chain, in ITensor's
    integer 2*Sz units -- the same units the DMRG backends take, which is
    what the factor of two in Spin_Chain.get_sector_charge_operators is
    there for."""
    n = 8
    sc, h, sz = heisenberg_chain(n)
    sc.set_conserved_sector(Sz=two_sz)
    e = sc.gs_energy(mode="ED")
    assert e == pytest.approx(ed_sector_energy(sc, h, sz, two_sz / 2.0), abs=TOL)
    assert sc.vev(sz, mode="ED").real == pytest.approx(two_sz / 2.0, abs=TOL)


@pytest.mark.parametrize("two_sz", [0, 2, 4])
def test_spin_one_sector_matches_ed(two_sz):
    """A spin-1 chain: three local states, each carrying 2*Sz = +2/0/-2,
    so the units matter as well as the bookkeeping."""
    n = 5
    sc, h, sz = heisenberg_chain(n, sites=3)
    sc.set_conserved_sector(Sz=two_sz)
    assert sc.gs_energy(mode="ED") == pytest.approx(
        ed_sector_energy(sc, h, sz, two_sz / 2.0), abs=TOL)
    assert sc.vev(sz, mode="ED").real == pytest.approx(two_sz / 2.0, abs=TOL)


def test_ed_and_itensor_v3_agree_sector_by_sector():
    """The two implementations of the same public API, on the same chain:
    block-sparse QN tensors versus a submatrix of the ED basis."""
    n = 6
    for nf in range(n + 1):
        fc, h, number = tV_chain(n)
        fc.mode = None  # let DMRG answer this one
        fc.set_conserved_sector(Nf=nf)
        e_dmrg = fc.gs_energy()
        fc_ed, _, number_ed = tV_chain(n)
        fc_ed.set_conserved_sector(Nf=nf)
        assert fc_ed.gs_energy(mode="ED") == pytest.approx(e_dmrg, abs=1e-6)
        assert fc_ed.vev(number_ed, mode="ED").real == pytest.approx(nf, abs=TOL)


def test_excited_states_stay_inside_the_sector():
    """A sector confines the whole spectrum, not just its lowest state:
    the ED excited states must be the sector's, and in ascending order
    from the sector's ground state rather than the global one."""
    n = 6
    fc, h, number = tV_chain(n)
    fc.set_conserved_sector(Nf=2)
    es = fc.get_excited(n=4, mode="ED")
    assert es[0] == pytest.approx(ed_sector_energy(fc, h, number, 2), abs=TOL)
    assert np.all(np.diff(es) >= -TOL)
    _, wfs = fc.get_excited_states(n=4, mode="ED")
    for wf in wfs:
        assert wf.aMb(number, wf).real == pytest.approx(2.0, abs=1e-7)


def test_static_correlator_inside_a_sector():
    """A charge-conserving correlator is measured in the sector; the same
    numbers as the reference computed on the restricted eigenvector."""
    n = 6
    fc, h, number = tV_chain(n)
    fc.set_conserved_sector(Nf=3)
    fc.gs_energy(mode="ED")
    # <N_i N_j> summed over the chain must equal Nf^2 - <N> for these
    # densities only if correlations are trivial, so compare against a
    # direct restricted-ED evaluation instead
    ed = fc.get_ED_obj()
    H = np.array(MO2matrix(h, ed).todense())
    N = np.array(MO2matrix(number, ed).todense())
    sel = np.abs(np.diag(N).real - 3.0) < 1e-9
    w = np.linalg.eigh(H[np.ix_(sel, sel)])[1][:, 0]
    for i in range(n - 1):
        op = np.array(MO2matrix(fc.N[i] * fc.N[i + 1], ed).todense())
        ref = complex(np.conjugate(w) @ op[np.ix_(sel, sel)] @ w).real
        got = fc.vev(fc.N[i] * fc.N[i + 1], mode="ED").real
        assert got == pytest.approx(ref, abs=1e-7)


def test_interleaved_spinful_chain_gets_sz_only_from_ed():
    """The interleaved spinful chain is 2n *spinless* fermionic sites on
    the DMRG side, which know about Nf and nothing about spin -- so
    mode="ED" is the only way to target an Sz sector on it, and DMRG must
    say so instead of answering."""
    n = 3
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    fc.mode = "ED"
    h = 0
    for i in range(n - 1):
        for cdag, c in [(fc.Cdagup, fc.Cup), (fc.Cdagdn, fc.Cdn)]:
            h = h + cdag[i] * c[i + 1] + cdag[i + 1] * c[i]
    h = h + sum(4.0 * fc.Nup[i] * fc.Ndn[i] for i in range(n))
    fc.set_hamiltonian(h)
    number, sz = sum(fc.Ntot), sum(fc.Nup) + (-1) * sum(fc.Ndn)
    fc.set_conserved_sector(Nf=3, Sz=1)
    assert fc.gs_energy(mode="ED") == pytest.approx(
        ed_sector_energy(fc, h, number + 100 * sz, 3 + 100 * 1), abs=TOL)
    assert fc.vev(number, mode="ED").real == pytest.approx(3.0, abs=TOL)
    assert fc.vev(sz, mode="ED").real == pytest.approx(1.0, abs=TOL)
    # ...and DMRG refuses rather than answering with an Nf-only sector
    fc.mode = None
    with pytest.raises(NotImplementedError):
        fc.gs_energy()


def test_native_hubbard_conserves_a_chosen_subset():
    """A native Hubbard chain offers Nf and Sz independently: asking for
    Nf alone leaves Sz free, asking for both pins both."""
    n = 4
    hc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    hc.mode = "ED"
    h = 0
    for i in range(n - 1):
        for c, cdag in [(hc.Cup, hc.Cdagup), (hc.Cdn, hc.Cdagdn)]:
            h = h + cdag[i] * c[i + 1] + cdag[i + 1] * c[i]
    h = h + sum(4.0 * hc.Nup[i] * hc.Ndn[i] for i in range(n))
    hc.set_hamiltonian(h)
    number, sz = sum(hc.Ntot), sum(hc.Nup) + (-1) * sum(hc.Ndn)
    hc.set_conserved_sector(Nf=4)
    e_nf = hc.gs_energy(mode="ED")
    assert e_nf == pytest.approx(ed_sector_energy(hc, h, number, 4), abs=TOL)
    assert hc.vev(number, mode="ED").real == pytest.approx(4.0, abs=TOL)
    hc.set_conserved_sector(Nf=4, Sz=2)
    e_both = hc.gs_energy(mode="ED")
    assert hc.vev(sz, mode="ED").real == pytest.approx(2.0, abs=TOL)
    # a strictly smaller space, so its ground state cannot be lower
    assert e_both >= e_nf - TOL


def test_boson_number_sector():
    """Boson chains conserve Nb, by exactly the same route."""
    n = 4
    bc = bosonchain.Bosonic_Chain(n, maxnb=[3] * n)
    bc.mode = "ED"
    h = 0
    for i in range(n - 1):
        h = h + bc.Adag[i] * bc.A[i + 1] + bc.Adag[i + 1] * bc.A[i]
    bc.set_hamiltonian(h)
    number = sum(bc.N)
    for nb in [0, 1, 2, 3]:
        bc.set_conserved_sector(Nb=nb)
        assert bc.gs_energy(mode="ED") == pytest.approx(
            ed_sector_energy(bc, h, number, nb), abs=TOL)
        assert bc.vev(number, mode="ED").real == pytest.approx(nb, abs=TOL)


def test_a_charge_changing_operator_is_refused_not_silently_zeroed():
    """Restricting an operator to the sector is exact for a static
    expectation value and identically *zero* for one that changes the
    charge -- so a charge-changing operator has to raise, exactly as it
    does on the DMRG backends, rather than hand back a clean wrong zero
    (which a dynamical correlator of C or S+ would then report as a flat
    spectrum)."""
    n = 6
    sc, _, _ = heisenberg_chain(n)
    sc.set_conserved_sector(Sz=0)
    sc.gs_energy(mode="ED")
    with pytest.raises(ValueError, match="conserve"):
        sc.vev(sc.Sx[0], mode="ED")
    with pytest.raises(ValueError, match="conserve"):
        sc.vev(sc.Sx[0] * sc.Sx[1], mode="ED")  # S+S+ / S-S- pieces
    # ...while the conserving ones go through
    assert sc.vev(sc.Sz[0] * sc.Sz[1], mode="ED").real == pytest.approx(
        sc.vev(sc.Sz[n - 2] * sc.Sz[n - 1], mode="ED").real, abs=1e-7)
    fc, _, _ = tV_chain(n)
    fc.set_conserved_sector(Nf=3)
    fc.gs_energy(mode="ED")
    with pytest.raises(ValueError, match="conserve"):
        fc.vev(fc.C[0], mode="ED")


def test_a_non_conserving_hamiltonian_is_refused():
    """The one operator whose restriction would be silently *wrong*: a
    Hamiltonian with matrix elements between sectors is not the sector's
    Hamiltonian, and its submatrix's ground state is not an eigenstate of
    anything."""
    n = 4
    sc = spinchain.Spin_Chain([2] * n)
    sc.mode = "ED"
    sc.set_hamiltonian(sum(sc.Sx[i] for i in range(n)))  # a transverse field
    sc.set_conserved_sector(Sz=0)
    with pytest.raises(ValueError, match="conserve"):
        sc.gs_energy(mode="ED")


def test_invalid_sector_requests_raise():
    n = 6
    fc, _, _ = tV_chain(n)
    with pytest.raises(ValueError, match="empty"):
        fc.set_conserved_sector(Nf=99)  # unreachable target
    with pytest.raises(ValueError):
        fc.set_conserved_sector(Sz=0)  # spinless fermions carry no Sz
    with pytest.raises(ValueError):
        fc.set_conserved_sector(Nf=1.5)  # not an integer
    pf = parafermionchain.Parafermionic_Chain(4)
    pf.mode = "ED"
    with pytest.raises((ValueError, NotImplementedError)):
        pf.set_conserved_sector(Nf=2)  # conserves nothing at all


def test_an_unreachable_ed_only_sector_raises_at_the_first_calculation():
    """A quantum number only ED can represent is not validated when it is
    set (the ED object is built lazily), so an unreachable target has to
    surface at the first calculation rather than pass silently."""
    n = 3
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    fc.mode = "ED"
    h = sum(fc.Cdagup[i] * fc.Cup[i + 1] + fc.Cdagup[i + 1] * fc.Cup[i]
            for i in range(n - 1))
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Sz=99)  # accepted here...
    with pytest.raises(ValueError, match="empty"):
        fc.gs_energy(mode="ED")  # ...and rejected here


def test_a_rejected_request_leaves_the_chain_untouched():
    """A request that turns out to be impossible must roll back, not leave
    the chain half switched into a sector every later call would fail in."""
    n = 8
    fc, h, number = tV_chain(n)
    fc.set_conserved_sector(Nf=3)
    e3 = fc.gs_energy(mode="ED")
    with pytest.raises(ValueError):
        fc.set_conserved_sector(Nf=99)
    assert fc.conserved_sector == {"Nf": 3}
    assert fc.gs_energy(mode="ED") == pytest.approx(e3, abs=TOL)


def test_site_and_pair_entropies_work_inside_a_sector():
    """The projector-based reduced density matrix (densitymatrix.py) builds
    products like S+S- that change the charge, so under a sector some of
    its entries are unrepresentable -- and are exactly zero anyway in a
    fixed-charge state, which is what its `in_sector` branch relies on.
    That branch was written against the DMRG backends' ValueError; ED
    raises the same one, so the entropies have to come out finite rather
    than propagating the error."""
    n = 6
    sc, _, sz = heisenberg_chain(n)
    sc.set_conserved_sector(Sz=2)
    wf = sc.get_gs(mode="ED")
    assert wf.aMb(sz, wf).real == pytest.approx(1.0, abs=1e-7)
    site = np.array([complex(sc.get_site_entropy(wf, i)).real
                     for i in range(n)])
    assert np.all(np.isfinite(site))
    assert np.all(site > -1e-9) and np.all(site < np.log(2.0) + 1e-9)
    mi = complex(sc.get_mutual_information(wf, 0, 1)).real
    assert np.isfinite(mi) and mi > -1e-9
