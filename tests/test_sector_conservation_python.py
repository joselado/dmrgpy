"""Conserved-sector (quantum-number) targeting on the pure-Python backend:
Many_Body_Chain.set_conserved_sector() with itensor_version="python".

The user-facing feature is the same one test_sector_conservation.py checks
on mpscpp3, and the reference is the same: ED restricted to the sector
*exactly* (a submatrix in the ED product basis), never ED filtered by an
expectation value -- see that file's docstring for why the difference
matters.

What is not the same is the implementation, and one test here exists
specifically because of it. mpscpp3 confines the calculation structurally
(QN-block-sparse tensors: an amplitude outside the sector has nowhere to
be stored), while pyitensor keeps dense storage and adds a charge penalty
to the variational solves, because a dense LAPACK SVD leaks ~1e-16 across
charge blocks on every truncation and a sweep amplifies that leak toward
whichever sector is lower in energy. See pyitensor/sector.py, and
test_the_charge_penalty_is_load_bearing below, which turns the penalty off
and shows the state leaving the sector entirely.
"""
import numpy as np
import pytest

from dmrgpy import bosonchain, fermionchain, spinchain
from dmrgpy.multioperator import MO2matrix

DMRG_TOL = 1e-6


def ed_sector_energy(chain, hamiltonian, charge, target):
    """Lowest eigenvalue of `hamiltonian` restricted to the sector where
    the (diagonal) operator `charge` equals `target`, by dense ED."""
    ed = chain.get_ED_obj()
    H = np.array(MO2matrix(hamiltonian, ed).todense())
    N = np.array(MO2matrix(charge, ed).todense())
    assert np.allclose(N, np.diag(np.diag(N))), "charge operator is not diagonal"
    sel = np.abs(np.diag(N).real - target) < 1e-9
    assert sel.any(), "empty sector in the ED reference"
    return float(np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0])


def tV_chain(n, v=2.0):
    """Spinless t-V chain plus its total-number operator, on the
    pure-Python backend."""
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + v * fc.N[i] * fc.N[i + 1]
    fc.maxm = 40
    fc.nsweeps = 8
    return fc, h, sum(fc.N)


def heisenberg_chain(n):
    """Heisenberg chain written the way dmrgpy's examples write it --
    Sx.Sx+Sy.Sy+Sz.Sz -- which is exactly the case where no individual
    term conserves Sz and only their sum does (see pyitensor/sector.py's
    expand_xy_terms/combine_terms)."""
    sc = spinchain.Spin_Chain([2] * n)
    sc.setup_python()
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.maxm = 40
    sc.nsweeps = 8
    return sc, h, sum(sc.Sz)


@pytest.mark.parametrize("nf", [0, 1, 3, 4, 6, 8])
def test_fermion_number_sector_matches_ed(nf):
    """Every particle-number sector of an 8-site t-V chain. Most of these
    are *not* the global ground state, which is the whole point: a solver
    that leaks out of the sector would answer with a lower energy."""
    n = 8
    fc, h, number = tV_chain(n)
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=nf)
    e = fc.gs_energy()
    assert e == pytest.approx(ed_sector_energy(fc, h, number, nf), abs=DMRG_TOL)
    # the state really is in that sector, not merely close in energy
    assert fc.vev(number).real == pytest.approx(nf, abs=1e-8)


def test_the_charge_penalty_is_load_bearing():
    """Turning the charge penalty off must break confinement -- otherwise
    the penalty is dead weight and this backend's sector mode could be
    simplified away.

    The setup is deliberately the worst case for dense storage: an
    attractive interaction, so the global ground state is the *full* band,
    and a target sector far from it. The start state is in the sector
    either way (Chain._default_mps is sector-aware), so what this isolates
    is exactly the SVD leakage pyitensor/sector.py describes, not a bad
    starting point.
    """
    n = 12
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h - 1.0 * fc.N[i] * fc.N[i + 1]
    number = sum(fc.N)
    fc.maxm = 40
    fc.nsweeps = 20
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=2)

    e = fc.gs_energy()
    assert e == pytest.approx(ed_sector_energy(fc, h, number, 2), abs=DMRG_TOL)
    assert fc.vev(number).real == pytest.approx(2, abs=1e-8)

    fc.restart()
    fc.set_sector_penalty(0.0)  # python-backend-only knob
    with pytest.raises(RuntimeError, match="not in the requested sector"):
        fc.gs_energy()


@pytest.mark.parametrize("n", [2, 3, 4])
def test_short_chains_and_the_symmetric_start_state(n):
    """Short chains, where the in-sector start state has few enough terms
    to be an exact eigenstate by accident -- and was one.

    A 2-site hopping chain at Nf=1 has exactly two in-sector product
    states, and their *equal-weight* sum (|10>+|01>)/sqrt(2) is the +1
    eigenvector. A noise-free DMRG started there never leaves it, so the
    first version of this backend's sector mode returned the sector's
    highest state instead of its lowest. The start state's coefficients
    are random for exactly this reason (see sector.py's sector_mps), and
    this is the regression guard. Note these sizes also exercise the
    n<3 chains that itensor_version=3 cannot run at all.
    """
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    h = sum(fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
            for i in range(n - 1))
    number = sum(fc.N)
    fc.maxm, fc.nsweeps = 20, 8
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=1)
    # one particle hopping on an open chain: E = -2*cos(pi/(n+1))
    assert fc.gs_energy() == pytest.approx(
        ed_sector_energy(fc, h, number, 1), abs=DMRG_TOL)
    assert fc.vev(number).real == pytest.approx(1.0, abs=1e-8)


def test_lowest_sector_reproduces_the_unconstrained_ground_state():
    """min over sectors == the ordinary, sector-free ground state, and
    switching sector targeting back off restores the same answer."""
    n = 8
    fc, h, number = tV_chain(n)
    fc.set_hamiltonian(h)
    e_dense = fc.gs_energy()
    energies = []
    for nf in range(n + 1):
        fc.set_conserved_sector(Nf=nf)
        energies.append(fc.gs_energy())
    assert min(energies) == pytest.approx(e_dense, abs=DMRG_TOL)
    fc.set_conserved_sector()
    assert fc.gs_energy() == pytest.approx(e_dense, abs=DMRG_TOL)


@pytest.mark.parametrize("two_sz", [-4, 0, 2, 6])
def test_total_sz_sector_matches_ed(two_sz):
    """Total-Sz sectors of a Heisenberg chain, in ITensor's 2*Sz units.
    Passing at all requires the Sx/Sy -> S+/S- expansion and the
    cancellation of the S+S+/S-S- strings, since the surviving strings are
    checked for conserving the sector one by one."""
    n = 8
    sc, h, sz = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=two_sz)
    e = sc.gs_energy()
    assert e == pytest.approx(ed_sector_energy(sc, h, sz, two_sz / 2.0), abs=DMRG_TOL)
    assert sc.vev(sz).real == pytest.approx(two_sz / 2.0, abs=1e-8)


def test_spin_one_sector_matches_ed():
    """A spin-1 chain: three local states, Sz in units where the site
    carries +2/0/-2, so the charge table's units matter as well as its
    order."""
    n = 6
    sc = spinchain.Spin_Chain([3] * n)
    sc.setup_python()
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.maxm, sc.nsweeps = 40, 10
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=2)
    e = sc.gs_energy()
    assert e == pytest.approx(ed_sector_energy(sc, h, sum(sc.Sz), 1.0), abs=DMRG_TOL)
    assert sc.vev(sum(sc.Sz)).real == pytest.approx(1.0, abs=1e-8)


def test_boson_number_sector_matches_ed():
    """Boson sites are dmrgpy's own site type, and their charge is the
    occupation number itself rather than a +-1 grading."""
    n, nb = 4, 3
    bc = bosonchain.Bosonic_Chain(n, maxnb=[3] * n)
    bc.setup_python()
    h = 0
    for i in range(n - 1):
        h = h + bc.Adag[i] * bc.A[i + 1] + bc.Adag[i + 1] * bc.A[i]
    h = h + sum(1.5 * bc.N[i] * bc.N[i] for i in range(n))
    bc.maxm, bc.nsweeps = 40, 8
    bc.set_hamiltonian(h)
    bc.set_conserved_sector(Nb=nb)
    e = bc.gs_energy()
    assert e == pytest.approx(ed_sector_energy(bc, h, sum(bc.N), nb), abs=DMRG_TOL)
    assert bc.vev(sum(bc.N)).real == pytest.approx(nb, abs=1e-8)


def test_hubbard_conserves_a_chosen_subset_of_its_quantum_numbers():
    """A native Hubbard chain offers Nf and Sz independently: asking for
    Nf alone must leave Sz free (and an Sz-breaking term legal), while
    asking for both must pin both."""
    n = 4
    hc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    hc.setup_python()
    h = 0
    for i in range(n - 1):
        for c, cdag in [(hc.Cup, hc.Cdagup), (hc.Cdn, hc.Cdagdn)]:
            h = h + cdag[i] * c[i + 1] + cdag[i + 1] * c[i]
    h = h + sum(4.0 * hc.Nup[i] * hc.Ndn[i] for i in range(n))
    hflip = h + sum(0.3 * (hc.Cdagup[i] * hc.Cdn[i] + hc.Cdagdn[i] * hc.Cup[i])
                    for i in range(n))
    hc.maxm, hc.nsweeps = 40, 8

    hc.set_hamiltonian(hflip)
    hc.set_conserved_sector(Nf=n)
    hc.gs_energy()
    assert hc.vev(sum(hc.Ntot)).real == pytest.approx(n, abs=1e-8)

    # the same Sz-breaking Hamiltonian is not buildable once Sz is pinned
    hc.set_conserved_sector(Nf=n, Sz=0)
    with pytest.raises(ValueError):
        hc.gs_energy()

    # ...while the Sz-conserving one is, and matches ED in that sector
    hc.set_hamiltonian(h)
    hc.set_conserved_sector(Nf=n, Sz=0)
    e = hc.gs_energy()
    ed = hc.get_ED_obj()
    H = np.array(MO2matrix(h, ed).todense())
    N = np.array(MO2matrix(sum(hc.Ntot), ed).todense())
    S = np.array(MO2matrix(sum(hc.Nup[i] - hc.Ndn[i] for i in range(n)), ed).todense())
    sel = (np.abs(np.diag(N).real - n) < 1e-9) & (np.abs(np.diag(S).real) < 1e-9)
    assert e == pytest.approx(float(np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0]),
                              abs=DMRG_TOL)


def test_excited_states_stay_inside_the_sector():
    """The overlap-penalty excited states are confined to the sector too.

    What is asserted is confinement, not spectral accuracy. This backend's
    excited-state solver is variational with an overlap penalty and is not
    accurate to DMRG_TOL on this model even *without* a sector (measured:
    -1.977 against an exact -2.002 on the same chain, dense) -- that is a
    pre-existing property of the solver, not of sector mode. What sector
    mode has to guarantee is that every state it returns lives in the
    requested sector, which is checked two ways: each state's own <Sz>,
    and the fact that no energy falls below the sector's ED minimum. The
    second is the discriminating half -- Sz=2 is deliberately *not* the
    global ground-state sector, so an unconfined solver would find states
    below that bound (the global ground state is 0.49 lower).
    """
    n, two_sz = 6, 2
    sc, h, sz = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=two_sz)
    energies, wfs = sc.get_excited_states(n=3)
    energies = np.array(energies).real

    e_sector = ed_sector_energy(sc, h, sz, two_sz / 2.0)
    assert min(energies) == pytest.approx(e_sector, abs=DMRG_TOL)
    assert np.all(energies >= e_sector - DMRG_TOL)
    for wf in wfs:
        assert sc.vev(sz, wf=wf).real == pytest.approx(two_sz / 2.0, abs=1e-8)

    # ...and the bound really does exclude states an unconfined solver
    # would reach: the global ground state sits well below it
    sc.set_conserved_sector()
    assert sc.gs_energy() < e_sector - 0.1


def test_conserving_observables_work_and_others_raise_by_name():
    """Everything built on a sector-mode chain is charge-checked first: a
    conserving observable just works, a non-conserving one raises a
    ValueError naming the offending term."""
    n = 6
    sc, h, sz = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=0)
    sc.gs_energy()

    assert sc.vev(sc.Sz[0] * sc.Sz[1]).real == pytest.approx(
        sc.vev(sc.Sz[4] * sc.Sz[5]).real, abs=1e-6)  # translation invariance
    with pytest.raises(ValueError, match="Sz=2|Sz=-2|definite charge"):
        sc.vev(sc.Sx[0])
    # a Sz-conserving dynamical correlator is fine on a sector chain
    _, spectrum = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
                                              es=np.linspace(-1.0, 3.0, 20))
    assert np.all(np.isfinite(spectrum))
    # ...an Sx one is not, and says so instead of answering
    with pytest.raises(ValueError):
        sc.get_dynamical_correlator(name=(sc.Sx[0], sc.Sx[0]),
                                    es=np.linspace(-1.0, 3.0, 20))


def test_non_conserving_hamiltonian_is_rejected():
    """A transverse field breaks Sz; no expansion can rescue it."""
    n = 6
    sc = spinchain.Spin_Chain([2] * n)
    sc.setup_python()
    h = sum(sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1)) \
        + sum(0.3 * sc.Sx[i] for i in range(n))
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=0)
    with pytest.raises(ValueError):
        sc.gs_energy()


def test_a_rejected_request_leaves_the_chain_untouched():
    """Switching sector rebuilds the site set, so a request that turns out
    to be impossible must roll back rather than leave the chain half
    switched into a sector whose every later call would fail."""
    n = 8
    fc, h, _ = tV_chain(n)
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=3)
    e3 = fc.gs_energy()
    with pytest.raises(ValueError):
        fc.set_conserved_sector(Nf=99)
    assert fc.conserved_sector == {"Nf": 3}
    assert fc.gs_energy() == pytest.approx(e3, abs=DMRG_TOL)


def test_invalid_sector_requests_raise():
    n = 6
    fc, h, _ = tV_chain(n)
    fc.set_hamiltonian(h)
    with pytest.raises(ValueError, match="empty"):
        fc.set_conserved_sector(Nf=99)  # unreachable target
    with pytest.raises(ValueError):
        fc.set_conserved_sector(Sz=0)  # spinless fermions carry no Sz
    with pytest.raises(ValueError):
        fc.set_conserved_sector(Nq=1)  # not a quantum number at all
    with pytest.raises(ValueError):
        fc.set_conserved_sector(Nf=1.5)  # not an integer


def test_ed_and_other_backends_refuse_rather_than_answering_globally():
    """A sector-mode chain must never silently fall back to a solver with
    no quantum numbers: that would return the global ground state."""
    n = 6
    fc, h, _ = tV_chain(n)
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=3)
    with pytest.raises(NotImplementedError):
        fc.gs_energy(mode="ED")
    with pytest.raises(NotImplementedError):
        fc.setup_cpp(version=2)


def test_metts_refuses_instead_of_averaging_over_every_sector():
    """METTS samples product states drawn from every sector at once, so a
    sector cannot restrict its ensemble -- it must refuse rather than
    return a number that looks like a sector average and isn't."""
    n = 6
    sc, h, _ = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=0)
    with pytest.raises(ValueError, match="conserved-sector mode"):
        sc.metts_vev(sc.Sz[0], T=1.0, nsamples=1, nwarmup=1)
