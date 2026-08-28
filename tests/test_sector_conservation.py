"""Conserved-sector (quantum-number) targeting: Many_Body_Chain.
set_conserved_sector(), backed by mpscpp3's QN-carrying site sets.

The reference in every energy test here is ED restricted to the sector
*exactly*, not ED filtered by an expectation value: the conserved
operator is diagonal in the ED product basis, so a sector is a set of
basis states and the restriction is a submatrix (see ed_sector_energy
below). Filtering full-spectrum eigenvectors by <N> instead looks
equivalent and is not -- degenerate eigenvalues shared across sectors
come back as arbitrary superpositions, and numpy hands out mixtures whose
<N> matches no sector at all. That mistake produced a 0.48 "discrepancy"
against a DMRG answer that was in fact correct.

These tests all run on itensor_version=3, dmrgpy's block-sparse QN
implementation of sector mode. itensor_version="python" implements the
same API by a different mechanism (dense storage plus a charge penalty,
see pyitensor/sector.py) and is covered by
test_sector_conservation_python.py; the backends that have no quantum
numbers at all still refuse rather than silently answering with the
global ground state, which is what
test_ed_and_other_backends_refuse_rather_than_answering_globally checks.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain
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
    """Spinless t-V chain plus its total-number operator."""
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + v * fc.N[i] * fc.N[i + 1]
    return fc, h, sum(fc.N)


def heisenberg_chain(n):
    """Heisenberg chain written the way dmrgpy's examples write it --
    Sx.Sx+Sy.Sy+Sz.Sz -- which is exactly the case where no individual
    term conserves Sz and only their sum does (see mo_terms.h)."""
    sc = spinchain.Spin_Chain([2] * n)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    return sc, h, sum(sc.Sz)


@pytest.mark.parametrize("nf", [0, 1, 3, 4, 6, 8])
def test_fermion_number_sector_matches_ed(nf):
    """Every particle-number sector of an 8-site t-V chain."""
    n = 8
    fc, h, number = tV_chain(n)
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=nf)
    e = fc.gs_energy()
    assert e == pytest.approx(ed_sector_energy(fc, h, number, nf), abs=DMRG_TOL)
    # the state really is in that sector, not merely close in energy
    assert fc.vev(number).real == pytest.approx(nf, abs=1e-8)


def test_lowest_sector_reproduces_the_unconstrained_ground_state():
    """min over sectors == the ordinary, sector-free ground state, and the
    sector containing it returns exactly the same energy."""
    n = 8
    fc, h, number = tV_chain(n)
    fc.set_hamiltonian(h)
    e_dense = fc.gs_energy()
    energies = []
    for nf in range(n + 1):
        fc.set_conserved_sector(Nf=nf)
        energies.append(fc.gs_energy())
    assert min(energies) == pytest.approx(e_dense, abs=DMRG_TOL)
    # ...and switching sector targeting back off restores the same answer
    fc.set_conserved_sector()
    assert fc.gs_energy() == pytest.approx(e_dense, abs=DMRG_TOL)


@pytest.mark.parametrize("two_sz", [-6, -2, 0, 2, 8])
def test_total_sz_sector_matches_ed(two_sz):
    """Total-Sz sectors of a Heisenberg chain, in ITensor's 2*Sz units.
    Passing at all requires the Sx/Sy -> S+/S- expansion and the
    cancellation of the S+S+/S-S- strings, since AutoMPO aborts on a
    flux-violating term even when its coefficient cancels exactly."""
    n = 8
    sc, h, sz = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=two_sz)
    e = sc.gs_energy()
    assert e == pytest.approx(ed_sector_energy(sc, h, sz, two_sz / 2.0), abs=DMRG_TOL)
    assert sc.vev(sz).real == pytest.approx(two_sz / 2.0, abs=1e-8)


def test_hubbard_conserves_a_chosen_subset_of_its_quantum_numbers():
    """A native Hubbard chain offers Nf and Sz independently: asking for
    Nf alone must leave Sz free (and an Sz-breaking term legal), while
    asking for both must pin both."""
    n = 4
    hc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = 0
    for i in range(n - 1):
        for c, cdag in [(hc.Cup, hc.Cdagup), (hc.Cdn, hc.Cdagdn)]:
            h = h + cdag[i] * c[i + 1] + cdag[i + 1] * c[i]
    h = h + sum(4.0 * hc.Nup[i] * hc.Ndn[i] for i in range(n))
    # a spin-flip term: conserves Nf, breaks Sz
    hflip = h + sum(0.3 * (hc.Cdagup[i] * hc.Cdn[i] + hc.Cdagdn[i] * hc.Cup[i])
                    for i in range(n))

    hc.set_hamiltonian(hflip)
    hc.set_conserved_sector(Nf=n)
    e_nf = hc.gs_energy()
    assert hc.vev(sum(hc.Ntot)).real == pytest.approx(n, abs=1e-8)

    # the same Sz-breaking Hamiltonian is not buildable once Sz is pinned
    hc.set_conserved_sector(Nf=n, Sz=0)
    with pytest.raises(ValueError):
        hc.gs_energy()

    # ...while the Sz-conserving Hamiltonian is, and sits above the
    # Sz-free answer (a smaller variational space, same particle number)
    hc.set_hamiltonian(h)
    hc.set_conserved_sector(Nf=n, Sz=0)
    e_nf_sz = hc.gs_energy()
    hc.set_conserved_sector(Nf=n)
    assert e_nf_sz >= hc.gs_energy() - DMRG_TOL
    assert np.isfinite(e_nf)  # the Nf-only run above produced a real energy


def test_excited_states_stay_inside_the_sector():
    """The overlap-penalty excited states are confined to the sector too,
    so they reproduce the sector's own low-lying ED spectrum."""
    n = 8
    sc, h, sz = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=0)
    got = np.sort(np.array(sc.get_excited(n=3)).real)

    ed = sc.get_ED_obj()
    H = np.array(MO2matrix(h, ed).todense())
    S = np.array(MO2matrix(sz, ed).todense())
    sel = np.abs(np.diag(S).real) < 1e-9
    want = np.linalg.eigvalsh(H[np.ix_(sel, sel)])[:3]
    assert got == pytest.approx(want, abs=1e-5)


def test_conserving_observables_work_and_others_raise_by_name():
    """Everything built on a sector-mode chain is flux-checked first: a
    conserving observable just works, a non-conserving one raises a
    ValueError naming it -- rather than aborting the interpreter, which
    is what ITensor does on its own for both of these.
    """
    n = 8
    sc, h, sz = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=0)
    sc.gs_energy()

    assert sc.vev(sc.Sz[0] * sc.Sz[1]).real == pytest.approx(
        sc.vev(sc.Sz[6] * sc.Sz[7]).real, abs=1e-6)  # translation invariance
    with pytest.raises(ValueError, match="Sz=2|Sz=-2|definite charge"):
        sc.vev(sc.Sx[0])
    # a Sz-conserving dynamical correlator is fine on a sector chain
    _, spectrum = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
                                              es=np.linspace(-1.0, 3.0, 20))
    assert np.all(np.isfinite(spectrum))
    # ...an Sx one is not, and says so instead of crashing
    with pytest.raises(ValueError):
        sc.get_dynamical_correlator(name=(sc.Sx[0], sc.Sx[0]),
                                    es=np.linspace(-1.0, 3.0, 20))


def test_non_conserving_hamiltonian_is_rejected():
    """A transverse field breaks Sz; no expansion can rescue it."""
    n = 6
    sc = spinchain.Spin_Chain([2] * n)
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
    fc2, h2, _ = tV_chain(n)
    fc2.setup_cpp(version=2)
    with pytest.raises(NotImplementedError):
        fc2.set_conserved_sector(Nf=3)


def test_switching_to_the_python_backend_keeps_the_sector():
    """itensor_version="python" does have quantum numbers (it grades its
    site indices and penalizes charge excursions instead of storing blocks
    -- see pyitensor/sector.py), so switching a sector-mode chain to it is
    allowed, and must carry the sector across rather than silently
    dropping it."""
    n, nf = 8, 3
    fc, h, number = tV_chain(n)
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=nf)
    e_v3 = fc.gs_energy()
    fc.setup_python()
    assert fc.conserved_sector == {"Nf": nf}
    assert fc.gs_energy() == pytest.approx(e_v3, abs=DMRG_TOL)
    assert fc.vev(number).real == pytest.approx(nf, abs=1e-8)


def test_dense_only_algorithms_refuse_instead_of_aborting():
    """METTS and TEBD assemble tensors by hand (dense links, bond gates
    summed before their fluxes agree). Under QN indices ITensor does not
    fail cleanly there -- it calls abort() -- so these paths must refuse
    up front. TDVP, the default time-evolution method, does work."""
    n = 6
    sc, h, _ = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=0)
    sc.gs_energy()

    with pytest.raises(ValueError, match="conserved-sector mode"):
        sc.metts_vev(sc.Sz[0], T=1.0, nsamples=1, nwarmup=1)

    sc.tevol_method = "TEBD"
    with pytest.raises(ValueError, match="conserved-sector mode"):
        sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), submode="TD",
                                    es=np.linspace(-1.0, 3.0, 8))

    sc.tevol_method = "TDVP"
    _, spectrum = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
                                              submode="TD",
                                              es=np.linspace(-1.0, 3.0, 8))
    assert np.all(np.isfinite(spectrum))


def test_sector_survives_a_clone():
    """clone()/deepcopy build a fresh session; it has to be told about the
    sector again or the clone would search the whole Hilbert space."""
    n = 8
    fc, h, number = tV_chain(n)
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=5)
    e = fc.gs_energy()
    clone = fc.clone()
    assert clone.conserved_sector == {"Nf": 5}
    assert clone.gs_energy() == pytest.approx(e, abs=DMRG_TOL)


def test_changing_sector_between_calls_actually_re_solves():
    """The session caches its energy across a skipped Hamiltonian re-send,
    so the sector has to be part of groundstate.py's send-cache key."""
    n = 8
    fc, h, _ = tV_chain(n)
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=3)
    e3 = fc.gs_energy()
    fc.set_conserved_sector(Nf=6)
    e6 = fc.gs_energy()
    assert abs(e3 - e6) > 1.0
