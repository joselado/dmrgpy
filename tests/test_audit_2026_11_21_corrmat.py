"""Regressions of docs/audit_2026_11_21_hole_hunt.md, findings 8 and 9: the
default route of get_correlation_matrix in a conserved sector.

In a sector the default was dmmode="full", the one-sweep matrix
<Cdag_i C_j> of the chain's own C, whatever operators= asked for: a
Hadamard rotation of two modes came back 0.507 off, C[:2] came back as the
whole (n, n) matrix, and the correlation entropy of a subsystem was that of
the whole pure chain, 0. And "full" ran with no site-type test, so on a
native spinful chain in an Nf sector the default call raised, where the
same call answers without a sector.

The fix is one predicate, correlationentropy._full_applies, read by both
branches of _default_dmmode: "full" only for no operators= (or exactly the
chain's own C, in order) on spinless fermion sites with a session or a
julia_live state that has the one-sweep matrix; otherwise "explicit" in a
sector and "fast" outside one. An explicit dmmode="full" with other
operators raises.

The anchors are ED restricted to the sector (a separate chain pinned to
mode="ED"), the exact rotation <X_a^dag X_b> = (conj(U) M0 U^T)_ab of the
chain's own matrix M0 for X_a = sum_j U_aj C_j, and, for the free chain,
the Slater determinant of the filled single-particle orbitals."""

import numpy as np
import pytest

from dmrgpy import cppext, fermionchain, spinchain
from dmrgpy.entropytk import correlationentropy as ce

from _helpers import julia_live_param, setup_backend

V3 = pytest.param(3, marks=pytest.mark.skipif(
    not cppext.available(3),
    reason="requires the compiled mpscpp3 (ITensor v3) extension"))
VERSIONS = ["python", V3]
VERSIONS_JL = ["python", V3, julia_live_param()]

R = float(2 ** -0.5)
HADAMARD = np.array([[R, R, 0, 0], [R, -R, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])


def _backend(chain, version):
    """Put `chain` on `version`, or pin it to ED for version "ED"."""
    if version == "ED":
        chain.mode = "ED"
    else:
        setup_backend(chain, version)


def _four_site_chain(version):
    """Hopping chain with three different onsite energies, no symmetry."""
    np.random.seed(0)
    fc = fermionchain.Fermionic_Chain(4)
    _backend(fc, version)
    h = 0
    for i in range(3):
        h = h - 1.0 * (fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i])
    h = h + 0.3 * fc.N[0] - 0.2 * fc.N[3] + 0.1 * fc.N[1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps, fc.cutoff = 16, 10, 1e-12
    return fc


def _six_site_chain(version):
    """Interacting hopping chain with no reflection symmetry."""
    np.random.seed(3)
    fc = fermionchain.Fermionic_Chain(6)
    _backend(fc, version)
    h = 0
    for i in range(5):
        h = h + fc.Cdag[i] * fc.C[i + 1]
    h = h + h.get_dagger()
    for i in range(5):
        h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    h = h + 0.4 * fc.N[0] - 0.2 * fc.N[2]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps, fc.cutoff = 40, 10, 1e-12
    return fc


def _combine(fc, U):
    """X_a = sum_j U[a, j] C_j, one MultiOperator per row of U."""
    ops = []
    for row in U:
        o = None
        for j, u in enumerate(row):
            if abs(u) < 1e-14:
                continue
            c = float(np.real(u)) if abs(np.imag(u)) < 1e-15 else complex(u)
            o = c * fc.C[j] if o is None else o + c * fc.C[j]
        ops.append(o)
    return ops


def _peschel(ns):
    """Correlation entropy of a set of occupations, as entanglement.py."""
    ns = np.clip(np.real(ns), 1e-12, 1. - 1e-12)
    return float(-np.sum(ns * np.log(ns) + (1. - ns) * np.log(1. - ns)))


def _mx(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.shape == b.shape, (a.shape, b.shape)
    return float(np.max(np.abs(a - b)))


@pytest.mark.parametrize("version", VERSIONS_JL)
def test_sector_default_reads_a_rotated_basis_and_a_subset(version):
    """Nf=2 on 4 sites: the default reads a Hadamard rotation of modes 0
    and 1 and a subset, as dmmode="explicit" and ED do."""
    fc = _four_site_chain(version)
    fc.set_conserved_sector(Nf=2)
    ed = _four_site_chain("ED")
    ed.set_conserved_sector(Nf=2)
    assert fc.get_mode() == "DMRG"
    wf = fc.get_gs()
    ops = _combine(fc, HADAMARD)
    # the resolved route: the one sweep only for the chain's own C
    assert ce._default_dmmode(fc, wf, "electron") == "full"
    assert ce._default_dmmode(fc, wf, "electron", list(fc.C)) == "full"
    assert ce._default_dmmode(fc, wf, "electron", ops) == "explicit"
    assert ce._default_dmmode(fc, wf, "electron", fc.C[:2]) == "explicit"
    assert ce._default_dmmode(fc, wf, "electron", fc.C[::-1]) == "explicit"
    M0 = np.asarray(fc.get_correlation_matrix(wf=wf))
    got = np.asarray(fc.get_correlation_matrix(wf=wf, operators=ops))
    explicit = np.asarray(fc.get_correlation_matrix(wf=wf, operators=ops,
                                                    dmmode="explicit"))
    ref = np.asarray(ed.get_correlation_matrix(operators=_combine(ed, HADAMARD)))
    rotated = np.conj(HADAMARD) @ M0 @ HADAMARD.T
    # the rotation moves this matrix, so returning M0 cannot pass
    assert _mx(rotated, M0) > 0.1
    assert _mx(got, explicit) < 1e-10
    assert _mx(got, rotated) < 1e-10
    assert _mx(got, ref) < 1e-6
    # a subset is the block of those modes, not the whole matrix
    sub = np.asarray(fc.get_correlation_matrix(wf=wf, operators=fc.C[:2]))
    assert sub.shape == (2, 2)
    assert _mx(sub, M0[:2, :2]) < 1e-10
    # the chain's own C, passed explicitly, is the default matrix
    own = np.asarray(fc.get_correlation_matrix(wf=wf, operators=fc.C))
    assert _mx(own, M0) < 1e-12


@pytest.mark.parametrize("version", VERSIONS)
def test_sector_default_reads_a_momentum_basis_and_a_reversed_list(version):
    """Nf=3 on 6 interacting sites: a momentum basis (complex
    coefficients) and the reversed site list, against ED in the sector."""
    n = 6
    fc = _six_site_chain(version)
    fc.set_conserved_sector(Nf=3)
    ed = _six_site_chain("ED")
    ed.set_conserved_sector(Nf=3)
    assert fc.get_mode() == "DMRG"
    wf = fc.get_gs()
    M0 = np.asarray(fc.get_correlation_matrix(wf=wf))
    U = np.exp(2j * np.pi * np.outer(np.arange(n), np.arange(n)) / n) / np.sqrt(n)
    got = np.asarray(fc.get_correlation_matrix(wf=wf, operators=_combine(fc, U)))
    ref = np.asarray(ed.get_correlation_matrix(operators=_combine(ed, U)))
    assert _mx(np.conj(U) @ M0 @ U.T, M0) > 0.1
    assert _mx(got, ref) < 1e-6
    assert _mx(got, np.conj(U) @ M0 @ U.T) < 1e-10
    rev = np.asarray(fc.get_correlation_matrix(wf=wf, operators=fc.C[::-1]))
    ref = np.asarray(ed.get_correlation_matrix(operators=ed.C[::-1]))
    assert _mx(M0[::-1, ::-1], M0) > 1e-2
    assert _mx(rev, ref) < 1e-6
    assert _mx(rev, M0[::-1, ::-1]) < 1e-10


@pytest.mark.parametrize("version", VERSIONS)
def test_sector_correlation_entropy_of_a_subsystem(version):
    """Free fermions on 6 sites at Nf=3: the ground state is one Slater
    determinant, so the whole chain has zero correlation entropy and the
    first three modes have the Peschel entropy of their block of the
    filled-orbital matrix."""
    n, nf, nsub = 6, 3, 3
    t = np.zeros((n, n))
    for i in range(n - 1):
        t[i, i + 1] = t[i + 1, i] = 1.0
    t[0, 0], t[2, 2] = 0.4, -0.2
    np.random.seed(1)
    fc = fermionchain.Fermionic_Chain(n)
    setup_backend(fc, version)
    # the same H, sum_ij t_ij Cdag_i C_j, written as the record writes it
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1]
    h = h + h.get_dagger()
    h = h + 0.4 * fc.N[0] - 0.2 * fc.N[2]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps, fc.cutoff = 40, 10, 1e-12
    fc.set_conserved_sector(Nf=nf)
    assert fc.get_mode() == "DMRG"
    # <Cdag_i C_j> = sum over the nf lowest orbitals of conj(v_i) v_j
    v = np.linalg.eigh(t)[1][:, :nf]
    exact = np.conj(v) @ v.T
    s_exact = _peschel(np.linalg.eigvalsh(exact[:nsub, :nsub]))
    assert s_exact > 0.5
    wf = fc.get_gs()
    assert ce._default_dmmode(fc, wf, "electron", fc.C[:nsub]) == "explicit"
    s_sub = fc.get_correlation_entropy(operators=fc.C[:nsub])
    assert s_sub == pytest.approx(s_exact, abs=1e-6)
    assert fc.get_correlation_entropy() == pytest.approx(0., abs=1e-4)


def _native_chain(version):
    """Hubbard chain on native spinful sites, no spin or site symmetry."""
    np.random.seed(4)
    fc = fermionchain.Spinful_Fermionic_Chain_Native(4)
    _backend(fc, version)
    h = 0
    for i in range(3):
        h = h + fc.Cdagup[i] * fc.Cup[i + 1] + fc.Cdagdn[i] * fc.Cdn[i + 1]
    h = h + h.get_dagger()
    h = h + 2.0 * sum(fc.Nup[i] * fc.Ndn[i] for i in range(4))
    h = h + 0.3 * fc.Nup[0] - 0.1 * fc.Ndn[2]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps, fc.cutoff = 40, 10, 1e-12
    return fc


@pytest.mark.parametrize("version", VERSIONS)
def test_native_spinful_chain_in_a_sector_answers_the_default_call(version):
    """Nf=4 on 4 native spinful sites: get_correlation_matrix() and
    get_correlation_entropy() with no keywords answer, as ED in the
    sector does, where "full" raised on these sites."""
    fc = _native_chain(version)
    fc.set_conserved_sector(Nf=4)
    ed = _native_chain("ED")
    ed.set_conserved_sector(Nf=4)
    assert fc.get_mode() == "DMRG"
    wf = fc.get_gs()
    assert ce._default_dmmode(fc, wf, "electron") == "explicit"
    ref = np.asarray(ed.get_correlation_matrix())
    assert np.trace(ref).real == pytest.approx(4., abs=1e-10)
    got = np.asarray(fc.get_correlation_matrix())
    assert got.shape == (8, 8)
    assert _mx(got, ref) < 1e-6
    s = fc.get_correlation_entropy()
    assert s == pytest.approx(_peschel(np.linalg.eigvalsh(ref)), abs=1e-6)


@pytest.mark.parametrize("version", VERSIONS)
def test_spin_chain_in_a_sector_answers_with_operators(version):
    """An XXZ chain at Sz=0 with operators=[Sz_i], which conserve Sz:
    the default answers <Sz_i Sz_j>, as ED in the sector."""
    def spin_chain(v):
        np.random.seed(2)
        sc = spinchain.Spin_Chain(["S=1/2"] * 6)
        _backend(sc, v)
        h = 0
        for i in range(5):
            h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
                + 0.7 * sc.Sz[i] * sc.Sz[i + 1]
        h = h + 0.3 * sc.Sz[0]
        sc.set_hamiltonian(h)
        sc.maxm, sc.nsweeps, sc.cutoff = 40, 10, 1e-12
        sc.set_conserved_sector(Sz=0)
        return sc
    sc = spin_chain(version)
    ed = spin_chain("ED")
    assert sc.get_mode() == "DMRG"
    got = np.asarray(sc.get_correlation_matrix(operators=sc.Sz))
    ref = np.asarray(ed.get_correlation_matrix(operators=ed.Sz))
    assert got.shape == (6, 6)
    assert _mx(got, ref) < 1e-6


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("nf", [None, 2])
def test_explicit_full_with_other_operators_raises_before_the_solve(version, nf):
    """dmmode="full" never reads operators=, so for operators that are not
    the chain's own C it raises, in or out of a sector, before any ground
    state is solved; for the chain's own C, rebuilt, it still answers."""
    fc = _four_site_chain(version)
    if nf is not None:
        fc.set_conserved_sector(Nf=nf)
    assert getattr(fc, "computed_gs", False) is False
    for other in (_combine(fc, HADAMARD), fc.C[:2], fc.C[::-1],
                  [2.0 * c for c in fc.C]):
        with pytest.raises(ValueError, match="explicit"):
            fc.get_correlation_matrix(operators=other, dmmode="full")
    assert getattr(fc, "computed_gs", False) is False
    own = [fc.get_operator("C", i) for i in range(4)]
    wf = fc.get_gs()
    full = np.asarray(fc.get_correlation_matrix(wf=wf, operators=own,
                                                dmmode="full"))
    explicit = np.asarray(fc.get_correlation_matrix(wf=wf, operators=own,
                                                    dmmode="explicit"))
    assert _mx(full, explicit) < 1e-10


@pytest.mark.parametrize("version", VERSIONS)
def test_no_sector_default_is_unchanged(version):
    """Without a sector the default takes the routes it took before: "full"
    for the chain's own C, "fast" for other operators and for Nambu."""
    fc = _four_site_chain(version)
    assert fc.get_mode() == "DMRG"
    wf = fc.get_gs()
    ops = _combine(fc, HADAMARD)
    assert ce._default_dmmode(fc, wf, "electron") == "full"
    assert ce._default_dmmode(fc, wf, "electron", ops) == "fast"
    assert ce._default_dmmode(fc, wf, "electron", fc.C[:2]) == "fast"
    assert ce._default_dmmode(fc, wf, "Nambu") == "fast"
    M0 = np.asarray(fc.get_correlation_matrix(wf=wf))
    assert _mx(M0, fc.get_correlation_matrix(wf=wf, dmmode="full")) < 1e-10
    got = np.asarray(fc.get_correlation_matrix(wf=wf, operators=ops))
    assert _mx(got, fc.get_correlation_matrix(wf=wf, operators=ops,
                                              dmmode="fast")) < 1e-10
    assert _mx(got, np.conj(HADAMARD) @ M0 @ HADAMARD.T) < 1e-10
    sub = np.asarray(fc.get_correlation_matrix(wf=wf, operators=fc.C[:2]))
    assert _mx(sub, M0[:2, :2]) < 1e-10
