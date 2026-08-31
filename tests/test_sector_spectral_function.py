"""Sector-resolved dynamical correlators and spectral functions:
`get_dynamical_correlator(submode="SECTOR")`, `get_spectral_function` and
`get_spin_spectral_function` (src/dmrgpy/sectordc.py).

The method sums a Lehmann representation over the lowest eigenstates of
the one quantum-number sector `B|gs>` lands in -- N+1 for `Cdag`, Sz+2 for
`S+` -- instead of over excited states of the whole Hilbert space. So the
references here are per-sector: ED *restricted* to a sector (a submatrix
in the ED product basis), never ED filtered by an expectation value, for
the reason test_sector_conservation.py spells out -- sector-degenerate
eigenvalues come back from a full diagonalization as arbitrary
superpositions.

Two of these tests exist for specific, easy-to-get-wrong reasons rather
than for coverage:

  * test_ladder_pair_matches_ed uses (S+,S-), a *non*-self-adjoint
    operator pair. A misplaced conjugation in the matrix elements is
    invisible for a pair like (C,Cdag) and shows up immediately here.
  * test_spectral_function_obeys_the_anticommutator_sum_rule checks
    sum_n w_n = <0|{c_i,c_i^dag}|0> = 1 across the particle *and* hole
    parts together -- an exact identity that no single sector can satisfy
    on its own, so it tests the assembly rather than either half.

Chains are kept small and the ground state non-degenerate on purpose: the
ED reference `dynamical_correlator_ED` averages over a near-degenerate
ground manifold with equal weights below its `dex` cutoff, while this
method uses one DMRG ground state, so on a degenerate chain the two would
disagree by construction rather than by bug.
"""
import warnings

import numpy as np
import pytest

from dmrgpy import fermionchain, sectordc, spinchain
from dmrgpy.multioperator import MO2matrix
from dmrgpy.multioperatortk import charge as chargetk
from dmrgpy.operatornames import name2MO

TOL = 1e-5


# -- models --------------------------------------------------------------

def tV_chain(n=6, v=1.5, backend="python"):
    """Spinless t-V chain with a site-0 potential, which lifts the
    accidental degeneracies a uniform chain has."""
    fc = fermionchain.Fermionic_Chain(n)
    if backend == "python":
        fc.setup_python()
    else:
        fc.setup_cpp(version=backend)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + v * fc.N[i] * fc.N[i + 1]
    h = h + 0.7 * fc.N[0]
    fc.set_hamiltonian(h)
    fc.maxm = 40
    fc.nsweeps = 8
    return fc, h


def heisenberg_chain(n=6, backend="python"):
    sc = spinchain.Spin_Chain([2] * n)
    if backend == "python":
        sc.setup_python()
    else:
        sc.setup_cpp(version=backend)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
            + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    sc.maxm = 40
    sc.nsweeps = 8
    return sc, h


def hubbard_chain(n=4, u=8.0):
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = 0
    for i in range(n - 1):
        for cd, c in ((fc.Cdagup, fc.Cup), (fc.Cdagdn, fc.Cdn)):
            h = h - cd[i] * c[i + 1] - cd[i + 1] * c[i]
    for i in range(n):
        h = h + u * fc.Nup[i] * fc.Ndn[i]
    fc.set_hamiltonian(h)
    fc.maxm = 60
    fc.nsweeps = 10
    return fc, h


# -- ED references -------------------------------------------------------

def ed_sector_lehmann(chain, h, charges, reference, target, A, B):
    """(energies, weights) of the exact Lehmann sum between two sectors:
    the ground state of `reference` and every eigenstate of `target`,
    both obtained by diagonalizing H restricted to that sector.

    `charges` maps a quantum-number name to its (diagonal) total-charge
    operator; `reference`/`target` map the same names to their targets.
    """
    ed = chain.get_ED_obj()
    H = np.array(MO2matrix(h, ed).todense())
    diag = {}
    for name, op in charges.items():
        Q = np.array(MO2matrix(op, ed).todense())
        assert np.allclose(Q, np.diag(np.diag(Q))), "charge is not diagonal"
        diag[name] = np.diag(Q).real

    def select(sector):
        sel = np.ones(H.shape[0], dtype=bool)
        for name, target_value in sector.items():
            sel &= np.abs(diag[name] - target_value) < 1e-9
        assert sel.any(), "empty sector in the ED reference"
        return sel

    sel0, sel1 = select(reference), select(target)
    e0s, v0s = np.linalg.eigh(H[np.ix_(sel0, sel0)])
    e1s, v1s = np.linalg.eigh(H[np.ix_(sel1, sel1)])
    v0 = np.zeros(H.shape[0], dtype=complex)
    v0[sel0] = v0s[:, 0]
    Am = np.array(MO2matrix(A, ed).todense())
    Bm = np.array(MO2matrix(B, ed).todense())
    energies, weights = [], []
    for k in range(len(e1s)):
        vk = np.zeros(H.shape[0], dtype=complex)
        vk[sel1] = v1s[:, k]
        a = np.conjugate(v0) @ (Am @ vk)
        b = np.conjugate(vk) @ (Bm @ v0)
        energies.append(e1s[k] - e0s[0])
        weights.append(a * b)
    return np.array(energies), np.array(weights)


def compare_poles(got_e, got_w, ref_e, ref_w, tol=TOL):
    """Every ED pole with appreciable weight is reproduced, at the same
    energy and with the same weight."""
    keep = np.abs(ref_w) > 1e-7
    assert keep.any(), "the ED reference itself carries no weight"
    for e, w in zip(ref_e[keep], ref_w[keep]):
        near = np.abs(got_e - e) < 1e-4
        assert near.any(), "missing pole at E=%.6f" % e
        assert complex(np.sum(got_w[near])) == pytest.approx(complex(w), abs=tol)


# -- the core submode ----------------------------------------------------

def test_free_fermion_poles_match_sector_restricted_ed():
    """The whole pipeline -- promotion, cross-sector matrix elements, the
    Jordan-Wigner strings, the Lehmann assembly -- on a chain small enough
    for `nex` to cover the target sector completely, so the sum is exact
    rather than truncated."""
    n, nf = 6, 3
    fc, h = tV_chain(n, v=0.0)  # non-interacting: a pure JW/free-fermion check
    fc.set_conserved_sector(Nf=nf)
    e, w, info = sectordc.sector_poles(fc, name=(fc.C[0], fc.Cdag[0]),
                                       nex=100, quiet=True)
    assert info["target_sector"] == {"Nf": nf + 1}
    assert info["captured"] == pytest.approx(1.0, abs=1e-6)
    ref_e, ref_w = ed_sector_lehmann(
        fc, h, {"Nf": sum(fc.N)}, {"Nf": nf}, {"Nf": nf + 1},
        fc.C[0], fc.Cdag[0])
    compare_poles(e, w, ref_e, ref_w)


def test_ladder_pair_matches_ed():
    """(S+,S-) -- a non-self-adjoint pair, where a misplaced conjugation
    in the matrix elements would show up and (C,Cdag) would hide it."""
    n = 6
    sc, h = heisenberg_chain(n)
    Sp, Sm = name2MO("Sp", sc), name2MO("Sm", sc)
    es = np.linspace(-0.2, 4.0, 300)
    for A, B, target in ((Sm[0], Sp[0], 2), (Sp[0], Sm[0], -2)):
        x, y, info = sc.get_dynamical_correlator(
            submode="SECTOR", name=(A, B), nex=100, es=es, delta=0.1,
            return_poles=True, quiet=True)
        assert info["target_sector"] == {"Sz": target}
        ref_e, ref_w = ed_sector_lehmann(
            sc, h, {"Sz": 2 * sum(sc.Sz)}, {"Sz": 0}, {"Sz": target}, A, B)
        compare_poles(info["energies"], info["weights"], ref_e, ref_w)


@pytest.mark.parametrize("backend", ["python", 3])
def test_curve_matches_the_exact_ed_correlator(backend):
    """The broadened curve, in the convention every other submode returns
    -- the real part of what SECTOR gives is what mode="ED",
    submode="ED" gives."""
    n = 6
    sc, h = heisenberg_chain(n, backend=backend)
    Sp, Sm = name2MO("Sp", sc), name2MO("Sm", sc)
    es = np.linspace(-0.2, 4.0, 300)
    x, y = sc.get_dynamical_correlator(
        submode="SECTOR", name=(Sm[0], Sp[0]), nex=100, es=es, delta=0.1,
        quiet=True)
    ref = spinchain.Spin_Chain([2] * n)
    ref.set_hamiltonian(h)
    xe, ye = ref.get_dynamical_correlator(
        mode="ED", submode="ED", name=(Sm[0], Sp[0]), es=es, delta=0.1)
    assert np.max(np.abs(y.real - ye)) < 1e-3 * max(1.0, np.max(np.abs(ye)))


def test_nex_is_capped_to_the_sector_dimension():
    """Asking for more states than the target sector holds is clipped and
    said out loud, instead of sending the overlap-penalty search after
    states that do not exist."""
    n = 6
    sc, _ = heisenberg_chain(n)
    Sp, Sm = name2MO("Sp", sc), name2MO("Sm", sc)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        x, y, info = sc.get_dynamical_correlator(
            submode="SECTOR", name=(Sm[0], Sp[0]), nex=500,
            es=np.linspace(0, 4, 50), return_poles=True)
    assert info["nex"] == 15  # dim of the Sz=+2 sector of 6 spins-1/2
    assert any("reduced to 15" in str(c.message) for c in caught)


def test_truncating_the_sum_warns_with_the_captured_weight():
    n = 6
    sc, _ = heisenberg_chain(n)
    Sp, Sm = name2MO("Sp", sc), name2MO("Sm", sc)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        x, y, info = sc.get_dynamical_correlator(
            submode="SECTOR", name=(Sm[0], Sp[0]), nex=1,
            es=np.linspace(0, 4, 50), return_poles=True)
    assert info["captured"] < 0.9
    assert any("spectral weight" in str(c.message) for c in caught)


def test_charge_indefinite_operators_are_decomposed():
    """Sx raises *and* lowers Sz, so it has no single target sector. It is
    split into its two charge components -- Sx=(S+ + S-)/2 -- and the two
    channels that survive the pairing are summed, which is also the
    physically right decomposition. Checked against the exact ED
    correlator, since a wrong split would still produce a plausible
    spectrum."""
    n = 6
    sc, h = heisenberg_chain(n)
    es = np.linspace(-0.2, 4.0, 300)
    ref = spinchain.Spin_Chain([2] * n)
    ref.set_hamiltonian(h)
    for A, B in ((sc.Sx[0], sc.Sx[0]), (sc.Sy[0], sc.Sy[0])):
        x, y, info = sc.get_dynamical_correlator(
            submode="SECTOR", name=(A, B), nex=100, es=es, delta=0.1,
            return_poles=True, quiet=True)
        assert info["target_sector"] == [{"Sz": -2}, {"Sz": 2}]
        assert info["captured"] == pytest.approx(1.0, abs=1e-6)
        xe, ye = ref.get_dynamical_correlator(mode="ED", submode="ED",
                                              name=(A, B), es=es, delta=0.1)
        assert np.max(np.abs(y.real - ye)) < 1e-4


def test_charge_components_sum_back_to_the_operator():
    """The decomposition is exact, not a projection: the pieces add up to
    the original operator, and each piece's *measured* charge is the one
    it was filed under."""
    n = 6
    sc, _ = heisenberg_chain(n)
    ed = sc.get_ED_obj()
    for op in (sc.Sx[0], sc.Sy[1], sc.Sz[0], sc.Sx[0] * sc.Sx[2]):
        parts = chargetk.charge_components(sc, op, ("Sz",))
        M = np.array(MO2matrix(op, ed).todense())
        S = sum(np.array(MO2matrix(v, ed).todense()) for v in parts.values())
        assert np.max(np.abs(M - S)) < 1e-12
        for q, piece in parts.items():
            assert chargetk.multioperator_charge(sc, piece, ("Sz",)) == q


def test_the_single_sector_primitive_refuses_a_multi_channel_pair():
    """sector_poles describes ONE sector -- its info dict names it and
    carries that sector's matrix elements -- so a pair spanning two
    channels has to go through sector_lehmann instead of being silently
    reduced to one of them."""
    sc, _ = heisenberg_chain(6)
    with pytest.raises(ValueError, match="charge channels"):
        sectordc.sector_poles(sc, name=(sc.Sx[0], sc.Sx[0]), nex=4)


def test_charges_that_do_not_cancel_are_rejected():
    """<0|A|n><n|B|0> is then zero for every n. Returning a zero spectrum
    for what is almost always a typo'd operator order would be worse."""
    sc, _ = heisenberg_chain(6)
    Sp = name2MO("Sp", sc)
    with pytest.raises(ValueError, match="do not cancel"):
        sc.get_dynamical_correlator(submode="SECTOR",
                                    name=(Sp[0], Sp[0]), nex=4)


def test_unsupported_backends_are_rejected_not_silently_substituted():
    sc, _ = heisenberg_chain(6)
    Sp, Sm = name2MO("Sp", sc), name2MO("Sm", sc)
    with pytest.raises(NotImplementedError, match="DMRG-only"):
        sc.get_dynamical_correlator(mode="ED", submode="SECTOR",
                                    name=(Sm[0], Sp[0]), nex=4)
    v2 = sc.copy()
    v2.setup_cpp(version=2)
    with pytest.raises(NotImplementedError, match="conserved-sector"):
        v2.get_dynamical_correlator(submode="SECTOR",
                                    name=(Sm[0], Sp[0]), nex=4)


def test_the_callers_chain_is_left_alone():
    """Every sector switch happens on an internal clone: the pipeline ends
    outside sector mode, which must not become a side effect."""
    n = 6
    fc, _ = tV_chain(n)
    fc.set_conserved_sector(Nf=3)
    fc.get_dynamical_correlator(submode="SECTOR", name=(fc.C[0], fc.Cdag[0]),
                                nex=4, es=np.linspace(0, 4, 50), quiet=True)
    assert fc.conserved_sector == {"Nf": 3}
    fc2, _ = tV_chain(n)
    fc2.get_dynamical_correlator(submode="SECTOR", name=(fc2.C[0], fc2.Cdag[0]),
                                 nex=4, es=np.linspace(0, 4, 50), quiet=True)
    assert fc2.conserved_sector is None


def test_the_reference_sector_is_measured_when_it_is_not_set():
    """A chain with no sector set still works: the reference sector is
    measured on the unconstrained ground state and rounded."""
    n = 6
    fc, _ = tV_chain(n)
    x, y, info = fc.get_dynamical_correlator(
        submode="SECTOR", name=(fc.C[0], fc.Cdag[0]), nex=4,
        es=np.linspace(0, 4, 50), return_poles=True, quiet=True)
    nf = info["reference_sector"]["Nf"]
    assert 0 < nf < n
    assert info["target_sector"] == {"Nf": nf + 1}


def test_sector_solves_are_shared_across_site_pairs():
    """One set of sector solves serves the whole matrix of (i,j) pairs --
    the property that makes a site sweep (and so S(q,w)/A(k,w)) affordable
    with this method."""
    n = 6
    fc, _ = tV_chain(n)
    fc.set_conserved_sector(Nf=3)
    kwargs = dict(submode="SECTOR", nex=6, es=np.linspace(0, 4, 50),
                  quiet=True)
    fc.get_dynamical_correlator(name=(fc.C[0], fc.Cdag[0]), **kwargs)
    cached = dict(fc._sector_states_cache)
    for i in range(n):
        fc.get_dynamical_correlator(name=(fc.C[i], fc.Cdag[i]), **kwargs)
    assert set(fc._sector_states_cache) == set(cached)
    fc.set_hamiltonian(fc.hamiltonian)  # calls restart()
    assert not fc._sector_states_cache


# -- the spectral functions ---------------------------------------------

def test_spectral_function_obeys_the_anticommutator_sum_rule():
    """sum over particle *and* hole poles of <0|c_i|n><n|c_i^dag|0> +
    <0|c_i^dag|m><m|c_i|0> = <0|{c_i,c_i^dag}|0> = 1. Neither half can
    satisfy this alone, so it tests the assembly."""
    n, nf = 6, 3
    fc, _ = tV_chain(n)
    fc.set_conserved_sector(Nf=nf)
    x, A, info = fc.get_spectral_function(
        i=1, nex=100, es=np.linspace(-6, 6, 200), delta=0.1,
        return_poles=True, quiet=True)
    assert float(np.sum(info["weights"]).real) == pytest.approx(1.0, abs=1e-6)
    assert float(np.sum(np.abs(info["weights"].imag))) < 1e-8


def test_spectral_function_gap_and_mu_match_ed():
    """mu and the charge gap come out of two independently converged
    sector energies, so they check the *difference* of large numbers
    rather than either energy on its own."""
    n, nf = 6, 3
    fc, h = tV_chain(n)
    fc.set_conserved_sector(Nf=nf)
    x, A, info = fc.get_spectral_function(
        i=1, nex=6, es=np.linspace(-6, 6, 100), delta=0.1,
        return_poles=True, quiet=True)
    ed = fc.get_ED_obj()
    H = np.array(MO2matrix(h, ed).todense())
    diag = np.diag(np.array(MO2matrix(sum(fc.N), ed).todense())).real

    def e0(nn):
        sel = np.abs(diag - nn) < 1e-9
        return float(np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0])

    assert info["mu"] == pytest.approx(0.5 * (e0(nf + 1) - e0(nf - 1)), abs=TOL)
    assert info["gap"] == pytest.approx(e0(nf + 1) + e0(nf - 1) - 2 * e0(nf),
                                        abs=TOL)
    # With the mu shift the two families separate cleanly: every particle
    # pole sits at +gap/2 or above and every hole pole at -gap/2 or below,
    # which is the whole reason the shift is applied by default.
    npart = len(info["particle"]["energies"])
    part, hole = info["poles"][:npart], info["poles"][npart:]
    assert part.min() == pytest.approx(0.5 * info["gap"], abs=TOL)
    assert hole.max() == pytest.approx(-0.5 * info["gap"], abs=TOL)
    assert part.min() > 0 > hole.max()


def test_hubbard_mott_gap_on_itensor_v3():
    """The headline case, on the compiled backend: a half-filled Hubbard
    chain's spectral function, whose gap must match E_0^{N+1}+E_0^{N-1}
    -2E_0^N from sector-restricted ED, and whose mu is U/2 by
    particle-hole symmetry."""
    n, u = 4, 8.0
    fc, h = hubbard_chain(n, u=u)
    assert fc.itensor_version == 3
    fc.set_conserved_sector(Nf=n, Sz=0)
    x, A, info = fc.get_spectral_function(
        i=1, spin="up", nex=4, es=np.linspace(-10, 10, 200), delta=0.3,
        return_poles=True, quiet=True)
    assert info["particle"]["sector"] == {"Nf": n + 1, "Sz": 1}
    assert info["hole"]["sector"] == {"Nf": n - 1, "Sz": -1}
    assert info["mu"] == pytest.approx(u / 2, abs=1e-4)  # particle-hole symmetric
    ed = fc.get_ED_obj()
    H = np.array(MO2matrix(h, ed).todense())
    dn = np.diag(np.array(MO2matrix(sum(fc.Ntot), ed).todense())).real
    dz = np.diag(np.array(MO2matrix(2 * sum(fc.Sz), ed).todense())).real

    def e0(nf, sz):
        sel = (np.abs(dn - nf) < 1e-9) & (np.abs(dz - sz) < 1e-9)
        return float(np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0])

    gap = e0(n + 1, 1) + e0(n - 1, -1) - 2 * e0(n, 0)
    assert info["gap"] == pytest.approx(gap, abs=1e-4)


def test_spin_spectral_function_channels_match_ed():
    """All three Sz channels against the exact ED correlator, and their
    documented combination S^zz+(S^{+-}+S^{-+})/2."""
    n = 6
    sc, h = heisenberg_chain(n)
    es = np.linspace(-0.2, 4.0, 300)
    x, S, info = sc.get_spin_spectral_function(
        i=n // 2, nex=100, es=es, delta=0.1, return_poles=True, quiet=True)
    assert info["channels"]["pm"]["sector"] == {"Sz": 2}
    assert info["channels"]["mp"]["sector"] == {"Sz": -2}
    assert info["channels"]["zz"]["sector"] == {"Sz": 0}
    Sp, Sm, Sz = (name2MO(k, sc) for k in ("Sp", "Sm", "Sz"))
    ref = spinchain.Spin_Chain([2] * n)
    ref.set_hamiltonian(h)
    total = 0
    for key, (A, B) in (("pm", (Sm[n // 2], Sp[n // 2])),
                        ("mp", (Sp[n // 2], Sm[n // 2])),
                        ("zz", (Sz[n // 2], Sz[n // 2]))):
        xe, ye = ref.get_dynamical_correlator(mode="ED", submode="ED",
                                              name=(A, B), es=es, delta=0.1)
        assert np.max(np.abs(info["curves"][key].real - ye)) < 1e-4
        total = total + (ye if key == "zz" else 0.5 * ye)
    assert np.max(np.abs(S.real - total)) < 1e-4


def test_longitudinal_channel_elastic_pole_is_optional():
    """S^zz's target sector is the reference sector, so its n=0 state is
    the ground state itself and the raw Lehmann sum carries an elastic
    w=0 pole (which the exact ED reference has too). connected=True drops
    it."""
    n = 6
    sc, _ = heisenberg_chain(n)
    es = np.linspace(-0.2, 4.0, 100)
    common = dict(i=0, nex=20, es=es, delta=0.1, return_poles=True, quiet=True)
    raw = sc.get_spin_spectral_function(**common)[2]
    con = sc.get_spin_spectral_function(connected=True, **common)[2]
    assert np.min(np.abs(raw["channels"]["zz"]["energies"])) < 1e-8
    assert np.min(np.abs(con["channels"]["zz"]["energies"])) > 1e-8
