"""Regression tests for the `pyitensor` cluster of the third 2026-09-24
hole hunt (docs/audit_2026_09_24c_hole_hunt.md, findings 10 and 13).

  10. The reduced-effort upper band edge, pyitensor Chain._maximum_energy,
      ran its local solves under the ground-state solver's floor of 200
      Krylov steps, and at the SU(2)-symmetric top of a Heisenberg spectrum
      the value-criterion Lanczos used up to 193 of them, so the bound cost
      5 to 9 full ground-state solves on 24 sites. It is capped at 20 now,
      which gives Emax to 1e-10.
  13. to_mpo's first truncating sweep ran on the uncanonicalized automaton,
      so any operator whose coefficients were all below about 2e-7 became a
      different operator: a lone 1e-7*Sz0 was the zero MPO, and a Heisenberg
      chain in units of 3e-7 was built 89 to 95 per cent wrong. The first
      sweep is exact now.

The MPO checks compare against AutoMPO.dense_matrix(), the independent
Kronecker reference; the band-edge checks count the Krylov dimension asked
for rather than timing anything, and compare Emax against ED.
"""

import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain
from dmrgpy.pyitensor import dmrg as dmrgmod
from dmrgpy.pyitensor import mpobuilder as mb
from dmrgpy.pyitensor.autompo import AutoMPO
from dmrgpy.pyitensor.sites import SiteX
from dmrgpy.pyitensor.tensor import contract_many


# ------------------------------------------------------------ finding 13

def _mpo_dense(mpo, sites):
    n = mpo.length()
    T = contract_many([mpo.A(i) for i in range(1, n + 1)])
    si = [sites.si(i) for i in range(1, n + 1)]
    arr = np.asarray(T.transpose_to([i.prime(1) for i in si] + si))
    dim = int(np.prod([i.dim for i in si]))
    return arr.reshape(dim, dim)


def _bond_dims(mpo):
    return [mpo.A(i).inds[-1].dim for i in range(1, mpo.length())]


def _heis(n, s):
    return [(s, [(op, i), (op, i + 1)]) for i in range(1, n) for op in ("Sx", "Sy", "Sz")]


def _relerr(n, terms):
    sites = SiteX([2]*n)
    a = AutoMPO.from_terms(sites, terms)
    ref = a.dense_matrix()
    m = mb.to_mpo(a, cutoff=1e-14)
    return np.linalg.norm(_mpo_dense(m, sites) - ref)/np.linalg.norm(ref), m


@pytest.mark.parametrize("site", [1, 4])
@pytest.mark.parametrize("eps", [1e-6, 3e-7, 2e-7, 1e-7, 3e-8, 1e-10])
def test_lone_small_term_is_built_exactly(site, eps):
    """A lone eps*Sz was the zero MPO for eps <= 2e-7 (relative error 1.0)."""
    err, _ = _relerr(6, [(eps, [("Sz", site)])])
    assert err < 1e-12


@pytest.mark.parametrize("n", [6, 8])
@pytest.mark.parametrize("s", [1.0, 1e-6, 3e-7, 1e-7])
def test_hamiltonian_in_small_units_is_built_exactly(n, s):
    """A Heisenberg chain at s <= 3e-7 was built 89 to 95 per cent wrong, at
    bond dimension 3 instead of 5."""
    err, m = _relerr(n, _heis(n, s))
    assert err < 1e-12
    assert max(_bond_dims(m)) == 5


@pytest.mark.parametrize("field", [False, True])
def test_bond_dimension_is_unchanged(field):
    """The exact first sweep leaves the final compression as it was: the
    same bond dimensions as the reference sum-of-terms builder."""
    n = 14
    terms = _heis(n, 1.0) + ([(0.3, [("Sz", 1)])] if field else [])
    a = AutoMPO.from_terms(SiteX([2]*n), terms)
    assert _bond_dims(mb.to_mpo(a, cutoff=1e-14)) \
        == _bond_dims(mb._sum_of_term_mpos(a, cutoff=1e-14))


def test_relative_truncation_inside_a_sum_is_unchanged():
    """A term 1e-7 below the largest one in a sum is still dropped, as by
    the builder before e448699: a relative weight cutoff of 1e-14 is a
    relative amplitude cutoff of 1e-7, by design."""
    err, _ = _relerr(6, [(1.0, [("Sz", 2)]), (1e-8, [("Sz", 1)])])
    assert err == pytest.approx(1e-8, rel=0.1)


@pytest.fixture(scope="module")
def heisenberg6():
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version="python")
    h = 0
    for i in range(5):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h + 0.3*sc.Sz[0])
    sc.maxm, sc.nsweeps = 40, 12
    np.random.seed(0)
    sc.gs_energy()
    return sc


def test_vev_of_a_small_operator_scales(heisenberg6):
    """vev(1e-7*Sz0) was exactly 0 on "python" (and -0.189361e-7 on v2, v3
    and ED)."""
    sc = heisenberg6
    ref = sc.vev(sc.Sz[0]).real
    assert ref == pytest.approx(-0.189361, abs=1e-5)
    assert sc.vev(1e-7*sc.Sz[0]).real/1e-7 == pytest.approx(ref, rel=1e-10)


def test_kpm_correlator_of_small_operators_scales(heisenberg6):
    """C[1e-7*Sz0, 1e-7*Sz3] was identically 0, with no raise, since the
    moment guard reads a zero operator as zero moments."""
    sc = heisenberg6
    es = np.linspace(-0.5, 3.0, 50)
    _, y = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[3]), es=es, delta=0.2)
    _, ys = sc.get_dynamical_correlator(name=(1e-7*sc.Sz[0], 1e-7*sc.Sz[3]),
                                        es=es, delta=0.2)
    peak = np.max(np.abs(y))
    assert peak > 0.1
    assert np.max(np.abs(ys/1e-14 - y)) < 1e-8*peak


# ------------------------------------------------------------ finding 10

def _niter_record(monkeypatch):
    """Record the Krylov dimension every local ground-state solve asks for."""
    asked = []
    orig = dmrgmod._lanczos_ground_state

    def spy(matvec, v0, niter=30, **kw):
        asked.append(niter)
        return orig(matvec, v0, niter=niter, **kw)
    monkeypatch.setattr(dmrgmod, "_lanczos_ground_state", spy)
    return asked


def _upper_edge(sc):
    sess = sc._session
    sess._bandwidth_max = None
    return sess._maximum_energy()


def _ed_top(sc):
    return np.max(np.linalg.eigvalsh(sc.get_ED_obj().get_hamiltonian().toarray()))


def test_upper_edge_is_capped_and_ground_state_floor_is_kept(monkeypatch):
    """The bound solve asks for 20 Krylov steps at most; the ground-state
    solve still asks for 200. On the 12-site Heisenberg chain Emax is the
    ferromagnetic top 0.25*(n-1) = 2.75."""
    n = 12
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 6
    np.random.seed(1)
    asked = _niter_record(monkeypatch)
    sc.gs_energy()
    assert asked and min(asked) == 200
    del asked[:]
    emax = _upper_edge(sc)
    assert asked and max(asked) == 20
    assert emax == pytest.approx(0.25*(n - 1), abs=1e-8)


@pytest.mark.parametrize("model", ["spin1", "field_heisenberg", "fermions"])
def test_upper_edge_matches_ed(model):
    """The capped bound stays exact on models without the SU(2) top."""
    np.random.seed(2)
    if model == "spin1":
        sc = spinchain.Spin_Chain(["S=1"]*5, itensor_version="python")
        h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
                for i in range(4)) + 0.2*sc.Sz[0]*sc.Sz[0]
    elif model == "field_heisenberg":
        sc = spinchain.Spin_Chain(["S=1/2"]*8, itensor_version="python")
        h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
                for i in range(7)) + 0.5*sum(sc.Sz[i] for i in range(8))
    else:
        sc = fermionchain.Fermionic_Chain(8, itensor_version="python")
        h = sum(sc.Cdag[i]*sc.C[i+1] + sc.Cdag[i+1]*sc.C[i] for i in range(7)) \
            + 0.7*sum(sc.N[i]*sc.N[i+1] for i in range(7)) + 0.2*sc.N[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 6
    sc.gs_energy()
    assert _upper_edge(sc) == pytest.approx(_ed_top(sc), abs=1e-7)
