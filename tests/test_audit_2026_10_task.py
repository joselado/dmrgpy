"""Regressions of the 2026-10 hole hunt (docs/audit_2026_10_hole_hunt.md),
cluster task: submode="TASK" reports the sum-rule weight it drops
(finding 1) and checks its arguments (finding 6)."""
import warnings
import numpy as np
import pytest
from dmrgpy import fermionchain, spinchain, taskdmrg


def kitaev(mu, itensor_version, n=10, maxm=4, t=1.0, D=0.8):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h - t*(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]) \
              + D*(fc.C[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.Cdag[i])
    for i in range(n): h = h - mu*fc.N[i]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps, fc.cutoff = maxm, 40, 0.
    return fc


def heisenberg(n=10, maxm=4):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h + 0.3*sc.Sz[0])
    sc.maxm, sc.nsweeps, sc.cutoff = maxm, 20, 0.
    return sc


def poles_and_warnings(chain, A, B):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        p = taskdmrg.get_poles(chain, name=(A, B))
    hits = [x for x in w if issubclass(x.category, RuntimeWarning)
            and "sum rule" in str(x.message)]
    return p, hits


def odd_term_chain(itensor_version, n=10, maxm=4):
    """Interacting hopping chain with an odd term Cdag_4 + C_4: a unique
    ground state that mixes the two fermion parities whatever the start"""
    fc = fermionchain.Fermionic_Chain(n, itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] \
              + fc.N[i]*fc.N[i+1]
    fc.set_hamiltonian(h + 0.3*fc.N[0] + 0.5*(fc.Cdag[4] + fc.C[4]))
    fc.maxm, fc.nsweeps, fc.cutoff = maxm, 40, 0.
    return fc


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_fermion_operator_on_a_parity_mixture_reports_the_missing_weight(
        itensor_version):
    # at maxm=4 a quarter of the weight of C_9 lies outside the one-site
    # tangent space (0.107 missing of a sum rule of 0.426); the topological
    # Kitaev chain loses most of it (0.016 of 0.63), but which mixture of
    # its degenerate pair DMRG reaches depends on the random start, so
    # the test uses a ground state that is a mixture by itself
    fc = odd_term_chain(itensor_version)
    A, B = fc.Cdag[9], fc.C[9]
    p, hits = poles_and_warnings(fc, A, B)
    ab = fc.vev(A*B).real
    kept = (np.sum(p["weight"]) + p["elastic"]).real
    assert p["missing_weight"].real == pytest.approx(ab - kept, abs=1e-10)
    assert p["missing_weight"].real > 0.05
    assert len(hits) == 1


def test_parity_definite_state_drops_nothing_even_at_maxm_2():
    fc = kitaev(3.0, "python", maxm=2) # trivial phase, unique ground state
    p, hits = poles_and_warnings(fc, fc.Cdag[9], fc.C[9])
    assert abs(p["missing_weight"]) < 1e-8
    assert not hits


def test_one_site_spin_operator_drops_nothing_at_truncated_maxm():
    sc = heisenberg()
    p, hits = poles_and_warnings(sc, sc.Sz[5], sc.Sz[5])
    assert abs(p["missing_weight"]) < 1e-8
    assert not hits


def test_bond_operator_warns():
    sc = heisenberg()
    op = sc.Sx[4]*sc.Sx[5]
    p, hits = poles_and_warnings(sc, op, op)
    assert abs(p["missing_weight"]) > 1e-3
    assert len(hits) == 1


@pytest.mark.parametrize("nkry", [2.7, True, 0, "x"])
def test_nkry_must_be_a_positive_integer(nkry):
    sc = heisenberg(n=4, maxm=8)
    with pytest.raises(ValueError):
        sc.get_dynamical_correlator(submode="TASK", name=(sc.Sz[0], sc.Sz[0]),
                                    nkry=nkry)


def test_es_none_raises_and_a_scalar_es_is_one_frequency():
    sc = heisenberg(n=4, maxm=8)
    name = (sc.Sz[0], sc.Sz[0])
    with pytest.raises(ValueError):
        sc.get_dynamical_correlator(submode="TASK", name=name, es=None)
    with pytest.raises(ValueError):
        sc.get_dynamical_correlator(submode="TASK", name=name,
                                    es=np.ones((2, 2)))
    es, y = sc.get_dynamical_correlator(submode="TASK", name=name, es=0.5)
    _, y2 = sc.get_dynamical_correlator(submode="TASK", name=name,
                                        es=np.array([0.5]))
    assert y.shape == (1,) and y[0] == pytest.approx(y2[0], abs=1e-12)
