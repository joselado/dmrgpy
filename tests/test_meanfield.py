"""Mean-field decoupling of a spin chain (`dmrgpy.meanfield`).

These lock in a fix. `spinchain_meanfield` read its exchange couplings off
`sc.exchange`, a list that `Spin_Chain.set_exchange()` used to populate.
That builder was removed in favour of writing Hamiltonians out with
`SS(i,j)`/`set_hamiltonian()`, leaving `sc.exchange` permanently the
integer 0 -- so every call raised `TypeError: 'int' object is not
iterable`, on any chain, while the user guide documented the function as
working. It now reads the chain's actual `MultiOperator` Hamiltonian.

Everything here runs under `mode="ED"` on 4-6 sites, so the whole file is
a couple of seconds and needs no compiled backend.
"""
import numpy as np
import pytest

import dmrgpy.spinchain as spinchain
from dmrgpy import meanfield


def _ising_chain(n, J, transverse=0.0):
    """Ising chain with coupling J (J<0 ferromagnetic, J>0 antiferro)."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        h = h + J*sc.Sz[i]*sc.Sz[i+1]
        if transverse != 0.0:
            h = h + transverse*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1])
    sc.set_hamiltonian(h)
    return sc


def _sz(chain):
    return np.array(chain.get_magnetization(mode="ED")).T[:, 2].real


def test_meanfield_runs_at_all():
    """The regression proper: this used to raise TypeError on every call."""
    np.random.seed(0)
    sc = _ising_chain(4, -1.0)
    out = meanfield.spinchain_meanfield(sc, p=0.0, mix=0.5, maxerror=1e-8,
                                         maxite=200, mode="ED")
    assert out is not None


def test_decomposition_is_exact():
    """p=1 switches the Weiss field off, so the rebuilt mean-field
    Hamiltonian must be the original one -- the tightest check that
    reading the MultiOperator back loses nothing. Includes on-site fields
    and anisotropic exchange, neither of which the old sc.exchange list
    could represent faithfully."""
    n = 5
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        h = h + sc.SS(i, i+1)
    h = h + 0.7*sc.Sz[0] - 0.4*sc.Sx[2]     # on-site fields
    h = h + 0.35*sc.Sx[0]*sc.Sz[3]          # anisotropic, non-nearest-neighbour
    sc.set_hamiltonian(h)
    e_exact = sc.gs_energy(mode="ED")

    b, J, const = meanfield.decompose_spin_hamiltonian(sc)
    bonds = [(i, j, a, bb, J[i, j, a, bb])
             for i in range(n) for j in range(n)
             for a in range(3) for bb in range(3)
             if abs(J[i, j, a, bb]) > 1e-10]
    zero_field = np.zeros((n, 3), dtype=np.complex128)
    rebuilt = meanfield._mean_field_hamiltonian(sc, bonds, b, const,
                                                 zero_field, 1.0)
    sc2 = sc.copy()
    sc2.set_hamiltonian(rebuilt)
    assert sc2.gs_energy(mode="ED") == pytest.approx(e_exact, abs=1e-10)


def test_ferromagnetic_ising_saturates():
    np.random.seed(1)
    out = meanfield.spinchain_meanfield(_ising_chain(6, -1.0), p=0.0,
                                         mix=0.5, maxerror=1e-8, maxite=200,
                                         mode="ED")
    sz = _sz(out)
    # every site saturated and all pointing the same way
    assert np.allclose(np.abs(sz), 0.5, atol=1e-6)
    assert np.allclose(sz, sz[0], atol=1e-6)


def test_antiferromagnetic_ising_gives_neel_order():
    """From a staggered start, an antiferromagnetic Ising chain must
    mean-field decouple to perfect Neel order. This is what actually
    pins the sign convention of the Weiss field: get it backwards and
    this converges to the ferromagnet instead."""
    n = 6
    m0 = [[0.0, 0.0, 0.5*(-1)**i] for i in range(n)]
    out = meanfield.spinchain_meanfield(_ising_chain(n, 1.0), p=0.0,
                                         mix=0.5, m0=m0, maxerror=1e-8,
                                         maxite=200, mode="ED")
    sz = _sz(out)
    expected = np.array([0.5*(-1)**i for i in range(n)])
    assert np.allclose(sz, expected, atol=1e-6) or \
           np.allclose(sz, -expected, atol=1e-6)


def test_p_equals_one_reproduces_the_exact_ground_state():
    """At p=1 the many-body exchange is kept in full and the Weiss field
    is switched off, so the solved chain is the original model."""
    np.random.seed(2)
    n = 4
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        h = h + sc.SS(i, i+1)
    sc.set_hamiltonian(h)
    e_exact = sc.gs_energy(mode="ED")
    out = meanfield.spinchain_meanfield(sc, p=1.0, mix=0.5, maxerror=1e-8,
                                         maxite=50, mode="ED")
    assert out.gs_energy(mode="ED") == pytest.approx(e_exact, abs=1e-8)


def test_non_spin_hamiltonian_is_rejected():
    """A fermionic chain has no spin decoupling here, and should say so
    rather than producing a wrong number."""
    import dmrgpy.fermionchain as fermionchain
    fc = fermionchain.Fermionic_Chain(4)
    h = 0
    for i in range(3):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    fc.set_hamiltonian(h)
    with pytest.raises(ValueError):
        meanfield.decompose_spin_hamiltonian(fc)


def test_three_site_term_is_rejected():
    sc = spinchain.Spin_Chain(["S=1/2"]*4)
    sc.set_hamiltonian(sc.Sz[0]*sc.Sz[1]*sc.Sz[2])
    with pytest.raises(ValueError):
        meanfield.decompose_spin_hamiltonian(sc)


def test_removed_builders_are_gone():
    """set_fields() silently replaced the Hamiltonian instead of adding to
    it, and set_swave_pairing() took a one-argument function where every
    sibling builder takes two. Both were removed rather than fixed."""
    import dmrgpy.fermionchain as fermionchain
    sc = spinchain.Spin_Chain(["S=1/2"]*2)
    assert not hasattr(sc, "set_fields")
    assert not hasattr(fermionchain.Spinful_Fermionic_Chain(2),
                        "set_swave_pairing")
    assert not hasattr(fermionchain.Spinful_Fermionic_Chain_Native(2),
                        "set_swave_pairing")
