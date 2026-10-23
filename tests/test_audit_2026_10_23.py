"""Regressions for docs/audit_2026_10_23_hole_hunt.md, one test group per finding."""
import numpy as np
import pytest

from dmrgpy import spinchain, fermionchain

BACKENDS = [2, 3, "python"]
METTS_BACKENDS = [3, "python"]


def heisenberg(n, v, jz=1.0):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + jz*sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40
    sc.nsweeps = 10
    return sc


def hopping(n, v, mu=0.0):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=v)
    h = sum(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] for i in range(n-1))
    h = h + mu*sum(fc.N[i] for i in range(n))
    fc.set_hamiltonian(h)
    return fc


# Finding 1: METTS with a collapse basis that cannot change a conserved charge

@pytest.mark.parametrize("v", METTS_BACKENDS)
def test_metts_refuses_basis_frozen_by_conserved_charge(v):
    fc = hopping(4, v, mu=0.3)
    with pytest.raises(ValueError, match="conserves"):
        fc.metts_vev(fc.N[0], T=1.0, nsamples=5, nwarmup=2, seed=1, basis_ops=("N",))
    sc = heisenberg(4, v)
    with pytest.raises(ValueError, match="conserves"):
        sc.metts_vev(sc.Sz[0], T=1.0, nsamples=5, nwarmup=2, seed=1, basis_ops=("Sz",))


@pytest.mark.parametrize("v", METTS_BACKENDS)
def test_metts_refuses_non_hermitian_basis(v):
    fc = hopping(4, v)
    with pytest.raises(ValueError, match="not Hermitian"):
        fc.metts_vev(fc.N[0], T=1.0, nsamples=5, nwarmup=2, seed=1, basis_ops=("N", "Cdag"))


@pytest.mark.parametrize("v", METTS_BACKENDS)
def test_metts_accepts_basis_that_breaks_the_charge(v):
    # control: the default basis on a conserving chain, and a diagonal
    # basis on a chain that does not conserve it (transverse-field Ising)
    sc = heisenberg(4, v)
    sc.metts_vev(sc.Sz[0]*sc.Sz[1], T=1.0, nsamples=3, nwarmup=1, seed=1)
    h = sum(sc.Sz[i]*sc.Sz[i+1] for i in range(3)) + sum(sc.Sx[i] for i in range(4))
    sc.set_hamiltonian(h)
    sc.metts_vev(sc.Sz[0], T=1.0, nsamples=3, nwarmup=1, seed=1, basis_ops=("Sz",))


# Finding 2: vev and excited_vev dropping keywords

@pytest.mark.parametrize("v", BACKENDS)
def test_dmrg_vev_refuses_temperature_and_unknown_keywords(v):
    sc = heisenberg(4, v)
    O = sc.Sz[0]*sc.Sz[1]
    with pytest.raises(NotImplementedError):
        sc.vev(O, T=2.0)
    with pytest.raises(TypeError):
        sc.vev(O, bogus=1)
    # the ED route still takes T, and T=0 on DMRG is the ground state
    assert sc.vev(O, mode="ED", T=2.0).real == pytest.approx(-0.0344746681073044, abs=1e-8)
    assert sc.vev(O, T=0.).real == pytest.approx(sc.vev(O, mode="ED").real, abs=1e-6)


@pytest.mark.parametrize("v", BACKENDS)
def test_excited_vev_is_per_state(v):
    sc = heisenberg(4, v, jz=0.7) # XXZ: no degenerate multiplets among the lowest three
    O = sc.Sz[0]*sc.Sz[1]
    ed = np.real(sc.excited_vev(O, n=3, mode="ED"))
    dm = np.real(sc.excited_vev(O, n=3))
    assert len(ed) == 3 and len(dm) == 3
    assert ed[0] == pytest.approx(sc.vev(O, mode="ED").real, abs=1e-8)
    assert np.max(np.abs(ed[1:]-ed[0])) > 1e-2 # used to repeat the ground state
    assert np.max(np.abs(dm-ed)) < 1e-5


# Finding 3: Spin_Chain.get_correlator ignoring name=

@pytest.mark.parametrize("v", BACKENDS)
def test_get_correlator_refuses_name(v):
    sc = heisenberg(4, v)
    with pytest.raises(TypeError):
        sc.get_correlator([(0, 1)], name="SS")
    zz = sc.get_correlator([(0, 1)])
    assert zz[0].real == pytest.approx(sc.vev(sc.Sz[0]*sc.Sz[1]).real, abs=1e-8)


# Finding 4: bond entropy of a state that is not unit norm

@pytest.mark.parametrize("v", BACKENDS)
def test_bond_entropy_is_scale_invariant(v):
    sc = heisenberg(6, v, jz=0.7)
    wf = sc.get_gs()
    s1 = [wf.get_bond_entropy(b-1, b) for b in range(1, 6)]
    w2 = 2.0*wf
    s2 = [w2.get_bond_entropy(b-1, b) for b in range(1, 6)]
    assert np.max(np.abs(np.array(s2)-s1)) < 1e-10
    assert max(s1) > 0.1
