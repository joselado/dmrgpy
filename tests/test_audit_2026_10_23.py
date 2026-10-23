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


# Found by the guide check of the same day: evolve_and_measure without wf=

@pytest.mark.parametrize("v", BACKENDS)
def test_evolve_and_measure_solves_the_ground_state_when_none_is_given(v):
    from dmrgpy import timedependent
    sc = heisenberg(4, v)
    (ts, sz) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0]*sc.Sz[1], nt=5, dt=0.05)
    # it starts from the ground state, which is stationary; v2's MPO-Taylor
    # integrator drifts by about 4e-6 per step at this dt, the others less
    e = sc.vev(sc.Sz[0]*sc.Sz[1], mode="ED").real
    assert np.real(sz[0]) == pytest.approx(e, abs=1e-8)
    assert np.max(np.abs(np.real(sz)-e)) < 1e-4


# Found by the guide check: METTS on a chain with a conserved sector

@pytest.mark.parametrize("v", METTS_BACKENDS)
def test_metts_refuses_a_conserved_sector(v):
    sc = heisenberg(4, v)
    sc.set_conserved_sector(Sz=0)
    with pytest.raises(ValueError, match="conserved-sector mode"):
        sc.metts_vev(sc.Sz[0], T=1.0, nsamples=3, nwarmup=1, seed=1)
    with pytest.raises(ValueError, match="conserved-sector mode"):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), T=1.0, nt=3, nsamples=2, nwarmup=1, seed=1)


# Found while fixing the broken examples: the two missing modules, the
# CVMimag weight and Spinon_Chain.get_gs

def test_maxent_recovers_a_density_of_its_own_family():
    from dmrgpy.maxenttk.pymaxent import reconstruct
    x = np.linspace(-2, 3, 20001)
    p = np.exp(-0.5*x**2 + 0.3*x - 0.1*x**4)
    p /= np.trapezoid(p, x)
    mu = [np.trapezoid(p*x**k, x) for k in range(5)]
    sol, lam = reconstruct(mu, bnds=[-2, 3])
    assert np.max(np.abs(sol(x)-p)) < 1e-7
    assert sol(np.array([-3., 4.])).tolist() == [0., 0.]


def test_pade_reproduces_a_rational_function_near_the_axis():
    from dmrgpy.padetk import bruteforcepade
    E, W = np.array([0.5, 1.3, 2.2]), np.array([1., 0.5, 0.2])
    G = lambda z: np.sum(W/(z-E))
    zs = np.linspace(-0.5, 3, 20) + 2.0j
    p = bruteforcepade.fit(zs, [G(z) for z in zs])
    xs = np.linspace(0, 3, 7) + 0.05j
    assert max(abs(p(z)-G(z)) for z in xs) < 1e-6


def test_cvmimag_has_the_weight_of_the_exact_spectrum():
    # it returned -2 Im G/pi, twice the weight of every other route
    sc = heisenberg(6, 3)
    es = np.linspace(-0.5, 3, 200)
    name = (sc.Sx[0], sc.Sx[0])
    _, ed = sc.get_dynamical_correlator(mode="ED", submode="ED", name=name, es=es, delta=0.05)
    _, ci = sc.get_dynamical_correlator(submode="CVMimag", name=name, es=es, delta=0.05)
    w_ed, w_ci = np.trapezoid(np.real(ed), es), np.trapezoid(np.real(ci), es)
    assert w_ci == pytest.approx(w_ed, rel=2e-2)
    assert es[np.argmax(np.real(ci))] == pytest.approx(es[np.argmax(np.real(ed))], abs=0.05)


@pytest.mark.parametrize("v", [3, "python"])
def test_spinon_chain_reproduces_the_heisenberg_energy(v):
    n = 6
    fc = fermionchain.Spinon_Chain(n, itensor_version=v)
    fc.set_hamiltonian(sum(fc.Sx[i]*fc.Sx[i+1] + fc.Sy[i]*fc.Sy[i+1]
                           + fc.Sz[i]*fc.Sz[i+1] for i in range(n-1)))
    fc.maxm, fc.nsweeps = 40, 10
    sc = heisenberg(n, v)
    assert fc.gs_energy() == pytest.approx(sc.gs_energy(mode="ED"), abs=1e-6)
