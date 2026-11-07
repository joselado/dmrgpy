"""A spin chain on a superconducting island that leaks into a normal metal
(ysr.island_chain): the exact reductions, the sector energies and the STM
spectrum against ED on "python" and 3. Every number here was reached on
d22 with BLAS on one thread."""
import warnings

import numpy as np
import pytest

from dmrgpy import wilsonchain, ysr
from dmrgpy.ysrtk import paritystate, zaverage


def _heisenberg_energy(ns, Js, B):
    """Ground-state energy of an open Heisenberg chain of ns spins 1/2 in a
    field B along z, by dense diagonalization."""
    sx = np.array([[0, 1], [1, 0]])/2
    sy = np.array([[0, -1j], [1j, 0]])/2
    sz = np.diag([0.5, -0.5])
    def op(o, i):
        m = np.eye(1)
        for k in range(ns):
            m = np.kron(m, o if k == i else np.eye(2))
        return m
    H = np.zeros((2**ns, 2**ns), complex)
    for i in range(ns - 1):
        for o in (sx, sy, sz):
            H += Js*op(o, i) @ op(o, i + 1)
    for i in range(ns):
        H += B*op(sz, i)
    return np.linalg.eigvalsh(H)[0]


def _bdg_energy(eps, t, Delta):
    """Ground-state energy of sum h_ij f^dag_i f_j (both spins) plus
    sum_i Delta_i (f^dag_{i up} f^dag_{i dn} + h.c.): in the Nambu basis
    (f_up, f_dn^dag) it is the sum of the negative Bogoliubov energies
    plus Tr h, from anticommuting f_dn past f_dn^dag."""
    h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    P = np.diag(Delta)
    E = np.linalg.eigvalsh(np.block([[h, P], [P, -h]]))
    return np.sum(E[E < 0]) + np.trace(h), h


def _mode(v):
    return ("ED", 3) if v == "ED" else ("DMRG", v)


@pytest.mark.parametrize("version", ["ED", "python", 3])
def test_one_spin_on_one_site_island_is_kondo_chain(version):
    """ns=1, ni=1 with the island site taken as the first site of the
    Wilson chain and Delta=0 is kondo_chain(n+1) exactly. Run smoke1:
    0 (ED), 2.2e-15 ("python"), 8.9e-16 (3). The Upen defaults differ,
    which does not matter: no term moves charge on the spin site, and the
    one-electron block does not see Upen."""
    mode, v = _mode(version)
    n, J = 3, 0.4
    eps, t = wilsonchain.wilson_chain(n + 1, Lambda=2.0)
    a = ysr.kondo_chain(n + 1, J=J, Delta=0.0, itensor_version=v)
    b = ysr.island_chain(1, 1, n, JK=J, Delta=0.0, eps_island=eps[0],
                         t_c=t[0], wilson=(eps[1:], t[1:]), itensor_version=v)
    for c in (a, b):
        c.maxm, c.nsweeps = 64, 20
    assert b.gs_energy(mode=mode) == pytest.approx(a.gs_energy(mode=mode),
                                                   abs=1e-10)


@pytest.mark.parametrize("ns,B,eps_i", [(2, 0.0, 0.1), (3, 0.07, -0.2)])
def test_decoupled_energy_is_heisenberg_plus_bdg(ns, B, eps_i):
    """At JK=0 the chain is the Heisenberg chain times island plus metal,
    and the fermion part is solved by its Bogoliubov-de Gennes matrix,
    built here from the Hamiltonian, not from fc.bdg_params. It pins every
    sign of the pairing, the hoppings and t_c. Run smoke1: 4.4e-16 and
    2.7e-15."""
    ni, n, Delta, ti, tc = 2, 3, 0.15, 0.4, 0.3
    fc = ysr.island_chain(ns, ni, n, Js=0.6, JK=0.0, Delta=Delta,
                          t_island=ti, eps_island=eps_i, t_c=tc, B=B)
    ew, tw = wilsonchain.wilson_chain(n, Lambda=2.0)
    eps = np.concatenate([[eps_i]*ni, ew])
    t = np.concatenate([[ti]*(ni - 1), [tc], tw])
    Dl = np.concatenate([[Delta]*ni, [0.0]*n])
    e0, h = _bdg_energy(eps, t, 0*Dl) # the Nambu constant, at Delta = 0
    ek = np.linalg.eigvalsh(h)
    assert e0 == pytest.approx(2*np.sum(ek[ek < 0]), abs=1e-12)
    eb, _ = _bdg_energy(eps, t, Dl)
    e = fc.gs_energy(mode="ED")
    assert e == pytest.approx(_heisenberg_energy(ns, 0.6, B) + eb, abs=1e-10)
    assert fc.qp_gap == pytest.approx(ysr.bath_quasiparticle_gap(eps, t, Dl))
    assert fc.island_gap >= abs(Delta)


def test_island_density_of_states_at_zero_exchange_is_bdg():
    """At JK=0 the local density of states of spin up on an island site and
    on the first Wilson site, from ED, is bdg_poles(*fc.bdg_params) pole by
    pole (the same Lorentzians). Run smoke1: 3.7e-13 and 6.1e-13 on peaks
    of about 2.5."""
    ns, ni, n = 2, 2, 3
    fc = ysr.island_chain(ns, ni, n, Js=0.6, JK=0.0, Delta=0.15,
                          t_island=0.4, eps_island=0.1, t_c=0.3)
    es, d = np.linspace(-1.6, 1.6, 81), 0.05
    for k in (ns, ns + ni):
        _, ce = fc.get_dynamical_correlator(mode="ED", submode="CVM",
                    name=(fc.Cup[k], fc.Cdagup[k]), es=es, delta=d)
        _, ch = fc.get_dynamical_correlator(mode="ED", submode="CVM",
                    name=(fc.Cdagup[k], fc.Cup[k]), es=-es[::-1], delta=d)
        A = np.real(ce) + np.real(ch)[::-1]
        E, w = zaverage.bdg_poles(*fc.bdg_params, site=k - ns)
        R = np.sum(w[None, :]*d/np.pi/((es[:, None] - E[None, :])**2 + d*d),
                   axis=1)
        assert np.max(np.abs(A - R)) < 1e-10


def test_bdg_helpers_take_one_pairing_per_site():
    """A scalar Delta and the same value on every site give the same
    matrix, and island_chain's references see the pairing on the island
    only."""
    eps, t = wilsonchain.wilson_chain(6, Lambda=2.0)
    for Delta in (0.1, -0.1):
        assert ysr.bath_quasiparticle_gap(eps, t, Delta) == \
            ysr.bath_quasiparticle_gap(eps, t, np.full(6, Delta))
        a, b = zaverage.bdg_poles(eps, t, Delta, site=2), \
            zaverage.bdg_poles(eps, t, np.full(6, Delta), site=2)
        assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
    fc = ysr.island_chain(2, 2, 6, Delta=0.3)
    eps_f, t_f, D_f = fc.bdg_params
    assert list(D_f) == [0.3, 0.3] + [0.0]*6
    assert len(eps_f) == 8 and len(t_f) == 7 and t_f[1] == fc.t_c
    assert fc.qp_gap < fc.island_gap and not fc.reaches_gap
    assert not hasattr(fc, "impurities") # not a surrogate chain


_SMALL = dict(ns=2, ni=2, n=3, JK=0.5, Delta=0.2, t_c=0.3)
_E_SECTORS = {0: -3.016406563867054, 1: -2.999170588703454} # ED, smoke1


@pytest.mark.parametrize("version", ["ED", "python", 3])
def test_sector_energies_against_ed(version):
    """Both parity sectors of a 7-site chain, with the charge check on
    every spin site. Run smoke1: "python" 1.4e-13 and v3 1e-15 off ED; the
    parity splitting 0.0172360 sits below the JK=0 value qp_gap=0.020468,
    and the ground state is the even singlet."""
    mode, v = _mode(version)
    p = dict(_SMALL)
    fc = ysr.island_chain(p.pop("ns"), p.pop("ni"), p.pop("n"),
                          itensor_version=v, **p)
    fc.maxm, fc.nsweeps = 64, 20
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r = ysr.island_energies(fc, mode=mode)
    assert [str(x.message) for x in w
            if issubclass(x.category, RuntimeWarning)] == []
    for s, e in _E_SECTORS.items():
        assert r["sectors"][s] == pytest.approx(e, abs=1e-9)
    assert r["E_parity"] == pytest.approx(_E_SECTORS[1] - _E_SECTORS[0],
                                          abs=1e-9)
    assert 0 < r["E_parity"] < r["qp_gap"] == fc.qp_gap
    assert r["Sz_ground"] == 0 and r["parity_ground"] == "even"
    assert fc.conserved_sector is None # the chain itself is untouched


def test_kondo_tools_refuse_the_wrong_chain():
    fc = ysr.island_chain(2, 1, 2)
    with pytest.raises(ValueError, match="island_energies"):
        ysr.subgap_energies(fc, mode="ED")
    with pytest.raises(ValueError, match="subgap_energies"):
        ysr.island_energies(ysr.kondo_chain(3, Delta=0.1), mode="ED")


_ES, _DELTA = np.linspace(-1.2, 1.2, 25), 0.06
_STM_ED = {}


def _small_chain(**kw):
    p = dict(_SMALL)
    return ysr.island_chain(p.pop("ns"), p.pop("ni"), p.pop("n"), **p, **kw)


def _stm_ed():
    """The ED spectrum on the end spin of the 7-site chain, computed once."""
    if "a" not in _STM_ED:
        _STM_ED["a"] = ysr.get_stm_spectrum(_small_chain(), _ES,
            delta=_DELTA, mode="ED", submode="CVM")[1]
    return _STM_ED["a"]


def test_default_vertex_is_the_end_spin():
    """The tunneling vertex takes the end spin (fc.tip_spin), not the sum
    of the spins of the chain, and its electron at the first island site.
    On ns=2, ni=1, n=2, 5 sites, so that ED is cheap."""
    ns = 2
    fc = ysr.island_chain(ns, 1, 2, JK=0.5, Delta=0.2, t_c=0.3)
    es, d = _ES, _DELTA
    kw = dict(delta=d, mode="ED", submode="CVM")
    _, a = ysr.get_stm_spectrum(fc, es, **kw)
    _, b = ysr.get_stm_spectrum(fc, es, imp=ns - 1, **kw)
    _, c = ysr.get_stm_spectrum(fc, es, imp=0, **kw)
    assert np.max(np.abs(a - b)) < 1e-12
    assert np.max(np.abs(a - c)) > 1e-2*np.max(a)


@pytest.mark.parametrize("Lambda", [2.0, 4.0])
def test_hopping_into_the_metal_takes_one_square_root_of_lambda_factor(Lambda):
    """The width of a level at e_d << Gamma hopping by v onto the Wilson
    chain, read off its occupation 1/2 - arctan(e_d/Gamma)/pi and averaged
    over 8 twists: v sqrt(A_Lambda) on the chain gives the continuum
    Gamma = pi rho v^2, bare v gives Gamma/A_Lambda. This is the convention
    island_chain states for t_c. Run tc1: 0.991 and 0.994, bare 0.953 and
    0.861."""
    D, v, ed = 1.0, 0.1, 0.005
    A = ysr.lambda_factor(Lambda)
    def width(vc):
        occ = []
        for z in np.linspace(0, 1, 8, endpoint=False):
            eps, t = wilsonchain.wilson_chain(60, Lambda=Lambda, z=z, nint=90)
            tt = np.concatenate([[vc], t])
            h = np.diag(np.concatenate([[ed], eps])) + np.diag(tt, 1) \
                + np.diag(tt, -1)
            E, V = np.linalg.eigh(h)
            occ.append(np.sum(np.abs(V[0, E < 0])**2))
        return ed/np.tan(np.pi*(0.5 - np.mean(occ)))/(np.pi*v*v/(2*D))
    assert width(v*np.sqrt(A)) == pytest.approx(1.0, abs=0.015)
    assert width(v) == pytest.approx(1.0/A, abs=0.015)


@pytest.mark.parametrize("version,tol", [("python", 1e-5), (3, 1e-8)])
def test_stm_spectrum_on_the_parity_state_against_ed(version, tol):
    """TASK on the promoted parity state against ED, with the tangent space
    padded to the whole space (pad_bonds="full"). Run smoke2: 1.5e-6
    ("python") and 3.9e-10 (3) on a peak of 0.680. Without padding both
    backends are 3.3e-5 off, the one-site tangent-space approximation, the
    same on the plain state as on the parity state."""
    fc = _small_chain(itensor_version=version)
    fc.maxm, fc.nsweeps = 64, 20
    c, info = paritystate.parity_ground_state(fc)
    assert info["Sz_ground"] == 0 and "E_ysr" not in info
    assert info["E_parity"] == pytest.approx(_E_SECTORS[1] - _E_SECTORS[0],
                                             abs=1e-9)
    _, dm = ysr.get_stm_spectrum(c, _ES, delta=_DELTA, pad_bonds="full")
    assert np.max(np.abs(dm - _stm_ed())) < tol


_E_NS4 = {0: -3.99108120, 1: -3.77195671} # smoke1, v3 and "python" agree to 6e-9


@pytest.mark.parametrize("version", ["python", 3])
def test_four_spin_sector_solves_keep_one_electron_per_spin(version):
    """ns=4, ni=2, n=4: a sector solve on v3 starts from the seed that hops
    every spin onto the first island site, and on both backends every
    solved state is checked to hold one electron per spin site (a failed
    check raises). Run smoke1 (maxm=48, nsweeps=12): the two backends agree
    to 6e-9 in both sectors, and E_parity = 0.219124 against qp_gap =
    0.219040, the lowest Wilson level of this short chain."""
    fc = ysr.island_chain(4, 2, 4, JK=0.5, Delta=0.2, t_c=0.3,
                          itensor_version=version)
    fc.maxm, fc.nsweeps = 48, 12
    r = ysr.island_energies(fc)
    for s, e in _E_NS4.items():
        assert r["sectors"][s] == pytest.approx(e, abs=5e-8)
    assert r["Sz_ground"] == 0


@pytest.mark.parametrize("version", ["python", 3])
def test_spin_correlations_against_ed(version):
    """<S_a.S_b> over the two spins and the first island site, on the plain
    ground state, against ED. The diagonal is 3/4 on the spin sites."""
    ref = ysr.island_spin_correlations(_small_chain(), mode="ED")
    fc = _small_chain(itensor_version=version)
    fc.maxm, fc.nsweeps = 64, 20
    C = ysr.island_spin_correlations(fc)
    assert C.shape == (3, 3) and np.allclose(np.diag(ref)[:2], 0.75)
    assert ref[0, 1] < 0 and ref[1, 2] < 0 # antiferromagnetic Js and JK
    assert np.max(np.abs(C - ref)) < 1e-7


def test_field_sectors_reach_both_signs():
    """With a field on the spins the lowest state of a parity can sit at
    any |2Sz| <= ns+1, so island_energies solves all of them; on ED."""
    fc = ysr.island_chain(2, 1, 2, JK=0.5, Delta=0.2, t_c=0.3, B=0.3)
    r = ysr.island_energies(fc, mode="ED")
    assert sorted(r["sectors"]) == [-3, -2, -1, 0, 1, 2, 3]
    e = r["sectors"]
    assert e[-2] < e[2] # B > 0 favours Sz < 0
    assert r["E_parity"] == pytest.approx(min(e[-3], e[-1], e[1], e[3])
                                          - min(e[-2], e[0], e[2]))


@pytest.mark.parametrize("ns", [2, 3])
def test_parity_splitting_at_zero_exchange_is_one_quasiparticle(ns):
    """At JK=0 the lowest state of the other parity is the spin-chain
    ground state plus one Bogoliubov quasiparticle of island plus metal,
    so E_parity = +qp_gap for ns even and -qp_gap for ns odd, exactly."""
    fc = ysr.island_chain(ns, 1, 2, Js=0.5, JK=0.0, Delta=0.2, t_c=0.3)
    r = ysr.island_energies(fc, mode="ED")
    sign = 1.0 if ns % 2 == 0 else -1.0
    assert r["E_parity"] == pytest.approx(sign*fc.qp_gap, abs=1e-10)
    assert r["parity_ground"] == ("even" if ns % 2 == 0 else "odd")


@pytest.mark.parametrize("version,tol", [("python", 1e-7), (3, 1e-6)])
def test_sector_solves_on_sixteen_sites_keep_the_spin_charges(version, tol):
    """ns=4, ni=2, n=10: without the warm stage of _seed_impurity_charge a
    "python" sector solve stopped with a spin site empty or doubly occupied
    (run big1, weight 1.0 outside one electron per spin, which the charge
    check turns into a RuntimeError). With it the two sector energies of
    run big2 (maxm=64, nsweeps=16, "python") are -4.86655768 and
    -4.83771241, and v3 is 1.8e-7 and 3.2e-8 above them. v3 starts from a
    random MPS, so its pin is looser. At maxm=48 the even sector is 1.8e-5
    higher (run test4), a truncation of the bond dimension; at JK=0 both
    backends are about 2e-6 above the exact energy at maxm=64."""
    fc = ysr.island_chain(4, 2, 10, JK=0.5, Delta=0.2, t_c=0.3,
                          itensor_version=version)
    fc.maxm, fc.nsweeps = 64, 16
    r = ysr.island_energies(fc) # raises if a spin site lost its electron
    assert r["Sz_ground"] == 0
    assert r["sectors"][0] == pytest.approx(-4.86655768, abs=tol)
    assert r["sectors"][1] == pytest.approx(-4.83771241, abs=tol)


def _island_warnings(fc):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r = ysr.island_energies(fc, mode="ED")
    return r, [str(x.message) for x in w if "island_energies" in str(x.message)]


def test_parity_warning_fires_on_a_splitting_above_one_quasiparticle():
    """The warning is the trip of an odd sector solve that missed its
    ground state. Here it is forced by a tiny qp_gap on the 7-site chain,
    whose ground state is even."""
    fc = _small_chain()
    fc.qp_gap = 1e-6
    with pytest.warns(RuntimeWarning, match="odd sector"):
        ysr.island_energies(fc, mode="ED")


def test_no_parity_warning_on_an_odd_ground_state_bound_by_the_kondo_exchange():
    """On ns=2, ni=1, n=2 at JK=1 the end spin is screened, the ground
    state is odd and lies 1.06 qp_gap below the even one, exactly (ED).
    That is the Kondo binding, not a failed solve, so no warning."""
    fc = ysr.island_chain(2, 1, 2, JK=1.0, Delta=0.2, t_c=0.3)
    r, w = _island_warnings(fc)
    assert w == []
    assert r["parity_ground"] == "odd"
    assert r["E_parity"] < -1.05*fc.qp_gap


@pytest.mark.parametrize("version", ["python", 3])
def test_parity_at_small_pairing_against_ed(version):
    """At Delta=0.02 on the 7-site chain the parity splitting is 6e-6 and
    the pairing is the only term that changes the particle number. On
    "python" 20 sweeps left the even sector 2e-5 high and gave the odd
    parity; 40 sweeps give the even one, as ED does."""
    p = dict(_SMALL, Delta=0.02)
    ns, ni, n = p.pop("ns"), p.pop("ni"), p.pop("n")
    ed = ysr.island_energies(ysr.island_chain(ns, ni, n, **p), mode="ED")
    fc = ysr.island_chain(ns, ni, n, itensor_version=version, **p)
    fc.maxm, fc.nsweeps = 64, 40
    r = ysr.island_energies(fc)
    assert ed["parity_ground"] == "even"
    assert r["parity_ground"] == "even"
    for s in (0, 1):
        assert r["sectors"][s] == pytest.approx(ed["sectors"][s], abs=1e-6)
