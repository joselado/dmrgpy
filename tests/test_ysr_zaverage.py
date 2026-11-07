"""z-averaged STM spectrum on a superconductor (ysrtk/zaverage.py): the
log-Gaussian kernel, the gap-edge map, the z average of the BdG chain
against the BCS density of states, and the TASK poles against BdG."""
import numpy as np
import pytest

from dmrgpy import wilsonchain, ysr
from dmrgpy.ysrtk import zaverage as za


@pytest.mark.parametrize("b", [0.2, 0.5, 1.0])
def test_log_gaussian_has_unit_area(b):
    x = np.logspace(-8, 4, 200001)
    for xp in (1e-3, 1.0):
        assert np.trapezoid(za.log_gaussian(x, xp, b), x) == pytest.approx(1, abs=1e-7)


def test_log_gaussian_keeps_a_flat_spectrum_flat():
    """Poles on a logarithmic mesh with weights proportional to the energy
    are a flat density, and the shift b/4 keeps it flat."""
    xp = np.exp(np.arange(-30, 10, 0.1))
    for b in (0.3, 0.6):
        A = za._half_line_density(np.array([1e-3, 1e-2, 0.1]), xp, 0.1*xp, b, 1e-12)
        assert np.allclose(A, 1.0, atol=1e-9)


@pytest.mark.parametrize("xp", [0.0, 1e-6, 5e-3, 2e-2])
def test_crossover_keeps_the_weight(xp):
    """Below omega0 the kernel blends into a Gaussian, at the pole, so
    every pole keeps exactly its weight."""
    g = np.concatenate([np.logspace(-12, -4, 20000), np.linspace(0, 1, 400001)[1:]])
    A = za._half_line_density(g, np.array([xp]), np.array([1.0]), 0.5, 1e-2)
    assert np.trapezoid(A, g) == pytest.approx(1, abs=1e-6)


def test_gap_map_keeps_the_weight_out_of_the_gap():
    Delta = 0.1
    om = np.array([-0.5, -0.12, -0.03, 0.03, 0.1001, 0.3])
    wt = np.array([0.2, 0.1, 0.05, 0.05, 0.3, 0.3])
    xi = np.logspace(-7, 1, 20001)
    es = np.sqrt(xi**2 + Delta**2)
    Ap, sub = za.broaden_poles(es, om, wt, b=0.5, omega0=1e-3, gap=Delta)
    Am, _ = za.broaden_poles(-es, om, wt, b=0.5, omega0=1e-3, gap=Delta)
    # the weight above the gap, integrated in xi where the map is smooth
    jac = xi/es
    assert np.trapezoid(Ap*jac, xi) == pytest.approx(0.6, abs=1e-4)
    assert np.trapezoid(Am*jac, xi) == pytest.approx(0.3, abs=1e-4)
    inside, _ = za.broaden_poles(np.linspace(-0.0999, 0.0999, 101), om, wt,
                                 b=0.5, omega0=1e-3, gap=Delta)
    assert np.max(np.abs(inside)) == 0.0
    assert np.allclose(np.sort(sub["omega"]), [-0.03, 0.03])


def test_normal_state_keeps_the_weight_on_both_sides():
    """At gap = 0 the map is the identity; a pole at w = 0 goes half to
    each side."""
    om = np.array([0.0, 0.02, -0.3])
    wt = np.array([0.4, 0.5, 0.1])
    x = np.concatenate([np.logspace(-9, -4, 4000), np.linspace(1e-4, 5, 200001)[1:]])
    Ap, sub = za.broaden_poles(x, om, wt, b=0.4, omega0=1e-3)
    Am, _ = za.broaden_poles(-x, om, wt, b=0.4, omega0=1e-3)
    assert np.trapezoid(Ap, x) == pytest.approx(0.7, abs=1e-5) # quadrature
    assert np.trapezoid(Am, x) == pytest.approx(0.3, abs=1e-5)
    assert len(sub["omega"]) == 0


@pytest.mark.parametrize("nz,b,tol", [(4, 0.4, 1e-5), (2, 0.6, 2e-3)])
def test_z_average_of_the_bdg_chain_is_the_bcs_density(nz, b, tol):
    """The z-averaged first-site density of states of the pairing Wilson
    chain, times A_Lambda, is rho |w|/sqrt(w^2-Delta^2) per spin: measured
    3.9e-6 at nz=4, b=0.4 and 1.2e-3 at nz=2, b=0.6 (Lambda=4, n=20),
    between the resolution of the chain and D/Lambda."""
    Delta, n = 1e-3, 20
    xi = np.logspace(-4, -1, 200)
    es = np.sqrt(xi**2 + Delta**2)
    A = za.bdg_chain_dos(es, n, Delta, Lambda=4.0, nz=nz, b=b)
    exact = 0.5*ysr.bcs_dos(es, Delta)
    assert np.max(np.abs(A/exact - 1)) < tol


def _shiba_ldos(w, v, Delta, rho=0.5):
    """Spin-up density of states at a classical spin, exchange field v,
    in a flat BCS band: G = [(w - Delta tau_x)/(pi rho s) - v]^-1,
    s = sqrt(Delta^2 - w^2), continued above the gap."""
    s = -1j*np.sqrt(w*w - Delta**2 + 0j)*np.sign(w)
    a = 1/(np.pi*rho*s)
    return -np.imag((a*w - v)/((a*w - v)**2 - (a*Delta)**2))/np.pi


@pytest.mark.parametrize("nz,b,edge", [(8, 0.3, np.exp(1.5*0.09)),
                                       (16, 0.15, np.exp(1.5*0.0225))])
def test_z_average_of_a_classical_spin_is_the_shiba_continuum(nz, b, edge):
    """A static field v on the first chain site (v*A_Lambda on the chain)
    against the Shiba formula above the gap. Far from the edge they agree
    to 1e-2; at the edge the density goes as xi^2 and the log-Gaussian is
    high by exp(3 b^2/2), which we pin."""
    Delta, n, L, alpha = 1e-3, 24, 4.0, 0.5
    v = alpha/(np.pi*0.5)
    xi = np.array([3e-5, 1e-3, 3e-3, 1e-2])
    es = np.sqrt(xi**2 + Delta**2)
    A = np.zeros(len(es))
    for z in za.twists(nz):
        eps, t = wilsonchain.wilson_chain(n, Lambda=L, z=z)
        h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
        H = np.block([[h, Delta*np.eye(n)], [Delta*np.eye(n), -h]])
        H[0, 0] += v*ysr.lambda_factor(L)
        H[n, n] += v*ysr.lambda_factor(L)
        E, V = np.linalg.eigh(H)
        A += za.broaden_poles(es, E, np.abs(V[0])**2, b=b,
                              omega0=abs(t[-1]), gap=Delta)[0]
    r = A/nz*ysr.lambda_factor(L)/_shiba_ldos(es, v, Delta)
    assert r[0] == pytest.approx(edge, abs=3e-3)
    assert np.max(np.abs(r[1:] - 1)) < 1.2e-2


def _merge(E, wt, tol=1e-9):
    """Sum the weights of degenerate poles."""
    out = []
    for e, w in zip(E, wt):
        for o in out:
            if abs(o[0] - e) < tol:
                o[1] += w
                break
        else:
            out.append([e, w])
    return np.array(out)


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_task_poles_are_the_bdg_poles_at_zero_coupling(itensor_version):
    """At J=0 a potential vertex U f^dag/2 sees the local density of
    states of the first chain site, the BdG poles of the chain."""
    n, Delta, U, z = 6, 0.05, 0.7, 0.5
    fc = ysr.kondo_chain(n, J=0.0, Delta=Delta, Lambda=4.0, z=z,
                         itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 20
    p = za.stm_poles(fc, operators=[0.5*U*fc.Cdagup[1], 0.5*U*fc.Cdagdn[1]],
                     nkry=100)
    eps, t = wilsonchain.wilson_chain(n, Lambda=4.0, z=z)
    ref = _merge(*za.bdg_poles(eps, t, Delta))
    for e, w in ref:
        near = np.abs(p["omega"] - e) < 1e-6
        assert np.sum(p["weight"][near]) == pytest.approx(0.5*U*U*w, abs=1e-9)
    assert np.sum(p["weight"]) == pytest.approx(0.5*U*U, abs=1e-9)
    assert p["missing_weight"] < 1e-12


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_z_averaged_spectrum_at_zero_coupling_is_the_bdg_one(itensor_version):
    """The whole route (one chain and one TASK solve per z, the kernel,
    the average and A_Lambda) against bdg_chain_dos."""
    n, Delta, U = 6, 0.05, 1.0
    xi = np.logspace(-2, -0.7, 30)
    es = np.concatenate([-np.sqrt(xi**2 + Delta**2), [0.02],
                         np.sqrt(xi**2 + Delta**2)])
    ops = lambda fc: [0.5*U*fc.Cdagup[1], 0.5*U*fc.Cdagdn[1]]
    with pytest.warns(RuntimeWarning, match="does not reach"): # 6 sites
        _, A, info = za.zaverage_stm_spectrum(n, es, J=0.0, Delta=Delta,
                nz=2, b=0.6, omega0=1e-3, maxm=64, nsweeps=20, nkry=100,
                operators=ops, itensor_version=itensor_version)
    ref = 0.5*U*U*za.bdg_chain_dos(es, n, Delta, Lambda=4.0, nz=2, b=0.6,
                                   omega0=1e-3)
    assert np.max(np.abs(A - ref)) < 1e-9*np.max(ref)
    assert A[len(xi)] == 0.0 # inside the gap there are no poles at J=0
    assert np.all(np.array(info["missing_weight"]) < 1e-12)


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_subgap_pole_is_the_parity_gap(itensor_version):
    """With J the TASK spectrum has one subgap pole at +-|E_odd-E_even|,
    kept out of the log kernel."""
    n, Delta = 6, 0.05
    fc = ysr.kondo_chain(n, J=0.6, Delta=Delta, Lambda=4.0, z=0.25,
                         itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 20
    p = za.stm_poles(fc, U=0.3, nkry=100)
    E = abs(ysr.subgap_energies(fc)["E_ysr"])
    assert E < Delta
    es = np.linspace(-0.3, 0.3, 61)
    _, sub = za.broaden_poles(es, p["omega"], p["weight"], b=0.5, gap=Delta)
    big = np.abs(sub["weight"]) > 1e-8
    assert np.max(np.abs(np.abs(sub["omega"][big]) - E)) < 1e-8
    # electron and hole weights differ at U != 0
    we = np.sum(sub["weight"][sub["omega"] > 0])
    wh = np.sum(sub["weight"][sub["omega"] < 0])
    assert abs(we - wh) > 1e-3


def test_the_sign_of_the_gap_is_a_gauge():
    om = np.array([-0.3, -0.07, -0.03, 0.03, 0.0629, 0.2])
    wt = np.array([0.1, 0.2, 0.1, 0.1, 0.3, 0.2])
    es = np.linspace(-0.2, 0.2, 401)
    Ap, sp = za.broaden_poles(es, om, wt, b=0.4, gap=0.05, delta_sub=1e-3)
    Am, sm = za.broaden_poles(es, om, wt, b=0.4, gap=-0.05, delta_sub=1e-3)
    assert np.array_equal(Ap, Am)
    assert np.array_equal(sp["omega"], sm["omega"])
    assert np.allclose(np.sort(sp["omega"]), [-0.03, 0.03])
    Ac, _ = za.broaden_poles(es, om, wt, b=0.4, gap=-0.05) # no subgap peaks
    assert np.max(np.abs(Ac[np.abs(es) < 0.0499])) == 0.0


def test_default_omega0_makes_the_density_vanish_at_the_edge():
    """A pole whose xi is below the given omega0 takes the Gaussian part,
    a finite xi density at xi = 0, and the map |w|/xi makes it diverge at
    the edge; the default cap puts the pole on the log-Gaussian only."""
    Delta, om, wt = 0.05, np.array([0.0629]), np.array([1.0])
    es = Delta + np.array([1e-17, 1e-12, 1e-8, 1e-4])
    big, _ = za.broaden_poles(es, om, wt, b=0.4, omega0=0.0778, gap=Delta)
    assert big[0] > 1e6 # the documented divergence of a large omega0
    w0 = za._default_omega0(om, Delta, 0.0778)
    assert w0 == pytest.approx(np.sqrt(0.0629**2 - Delta**2))
    A, _ = za.broaden_poles(es, om, wt, b=0.4, omega0=w0, gap=Delta)
    assert np.max(A) < 1e-10 # 3e-14 at Delta + 1e-4, 0 closer


def test_operators_and_U_are_exclusive():
    with pytest.raises(ValueError, match="fold U"):
        za.zaverage_stm_spectrum(4, [0.1], U=0.3, operators=lambda fc: [])


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_short_chain_spectrum_is_finite_and_the_same_for_both_signs(
        itensor_version):
    """5 sites, rho J = 0.075: the subgap pole is at 0.84 Delta, the curve
    is finite on a grid that has a point within roundoff of the edge, and
    the chain, too short for the gap, raises the warning."""
    es = np.linspace(-0.2, 0.2, 401)
    out = []
    for Delta in (0.05, -0.05):
        with pytest.warns(RuntimeWarning, match="does not reach"):
            _, A, info = za.zaverage_stm_spectrum(5, es, J=0.15, Delta=Delta,
                    Lambda=4.0, nz=1, b=0.4, U=0.3, delta_sub=1e-3, maxm=64,
                    nsweeps=20, nkry=100, itensor_version=itensor_version)
        assert np.all(np.isfinite(A)) and np.max(A) < 50.0
        sub = info["subgap"][0]
        big = np.abs(sub["weight"]) > 1e-8
        assert np.allclose(np.sort(np.abs(sub["omega"][big])),
                           0.04189, atol=1e-4)
        out.append(A)
    assert np.max(np.abs(out[0] - out[1])) < 1e-6*np.max(out[0])
    assert 0.1 < np.trapezoid(out[0], es) < 2.0
