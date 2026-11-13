"""Tests of nrgspectral: the complete-basis NRG with the full density
matrix for the STM spectrum of a spin cluster on a Wilson chain
(milestones N0 and N1 of docs/nrg_mps_kondo_plan.md)."""

import numpy as np
import pytest

from dmrgpy import nrgspectral, nrgthermo, wilsonchain, ysr

A4 = ysr.lambda_factor(4.0)
ZS = np.array([0.3+0.05j, -0.2+0.01j, 0.01+0.003j, 1.5+0.2j, -0.05+0.02j,
               0.002+0.001j])


def dense_model(ns, site, Js, JK, B, U, eps, t):
    """H and the two tip vertices X^dag_sigma of ns spins 1/2 and a Wilson
    chain (eps, t), dense, with Jordan-Wigner modes (site, spin), up
    before dn, after the spins."""
    n, nm = len(eps), 2*len(eps)
    Zm, a = np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    ds, df = 2**ns, 2**nm
    c = []
    for j in range(nm):
        m = np.eye(ds)
        for o in [Zm]*j + [a] + [np.eye(2)]*(nm - j - 1):
            m = np.kron(m, o)
        c.append(m)
    sz, sp = np.diag([0.5, -0.5]), np.array([[0.0, 1.0], [0.0, 0.0]])

    def spin(o, i):
        return np.kron(np.kron(np.kron(np.eye(2**i), o),
                               np.eye(2**(ns - i - 1))), np.eye(df))

    Sz = [spin(sz, i) for i in range(ns)]
    Sp = [spin(sp, i) for i in range(ns)]
    up, dn = c[0], c[1]
    szl, spl = 0.5*(up.T @ up - dn.T @ dn), up.T @ dn
    H = 0.0
    for i in range(ns - 1):
        H = H + Js*(Sz[i] @ Sz[i+1] + 0.5*(Sp[i] @ Sp[i+1].T
                                            + Sp[i].T @ Sp[i+1]))
    H = H + JK*(Sz[site] @ szl + 0.5*(Sp[site] @ spl.T + Sp[site].T @ spl))
    for i in range(ns):
        H = H + B*Sz[i]
    for k in range(n):
        H = H + eps[k]*(c[2*k].T @ c[2*k] + c[2*k+1].T @ c[2*k+1])
    for k in range(n - 1):
        for s in (0, 1):
            h = c[2*k+s].T @ c[2*k+2+s]
            H = H + t[k]*(h + h.T)
    one = np.eye(H.shape[0])
    ops = [0.5*(up.T @ (U*one + Sz[site]) + dn.T @ Sp[site]),
           0.5*(up.T @ Sp[site].T + dn.T @ (U*one - Sz[site]))]
    return H, ops


def dense_poles(H, ops, T, tol):
    """Poles of sum_O sum_ab (p_a + p_b)|<a|O|b>|^2 delta(w - E_a + E_b)."""
    E, V = np.linalg.eigh(H)
    if T > 0:
        p = np.exp(-(E - E[0])/T)
    else:
        p = (E - E[0] <= tol).astype(float)
    p = p/np.sum(p)
    w, W = [], []
    for O in ops:
        M = V.T @ O @ V
        w.append((E[:, None] - E[None, :]).ravel())
        W.append((M*M*(p[:, None] + p[None, :])).ravel())
    return np.concatenate(w), np.concatenate(W)


def resolvent(w, W, zs=ZS):
    return np.array([np.sum(W/(z - w)) for z in zs])


@pytest.mark.parametrize("ns,site,n,B,T,U", [
    (1, 0, 3, 0.0, 0.0, 0.0), (1, 0, 4, 0.0, 0.05, 0.3),
    (1, 0, 3, 0.02, 0.0, 0.2), (2, 1, 3, 0.0, 0.0, 0.0),
    (2, 0, 3, 0.0, 0.0, 0.1), (3, 1, 2, 0.01, 0.1, 0.25)])
def test_untruncated_poles_equal_dense_ed(ns, site, n, B, T, U):
    """keep=None: the poles are those of a dense diagonalization of the
    same model (measured 1e-13 on the resolvent at Im z down to 1e-3)."""
    eps, t = wilsonchain.wilson_chain(n, Lambda=4.0, band=(-1.0, 1.0), z=0.3)
    r = nrgspectral.chain_poles(ns=ns, n=n, site=site, Js=0.3, JK=0.4, B=B,
                                U=U, wilson=(eps, t), keep=None, T=T,
                                bins_per_decade=None, Lambda=4.0)
    H, ops = dense_model(ns, site, 0.3, 0.4, B, U, eps, t)
    w, W = dense_poles(H, ops, T, nrgspectral.DEG_TOL*r.omega_n[-1])
    assert np.max(np.abs(resolvent(r.omega, r.weight) - resolvent(w, W))) < 1e-11
    assert r.total == pytest.approx(np.sum(W), abs=1e-13)
    assert r.sum_rule == pytest.approx(np.sum(W), abs=1e-13)
    assert r.trace == pytest.approx(1.0, abs=1e-13)
    assert r.E0[-1] == pytest.approx(np.linalg.eigvalsh(H)[0], abs=1e-12)


@pytest.mark.parametrize("ns,site,B,T,keep", [
    (1, 0, 0.0, 0.0, 33), (1, 0, 0.0, 0.3, 40), (1, 0, 0.02, 1.0, 33),
    (2, 1, 0.0, 0.5, 70), (2, 0, 0.0, 0.0, 64), (3, 1, 0.01, 0.4, 130)])
def test_decoupled_tail_truncated_equals_short_chain(ns, site, B, T, keep):
    """Three coupled Wilson sites and three with t = 0 and eps = 0. The
    truncation acts from the last coupled shell on, where every discarded
    state is an exact eigenstate with the energy of its shell, so the
    truncated complete-basis result must be the untruncated one of the
    three-site chain. At T > 0 up to 46% of rho sits on the discarded
    states of the earlier shells, which tests the backward pass and the
    kept-discarded terms (measured 1e-15)."""
    eps, t = wilsonchain.wilson_chain(3, Lambda=4.0, band=(-1.0, 1.0), z=0.2)
    eps6 = np.concatenate([eps, np.zeros(3)])
    t6 = np.concatenate([t, np.zeros(3)])
    kw = dict(ns=ns, site=site, Js=0.3, JK=0.4, B=B, U=0.2, T=T,
              bins_per_decade=None, Lambda=4.0)
    r3 = nrgspectral.chain_poles(n=3, wilson=(eps, t), keep=None, **kw)
    r6 = nrgspectral.chain_poles(n=6, wilson=(eps6, t6), keep=keep, **kw)
    assert np.any(r6.nkept[2:] < 4*r6.nkept[1:-1]) # it does truncate
    assert np.max(np.abs(resolvent(r3.omega, r3.weight)
                         - resolvent(r6.omega, r6.weight))) < 1e-12
    assert r6.total == pytest.approx(r6.sum_rule, abs=1e-13)
    assert r6.trace == pytest.approx(1.0, abs=1e-13)
    if T > 0.2:
        assert np.sum(r6.discarded_weight[:-1]) > 0.1


@pytest.mark.parametrize("ns,T", [(1, 0.0), (1, 1e-4), (2, 0.0)])
def test_sum_rule_at_small_keep(ns, T):
    """Truncated on a real Wilson chain: the weights still add up to
    sum_sigma Tr[rho {X, X^dag}] and Tr rho = 1, to roundoff."""
    r = nrgspectral.chain_poles(ns=ns, n=20, Js=0.05, JK=0.3, U=0.2,
                                Lambda=4.0, keep=60, T=T)
    assert r.total == pytest.approx(r.sum_rule, abs=1e-12)
    assert r.trace == pytest.approx(1.0, abs=1e-12)


def _lorentz(r, es, delta):
    return np.sum(r.weight[None, :]*delta/np.pi
                  /((es[:, None] - r.omega[None, :])**2 + delta**2), axis=1)


@pytest.mark.parametrize("J,B,U,z", [(0.4, 0.0, 0.0, 0.0),
                                     (0.4, 0.03, 0.3, 0.4)])
def test_kondo_chain_against_its_ed_stm_spectrum(J, B, U, z):
    """The conventions of ysr: the poles broadened by the Lorentzian of
    the dynamical correlator give ysr.get_stm_spectrum(mode="ED") of
    ysr.kondo_chain, whose impurity is a penalized fermion site (2e-14)."""
    es, delta = np.linspace(-1.5, 1.5, 201), 0.03
    fc = ysr.kondo_chain(3, J=J, Lambda=4.0, B=B, z=z,
                         itensor_version="python")
    _, A = ysr.get_stm_spectrum(fc, es, delta=delta, U=U, mode="ED",
                                submode="ED")
    r = nrgspectral.chain_poles(ns=1, n=3, JK=J, B=B, U=U, z=z, Lambda=4.0,
                                keep=None, bins_per_decade=None)
    assert np.max(np.abs(A - _lorentz(r, es, delta))) < 1e-12*np.max(A)


@pytest.mark.parametrize("ns,n,site", [(2, 3, 1), (2, 3, 0), (3, 2, 1)])
def test_spin_chain_on_metal_against_its_ed_stm_spectrum(ns, n, site):
    es, delta = np.linspace(-1.5, 1.5, 201), 0.03
    fc = ysr.spin_chain_on_metal(ns, n, site=site, Js=0.3, JK=0.5,
                                 Lambda=4.0, itensor_version="python")
    _, A = ysr.get_stm_spectrum(fc, es, delta=delta, U=0.2, mode="ED",
                                submode="ED")
    r = nrgspectral.chain_poles(ns=ns, n=n, site=site, Js=0.3, JK=0.5, U=0.2,
                                Lambda=4.0, keep=None, bins_per_decade=None)
    assert np.max(np.abs(A - _lorentz(r, es, delta))) < 1e-12*np.max(A)


@pytest.mark.parametrize("U,operators,factor", [(0.0, "tip", 3/8),
                                                (0.4, "tip", 0.08 + 3/8),
                                                (0.0, "ldos", 2.0)])
def test_free_spin_gives_the_density_of_states(U, operators, factor):
    """J = 0: the spin is a free moment, and A(w) is (U^2/2 + 3/8) times
    the density of states per spin of the first Wilson site, at any keep
    (2 for f^dag_sigma itself). The ground manifold of the 5-site chain
    holds the four states of its zero mode and the two of the spin."""
    eps, t = wilsonchain.wilson_chain(5, Lambda=4.0, band=(-1.0, 1.0))
    r = nrgspectral.chain_poles(ns=1, n=5, JK=0.0, U=U, wilson=(eps, t),
                                keep=None, bins_per_decade=None,
                                operators=operators, Lambda=4.0)
    from dmrgpy.ysrtk.zaverage import bdg_poles
    w, wt = bdg_poles(eps, t, 0.0)
    ref = factor*resolvent(w, wt)
    # relative to the largest value, about 50 on the ldos: 2.5e-12 absolute
    # on Triton, 5e-14 of it
    assert np.max(np.abs(resolvent(r.omega, r.weight) - ref)) \
        < 1e-12*np.max(np.abs(ref))
    rk = nrgspectral.chain_poles(ns=1, n=30, JK=0.0, U=U, Lambda=4.0,
                                 keep=100, operators=operators)
    assert rk.total == pytest.approx(factor, abs=1e-12)


def test_unitary_limit_of_the_kondo_resonance():
    """rho J = 0.2 at Lambda = 4 (T_K(Wilson) = 3.49e-3): the T-matrix
    t(w) = pi^2 rho J^2 A(w)/2 of the continuum, z-averaged over 4 meshes
    with b = 0.4, is 1 to 1% from 1e-9 to 1e-4 and falls through 1/2 near
    3 T_K (measured 0.9938 to 1.0097 on that window, 0.593 at 1e-2)."""
    J = 0.4*A4
    es = np.array([1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 1e-2, 0.0])
    _, A, info = nrgspectral.stm_spectrum(es, ns=1, n=40, nz=4, b=0.4, JK=J,
                                          Lambda=4.0, keep=400)
    tw = np.pi**2*0.5*(J/A4)**2*A/2.0
    assert np.max(np.abs(tw[:5] - 1.0)) < 0.015
    assert 0.5 < tw[5] < 0.7
    assert tw[6] == pytest.approx(1.0, abs=0.01) # w = 0, read at the floor
    _, Am, _ = nrgspectral.stm_spectrum(-es, ns=1, n=40, nz=4, b=0.4, JK=J,
                                        Lambda=4.0, keep=400)
    assert np.max(np.abs(Am - A)) < 1e-10*np.max(A) # particle-hole symmetric


def test_free_spin_at_finite_temperature():
    """J = 0 at T > 0: the spectrum does not depend on T. The NRG keeps the
    weight inside |w| < T but not its shape (250 times too much weight
    below 1e-7 at T = 3.5e-4), and the Fermi function of the tip averages
    over that window: dI/dV of the poles equals that of the exact density
    of states of the same chains, to 0.3% (measured 0.994 to 1.002)."""
    from dmrgpy.ysrtk.zaverage import bdg_poles, twists
    T = 3.5e-4
    es = np.array([0.0, 1e-5, 1e-4, 3.5e-4, 1e-3, 3e-3, 1e-2])
    _, A, info = nrgspectral.stm_spectrum(es, ns=1, n=60, nz=2, b=0.4,
                                          JK=0.0, Lambda=4.0, keep=400, T=T,
                                          operators="ldos")
    ref = np.zeros(len(es))
    for z, r in zip(twists(2), info["results"]):
        assert r.stop == "omega" and r.omega_n[-1] <= T/100
        assert r.total == pytest.approx(r.sum_rule, abs=1e-12)
        eps, t = wilsonchain.wilson_chain(60, Lambda=4.0, band=(-1.0, 1.0), z=z)
        w, wt = bdg_poles(eps, t, 0.0)
        ref += nrgspectral.thermal_didv(es, w, 2*wt, T, b=0.4)*A4/2
    assert np.max(np.abs(A/ref - 1.0)) < 0.005


def test_kondo_resonance_at_finite_temperature():
    """rho J = 0.2: dI/dV at zero bias, as the T-matrix of the continuum,
    falls from the unitary value with T (measured 1.001 at T = 0, 0.989 at
    0.1 T_K and 0.720 at T_K, Wilson's T_K = 3.489e-3)."""
    J, TK = 0.4*A4, 3.489e-3
    es = np.array([0.0])
    vals = []
    for T in (0.0, 0.1*TK, TK):
        _, A, _ = nrgspectral.stm_spectrum(es, ns=1, n=60, nz=4, b=0.4, JK=J,
                                           Lambda=4.0, keep=400, T=T,
                                           stop_tol=1e-4 if T == 0 else None)
        vals.append(np.pi**2*0.5*(J/A4)**2*A[0]/2.0)
    assert vals[0] == pytest.approx(1.0, abs=0.01)
    assert vals[1] == pytest.approx(0.989, abs=0.01)
    assert vals[2] == pytest.approx(0.720, abs=0.02)


def test_stop_rule_at_the_fixed_point():
    """The tail stops stop_extra shells after two consecutive shells agree
    with the shells two sites before (30 of 80 sites here), and the
    spectrum is that of the full chain to the precision of one mesh: on
    one mesh at Lambda = 4 and keep = 300 the T = 0 spectrum at a fixed
    energy moves by up to 0.7% with the length of the chain (measured on
    26 to 80 sites, both parities), since the density matrix is the
    ground manifold of the last shell."""
    J = 0.4*A4
    r = nrgspectral.chain_poles(ns=1, n=80, JK=J, Lambda=4.0, keep=300,
                                stop_tol=1e-4, stop_extra=4)
    assert r.stop == "fixed point" and 10 < r.nshells < 80
    full = nrgspectral.chain_poles(ns=1, n=80, JK=J, Lambda=4.0, keep=300)
    es = np.array([1e-8, 1e-6, 1e-4])
    a = nrgspectral.broaden(es, r.omega, r.weight, b=0.8,
                            omega0=r.omega_n[-1])
    f = nrgspectral.broaden(es, full.omega, full.weight, b=0.8,
                            omega0=full.omega_n[-1])
    assert np.max(np.abs(a/f - 1.0)) < 1.5e-2


def test_stm_spectrum_passes_the_number_of_spins():
    _, A1, i1 = nrgspectral.stm_spectrum(np.array([1e-3]), ns=1, n=10, nz=1,
                                         JK=0.4, keep=100)
    _, A3, i3 = nrgspectral.stm_spectrum(np.array([1e-3]), ns=3, n=10, nz=1,
                                         Js=0.1, JK=0.4, keep=100)
    assert i1["results"][0].nkept[0] == 8
    assert i3["results"][0].nkept[0] == 32


def test_kondo_resonance_splits_in_a_field():
    """One spin at rho J = 0.2 with B = 10 T_K on the impurity: the
    spin-summed resonance splits into two peaks just below +-B (measured
    at 9.06 T_K with four meshes) and the zero-bias value falls to 0.21,
    while B = T_K leaves a flat top at zero bias (0.90), with no peak
    above it."""
    J, TK = 0.4*A4, 3.489e-3
    es = np.concatenate([-np.logspace(-1, -5, 33), np.logspace(-5, -1, 33)])
    for B, zero in ((TK, 0.90), (10*TK, 0.21)):
        _, A, info = nrgspectral.stm_spectrum(es, ns=1, n=60, JK=J, B=B,
                                              Lambda=4.0, nz=2, b=0.5,
                                              keep=400, stop_tol=1e-4)
        assert info["results"][0].flip is False
        t = np.pi**2*0.5*(J/A4)**2*A/2.0
        assert t[np.argmin(np.abs(es - 1e-5))] == pytest.approx(zero, abs=0.03)
        pos = es > 0
        peak = es[pos][np.argmax(t[pos])]
        if B > 5*TK:
            assert 0.7*B < peak < 1.05*B
            assert -es[~pos][np.argmax(t[~pos])] == pytest.approx(peak, rel=1e-6)
        else:
            assert np.max(t[pos]) - t[np.argmin(np.abs(es - 1e-5))] < 0.03


def test_validation():
    with pytest.raises(ValueError):
        nrgspectral.spin_cluster_shell(ns=2, site=2)
    with pytest.raises(ValueError):
        nrgspectral.spin_cluster_shell(ns=0)
    with pytest.raises(ValueError):
        nrgspectral.spin_cluster_shell(operators="LDOS")
    with pytest.raises(ValueError):
        nrgspectral.chain_poles(n=4, keep=0)
    with pytest.raises(ValueError):
        nrgspectral.chain_poles(n=4, T=-1.0)
    with pytest.raises(ValueError):
        nrgspectral.chain_poles(n=4, wilson=(np.zeros(3), np.ones(2)))


@pytest.mark.parametrize("T", [1e-12, 1e-14])
def test_low_temperature_does_not_see_a_constant(T):
    """A constant added to the start shell does not move the result at low
    T: the Boltzmann weights use the ground energies relative to the last
    shell, sums of the small increments of the later shells, and not the
    absolute energies, whose roundoff of 1e-16 enters as 1e-16/T (it moved
    A by 1.3e-2 at T = 1e-14, found by an independent review)."""
    eps, t = wilsonchain.wilson_chain(60, Lambda=4.0, band=(-1.0, 1.0))
    es = np.array([-3.0, -1.0, 1.0, 3.0, 10.0])*T
    A = []
    for c in (0.0, 0.37, 100.0):
        start = nrgspectral.spin_cluster_shell(ns=1, JK=0.3*A4, U=0.1,
                                               eps0=eps[0], omega=2*abs(t[0]))
        start.E = {lab: e + c for lab, e in start.E.items()}
        r = nrgspectral.spectral_poles(start, eps[1:], t, keep=200, T=T)
        A.append(nrgspectral.broaden(es, r.omega, r.weight, b=0.6, omega0=T))
    assert np.max(np.abs(A[1]/A[0] - 1)) < 1e-8
    assert np.max(np.abs(A[2]/A[0] - 1)) < 1e-8


def test_a_tiny_field_is_not_averaged_away():
    """The spectra of (Q, 2Sz) and (Q, -2Sz) are averaged only at B = 0
    exactly, read off the start shell, not when they agree to a tolerance
    (a field of 5e-11 used to be erased)."""
    r = nrgspectral.chain_poles(ns=1, n=10, JK=-0.3*A4, B=5e-11, Lambda=4.0,
                                keep=100)
    assert r.flip is False
    r = nrgspectral.chain_poles(ns=1, n=10, JK=-0.3*A4, B=0.0, Lambda=4.0,
                                keep=100)
    assert r.flip is True
