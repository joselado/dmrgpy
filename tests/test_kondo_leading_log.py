"""The leading-log resummation of the STM Kondo spectrum
(kondospectrumtk/leadinglog.py, docs/new_methods/kondo_leading_log/).

Three checks, in the order of the method's validation ladder: the closed
form at S=1/2, B=0; the S=1/2 equations in a field of Rosch, Costi,
Paaske, Woelfle, cond-mat/0301106 eq. 6, integrated independently; and
the first order in Jrho_s against conductance.py's third order on the
models of test_kondo_spectrum_tmatrix.py. Plus the stop at the Kondo
scale, and the pieces of the joint flow (2026-11-11): the block form of
the vertex product, the cut to the columns of the occupied states, the
tip potential scattering that does not flow, the closed form below the
stop, and the values of the two-pass implementation it replaced.
"""
import warnings
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from dmrgpy.kondospectrumtk import conductance, leadinglog
from test_kondo_spectrum_tmatrix import MODELS, EVS, _ks, _spin

W0, G0 = 20e-3, 5e-6


def _K(x, D):
    return leadinglog._kernel(x, D, W0, G0)


def test_spin_half_at_zero_field_has_the_closed_form():
    """Every process is elastic, so the sample coupling runs as
    g(D) = g/(1 + g L(D)), L(D) = ln((D+omega0)/(D+Gamma0)), the on-shell
    tip coupling as g(D)/g, the exchange part of the tip vertex at eV as
    t = 1 - int g(D)^2/g K(eV,D) dD, U does not run, and
    dI/dV = 2 pi (3/4 t^2 + 4 U^2). With a sharp cut, t would be
    1/(1 - g ln(omega0/|eV|)), the leading log."""
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.zeros(2), p=np.array([.5, .5]), dim=2, T=0.,
                         kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    g, U = -0.05, 0.25
    eVs = np.array([-3e-3, -2e-4, 1e-5, 4e-4, 2e-3, 9e-3])
    got = leadinglog.leading_log_dIdV(ks, eVs, g, U=U, omega0=W0, Gamma0=G0)

    def gD(D): return g/(1. + g*np.log((D + W0)/(D + G0)))

    def t(v):
        f = lambda s: gD(np.exp(s))**2/g*_K(v, np.exp(s))*np.exp(s)
        val = quad(f, np.log(1e-12), np.log(1e6), limit=400, epsabs=1e-13,
                   points=[np.log(abs(v) + 1e-12), np.log(W0)])[0]
        return 1. - val
    ref = np.array([2*np.pi*(0.75*t(v)**2 + 4*U**2) for v in eVs])
    got2, info = leadinglog.leading_log_dIdV(ks, eVs, g, U=U, omega0=W0,
                                             Gamma0=G0, return_info=True)
    assert info["D_stop"] is None
    # the coupling ran from |g| to |g/(1+g ln(omega0/Gamma0))|
    assert info["max_coupling"] == pytest.approx(abs(g/(1 + g*np.log(W0/G0))), rel=1e-6)
    assert np.max(np.abs(got - ref)) < 1e-7*np.max(ref)
    # the Kondo sign: Jrho_s < 0 makes the tip vertex grow
    assert np.all(got > 2*np.pi*(0.75 + 4*U**2))


def test_spin_half_in_a_field_follows_the_published_equations():
    """cond-mat/0301106 eq. 6 at omega = 0 and +-B/2 is closed in g_z(0)
    and g_perp(B/2): dg_z/dlnD = -2 g_perp^2 Theta_B and
    dg_perp/dlnD = -g_perp g_z (Theta_B + Theta_0), with
    g_theirs = -g_ours/2 (their tau is the Pauli matrix and their AF sign
    is positive). With the sharp Theta_x replaced by the same smooth
    kernel, D*K(x,D), the flowed sample vertex must follow it."""
    B = 1e-3
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.array([0., B]), p=np.array([1., 0.]), dim=2, T=0.,
                         kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    # eigenbasis of H = B Sz: state 0 is spin down, state 1 spin up
    P = np.array([[0., 1.], [1., 0.]])
    ks.Sx, ks.Sy, ks.Sz = [P @ x @ P for x in S]
    ks.e = np.array([0., B])
    g = -0.08
    Dmax, Dmin = 1e10*(W0 + B), 1e-10*G0
    sol, D_stop = leadinglog.sample_flow(ks, g, omega0=W0, Gamma0=G0,
                                         Dmax=Dmax, Dmin=Dmin)
    assert D_stop is None

    def rhs(s, y):
        D = np.exp(s); gz, gp = y
        return [-2*gp**2*D*_K(B, D), -gp*gz*D*(_K(B, D) + _K(0., D))]
    ref = solve_ivp(rhs, (np.log(Dmax), np.log(Dmin)), [-g/2, -g/2],
                    method="DOP853", rtol=1e-11, atol=1e-14, dense_output=True)
    for s in np.linspace(np.log(Dmax), np.log(Dmin), 40):
        G = sol.sol(s).view(complex).reshape(3, 4, 2, 2)[0]
        gz0 = G[3, 0, 0]/(ks.Sz[0, 0]/2); gz1 = G[3, 1, 1]/(ks.Sz[1, 1]/2)
        gpx = G[1, 0, 1]/(ks.Sx[0, 1]/2); gpy = G[2, 0, 1]/(ks.Sy[0, 1]/2)
        gz_t, gp_t = ref.sol(s)
        for x in (gz0, gz1):
            assert abs(-x/2 - gz_t) < 1e-8
        for x in (gpx, gpy):
            assert abs(-x/2 - gp_t) < 1e-8
        assert np.max(np.abs(G[0])) < 1e-10 # no potential scattering at omega=0


def _paper_spin_half(B, om, g0, Dmax, Dmin):
    """cond-mat/0301106 eq. 6 (= cond-mat/0408506 eqs. 14, 15) at any
    omega, with the smooth kernel: the right-hand side holds only the
    Fermi-level couplings z = g_z(0) and p = g_perp(B/2). Returns
    g_z+(om), g_z-(om), g_perp(om-B/2), g_perp(om+B/2)."""
    def rhs(s, y):
        D = np.exp(s); z, p = y[0], y[1]
        return [-2*p*p*D*_K(B, D), -p*z*D*(_K(B, D) + _K(0., D)),
                -2*p*p*D*_K(om + B, D), -2*p*p*D*_K(om - B, D),
                -p*z*D*(_K(om, D) + _K(om - B, D)),
                -p*z*D*(_K(om + B, D) + _K(om, D))]
    sol = solve_ivp(rhs, (np.log(Dmax), np.log(Dmin)), [g0]*6, method="DOP853",
                    rtol=1e-11, atol=1e-16)
    return sol.y[2:, -1]


def test_spin_half_tip_vertex_follows_the_published_equations_at_any_omega():
    """The tip vertex at omega != 0 against _paper_spin_half, in the
    paper's convention (spin up at -B/2, g_paper = -Jrho_s/2,
    g_z_sigma = 4 g_paper sigma (A^0 + sigma A^z) on the up diagonal,
    g_perp(omega -+ B/2) from A^x below and above the diagonal)."""
    B, g = 1e-3, -0.08
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.array([-B/2, B/2]), p=np.array([1., 0.]), dim=2,
                         T=0., kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    oms = np.array([0.3e-3, 0.7e-3, -0.99e-3, 2.5e-3])
    A, _ = leadinglog.tip_vertex(ks, oms, 0.0, g, omega0=W0, Gamma0=G0)
    g0 = -g/2
    for a, om in zip(A, oms):
        Dmax, Dmin = 1e10*(W0 + B + abs(om)), 1e-10*G0
        ref = _paper_spin_half(B, om, g0, Dmax, Dmin)
        got = [4*g0*(a[0, 0, 0] + a[3, 0, 0]), -4*g0*(a[0, 0, 0] - a[3, 0, 0]),
               g0*a[1, 1, 0]/(S[0][1, 0]/2), g0*a[1, 0, 1]/(S[0][0, 1]/2)]
        assert np.max(np.abs(np.array(got) - ref)) < 1e-8


def test_spin_half_tip_vertex_on_shell_is_the_sample_vertex():
    """For S=1/2 in a field G(0) stays Hermitian, and the tip vertex at
    omega = e_f - e_i is the sample vertex G_if^* / g, the identity the
    published scheme keeps. (For a non-Hermitian G(0), as in an
    anisotropic S=1 in a tilted field, it holds only to O(g^2): see the
    README.)"""
    B, g = 1e-3, -0.08
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.array([-B/2, B/2]), p=np.array([1., 0.]), dim=2,
                         T=0., kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    sol, _ = leadinglog.onshell_flow(ks, g, omega0=W0, Gamma0=G0)
    G = sol.y[:, -1].view(complex).reshape(3, 4, 2, 2)[0]
    for f, i in [(1, 0), (0, 1)]:
        A, _ = leadinglog.tip_vertex(ks, [ks.e[f] - ks.e[i]], 0.0, g, omega0=W0, Gamma0=G0)
        assert np.max(np.abs(A[0, 1:, f, i] - np.conj(G[1:, i, f])/g)) < 1e-8


def test_a_large_tip_coupling_is_flagged():
    """S=3/2, D=1 meV, B=0.3 meV, Jrho_s=-0.1: the sample coupling stays
    small, but at the threshold eV=B the tip vertex is several times its
    bare value. The result must say so."""
    S = _spin(1.5)
    H = 1e-3*S[2]@S[2] + 0.3e-3*S[2]
    ks = _ks(H, S)
    eVs = np.array([0.3e-3 + 1e-6, 2.5e-3])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        got, info = leadinglog.leading_log_dIdV(ks, eVs, -0.1, omega0=W0,
                                                Gamma0=G0, return_info=True)
    assert info["coupling"][0] >= leadinglog.WARN_COUPLING
    assert any("running coupling" in str(x.message) for x in w)


def test_gamma0_must_be_positive():
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.zeros(2), p=np.array([.5, .5]), dim=2, T=0.,
                         kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    with pytest.raises(ValueError):
        leadinglog.leading_log_dIdV(ks, [1e-3], -0.05, Gamma0=0.0)


@pytest.mark.parametrize("name", list(MODELS))
def test_first_order_is_the_third_order(name):
    """Against second + third order of conductance.py, with U on: the
    difference is O(g^2), and its odd part in g is O(g^3), so the
    first order is exactly the third order (sign of g included)."""
    H, S = MODELS[name]
    ks = _ks(H, S)
    U = 0.25

    def LL(g): return leadinglog.leading_log_dIdV(ks, EVS, g, U=U, omega0=W0, Gamma0=G0)

    def third(g):
        return (conductance.third_order_kondo_dIdV(ks, EVS, g, omega0=W0, Gamma0=G0)
                + conductance.third_order_potential_dIdV(ks, EVS, g, U, omega0=W0, Gamma0=G0))
    d2 = conductance.second_order_dIdV(ks, EVS, U=U)
    r = {g: LL(g) - d2 - third(g) for g in (-0.01, -0.005, 0.005, 0.01)}
    ratio = np.max(np.abs(r[-0.01]))/np.max(np.abs(r[-0.005]))
    assert 3.5 < ratio < 4.6 # O(g^2)
    odd1 = (r[0.01] - r[-0.01])/2; odd2 = (r[0.005] - r[-0.005])/2
    assert np.max(np.abs(odd1)) < 0.05*np.max(np.abs(third(0.01)))
    assert 6. < np.max(np.abs(odd1))/np.max(np.abs(odd2)) < 10. # O(g^3)


def test_the_flow_stops_at_the_kondo_scale():
    """S=1/2 at B=0 with T_K above Gamma0: the coupling reaches order one,
    the biases near zero are NaN with a warning, the far ones are finite."""
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.zeros(2), p=np.array([.5, .5]), dim=2, T=0.,
                         kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    eVs = np.array([-1e-2, -1e-4, 0., 1e-4, 1e-2])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        got, info = leadinglog.leading_log_dIdV(ks, eVs, -0.2, omega0=W0,
                                                Gamma0=G0, return_info=True)
    assert info["D_stop"] is not None and info["D_stop"] > G0
    assert info["max_coupling"] >= 1.
    assert any("order one" in str(x.message) for x in w)
    assert np.all(np.isnan(got[1:4]))
    assert np.all(np.isfinite(got[[0, 4]]))


def test_finite_temperature_is_refused():
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.zeros(2), p=np.array([.5, .5]), dim=2, T=1.,
                         kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    with pytest.raises(NotImplementedError):
        leadinglog.leading_log_dIdV(ks, [1e-3], -0.05)


def test_get_kondo_spectrum_order_LL():
    """order="LL" on a 3-site S=1 chain with anisotropy and a tilted field
    is the module's result on the chain's KondoSpectrum, is close to
    order=3 at small Jrho_s, and refuses T>0 and mode="DMRG"."""
    from dmrgpy import spinchain
    from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
    mub = 5.7883818060e-5
    sc = spinchain.Spin_Chain(["1", "1/2", "1"])
    h = -1e-3*sc.Sz[0]*sc.Sz[0] + 2*mub*(sc.Sx[0] + 2*sc.Sz[0])
    h = h + 0.5e-3*(sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1])
    h = h + 0.3e-3*(sc.Sx[1]*sc.Sx[2] + sc.Sy[1]*sc.Sy[2] + sc.Sz[1]*sc.Sz[2])
    sc.set_hamiltonian(h)
    eVs = np.linspace(-3e-3, 3e-3, 13) + 1.3e-5
    _, got = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, U=0.25, T=0.,
                                   order="LL")
    ks = KondoSpectrum(sc, 0, 0.)
    ref = leadinglog.leading_log_dIdV(ks, eVs, -0.05, U=0.25)
    assert np.max(np.abs(got - ref)) < 1e-12*np.max(np.abs(ref))
    _, ll = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.002, U=0.25, T=0., order="LL")
    _, o3 = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.002, U=0.25, T=0., order=3)
    _, o2 = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.002, U=0.25, T=0., order=2)
    assert np.max(np.abs(ll - o3)) < 0.05*np.max(np.abs(o3 - o2))
    with pytest.raises(NotImplementedError):
        sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=1.0, order="LL")
    with pytest.raises(NotImplementedError):
        sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=0., order="LL", mode="DMRG")
    with pytest.raises(ValueError):
        sc.get_kondo_spectrum(eVs, site=0, order=4)


# ---------------------------------------------------------------------
# the joint flow (2026-11-11)

def _pauli_bracket(G, A, Kdirect, Kexchange):
    """The bracket of the flow written with the structure constants of
    the Pauli matrices, sixteen products per channel: the definition."""
    C = leadinglog._C
    part = np.einsum('mac,ncb,cb->mnab', G, A, Kdirect)
    hole = np.einsum('mac,ac,ncb->mnab', A, Kexchange, G)
    return (np.einsum('mnl,mnab->lab', C, part)
            - np.einsum('nml,mnab->lab', C, hole))


def test_the_block_product_is_the_pauli_contraction():
    """A vertex is a 2d x 2d matrix, and the contraction with the Pauli
    structure constants is its matrix product: _rate, _onshell_rate and
    the column-cut _tip_bracket against the sixteen-term definition."""
    rng = np.random.default_rng(3)
    d = 5
    rnd = lambda *sh: rng.normal(size=sh) + 1j*rng.normal(size=sh)
    G, A, K = rnd(4, d, d), rnd(4, d, d), rng.normal(size=(d, d))
    ref = _pauli_bracket(G, A, K, K)
    assert np.max(np.abs(leadinglog._rate(G, A, K) - ref)) < 1e-13*np.max(np.abs(ref))
    # on shell the kernel sits on G: direct K[a,c], exchange K[c,b]
    ref = (np.einsum('mnl,mac,ac,ncb->lab', leadinglog._C, G, K, A)
           - np.einsum('nml,mac,ncb,cb->lab', leadinglog._C, A, G, K))
    assert np.max(np.abs(leadinglog._onshell_rate(G, A, K) - ref)) < 1e-13*np.max(np.abs(ref))
    assert np.max(np.abs(leadinglog._unblock(leadinglog._block(A)) - A)) < 1e-15
    # the tip bracket at several energies, two columns, both kernels
    e = np.sort(rng.normal(size=d)); w = rng.normal(size=4); cols = np.array([1, 3])
    x = w[:, None, None] - (e[:, None] - e[None, :])[None]
    for D, Dlow in [(0.3, -1.), (0.7, 0.02)]:
        if Dlow < 0: Kw = leadinglog._kernel(x, D, 0.5, 0.01)
        else: Kw = np.array([[[quad(lambda D: leadinglog._kernel(xx, D, 0.5, 0.01), Dlow, D,
                                    epsabs=1e-13, epsrel=1e-13)[0] for xx in r] for r in xk] for xk in x])
        ref = np.array([_pauli_bracket(G, A, k, k) for k in Kw])[:, :, :, cols]
        got = leadinglog._cut_columns(leadinglog._tip_bracket(
            leadinglog._block(G), leadinglog._block(A), cols, w, e, D, Dlow, 0.5, 0.01))
        assert np.max(np.abs(got - ref)) < 1e-11*np.max(np.abs(ref))


def test_the_tip_potential_scattering_does_not_flow():
    """The direct and the exchange term of a potential vertex cancel at
    every cutoff, so A(omega; U) = A(omega; 0) + U and the on-shell tip
    vertex is A_on(0) + U. The flow of the tip is solved at U=0 on this."""
    ks = _ks(*MODELS["S=1, D=-1 meV, B=(1,0,2) T"])
    d, U, g = ks.dim, 0.3, -0.05
    sol, _ = leadinglog.onshell_flow(ks, g, U=U, omega0=W0, Gamma0=G0)
    for s in (sol.t[-1], sol.t[len(sol.t)//2]):
        X = sol.sol(s).view(complex).reshape(3, 4, d, d)
        dX = X[1] - X[2]
        dX[0] -= 2*U*np.eye(d)
        assert np.max(np.abs(dX)) < 1e-9
    # and the full rate of a potential vertex is zero for any G and kernel
    rng = np.random.default_rng(5)
    G = rng.normal(size=(4, d, d)) + 1j*rng.normal(size=(4, d, d))
    A = np.zeros((4, d, d), complex); A[0] = np.eye(d)
    K = np.full((d, d), 0.7) # the potential vertex is cut at omega for every state
    assert np.max(np.abs(leadinglog._rate(G, A, K))) < 1e-13
    assert np.max(np.abs(leadinglog._onshell_rate(G, A, rng.normal(size=(d, d))))) < 1e-13
    # tip_vertex adds each U at the end, any mix of them
    oms = np.array([1e-3, -2e-3, 0.4e-3])
    A0, _ = leadinglog.tip_vertex(ks, oms, 0., g, omega0=W0, Gamma0=G0)
    AU, _ = leadinglog.tip_vertex(ks, oms, [U, -U, 0.1], g, omega0=W0, Gamma0=G0)
    dA = AU - A0
    for k, u in enumerate([U, -U, 0.1]): dA[k, 0] -= u*np.eye(d)
    assert np.max(np.abs(dA)) < 1e-13


def test_the_occupied_columns_are_columns_of_the_full_vertex():
    """dI/dV reads only the columns of the occupied states; cut there, the
    tip vertex is the same numbers as the full one, and a repeated energy
    is solved once."""
    ks = _ks(*MODELS["S=1/2 dimer, J=1 meV, Bz=5 T, tip on 0"])
    oms = np.array([1e-3, -1e-3, 0.3e-3, 1e-3])
    full, _ = leadinglog.tip_vertex(ks, oms, 0.25, -0.05, omega0=W0, Gamma0=G0)
    cols = np.array([2, 0])
    cut, _ = leadinglog.tip_vertex(ks, oms, 0.25, -0.05, omega0=W0, Gamma0=G0, cols=cols)
    assert cut.shape == (4, 4, ks.dim, 2)
    assert np.max(np.abs(cut - full[:, :, :, cols])) < 1e-9
    assert np.array_equal(full[0], full[3])


def test_below_the_stop_the_quadrature_is_the_closed_form():
    """When the sample coupling reaches pole_threshold the couplings are
    held, and the rest of the tip quadrature is the integral of the
    kernel from Dmin to D_stop. Against the two-pass route: the on-shell
    flow first, then the tip vertex integrated down to Dmin with the
    couplings frozen below D_stop."""
    ks = _ks(*MODELS["S=1, D=-1 meV, B=(1,0,2) T"])
    d, g = ks.dim, -0.2
    oms = np.array([2.2e-3, -0.7e-3, 0.05e-3])
    got, info = leadinglog.tip_vertex(ks, oms, 0., g, omega0=W0, Gamma0=G0)
    assert info["D_stop"] is not None and info["D_stop"] > G0
    Dmax, Dmin = leadinglog._energy_range(ks, W0, oms, G0)
    sol, D_stop = leadinglog.onshell_flow(ks, g, omega0=W0, Gamma0=G0, Dmax=Dmax, Dmin=Dmin)
    assert abs(D_stop - info["D_stop"]) < 1e-9*D_stop
    frozen = sol.y[:, -1].view(complex).reshape(3, 4, d, d)
    x = oms[:, None, None] - (ks.e[:, None] - ks.e[None, :])[None]

    def rhs(s, y):
        X = frozen if s <= sol.t[-1] else sol.sol(s).view(complex).reshape(3, 4, d, d)
        D = np.exp(s)
        return (-D*leadinglog._rate(X[0], X[1], _K(x, D))).reshape(-1).view(float)
    A0 = np.broadcast_to(leadinglog._bare_vertex(leadinglog._spin_ops(ks), 1., 0.), (3, 4, d, d))
    ref = solve_ivp(rhs, (np.log(Dmax), np.log(Dmin)), A0.reshape(-1).view(float),
                    method="DOP853", rtol=1e-11, atol=1e-14).y[:, -1].view(complex).reshape(3, 4, d, d)
    assert np.max(np.abs(got - ref)) < 1e-8*np.max(np.abs(ref))


# leading_log_dIdV(ks, EVS[::5], -0.05, U=0.25, rtol=1e-12) of the two-pass
# implementation of 2026-11-06 (commit e4c8b28), which solved the on-shell
# flow with dense output and then one flow per bias for the full d x d tip
# vertex; the joint flow must give the same numbers
_TWO_PASS = {
    "S=1/2, Bx=10 T": [7.423359537055, 8.076148419175, 3.736654877698,
                       7.811601326244, 7.321681752024],
    "S=1, D=-1 meV, B=(1,0,2) T": [17.132506928825, 18.870761553211, 10.194985157709,
                                   18.283389856232, 16.913574725051],
    "S=1/2 dimer, J=1 meV, Bz=5 T, tip on 0": [7.451776352556, 8.17758329172,
                                               1.568310403449, 7.816335162193,
                                               7.321404128422],
}


@pytest.mark.parametrize("name", list(MODELS))
def test_the_joint_flow_gives_the_two_pass_values(name):
    ks = _ks(*MODELS[name])
    got = leadinglog.leading_log_dIdV(ks, EVS[::5], -0.05, U=0.25)
    assert got == pytest.approx(_TWO_PASS[name], rel=1e-9)
