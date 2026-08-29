"""High-order delta-Chebyshev (HODC) KPM reconstruction, arXiv:2512.03149.

Everything except the last test works on *exact* Chebyshev moments of a
known density, so the kernel is pinned independently of any DMRG backend
-- which is the point: HODC is pure post-processing of the same moments
ordinary KPM starts from (see algebra/kpm.py's block comment)."""
import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy.algebra import kpm


# a smooth test density on [-1,1], narrow enough that its value at the
# endpoints (~3e-4) does not pollute the interior convergence rates
_SIG = 0.25
def _rho(x): return np.exp(-x**2/(2*_SIG**2))/(np.sqrt(2*np.pi)*_SIG)

_GAUSS_CACHE = {}
def _gauss(ngauss):
    """Gauss-Legendre nodes/weights, cached: leggauss is O(n^2) and this
    file asks for the same rule a dozen times."""
    if ngauss not in _GAUSS_CACHE:
        _GAUSS_CACHE[ngauss] = np.polynomial.legendre.leggauss(ngauss)
    return _GAUSS_CACHE[ngauss]


_MOMENT_CACHE = {}
def _exact_moments(npol,ngauss=2000):
    """mu_k = int rho(x) T_k(x) dx, the convention generate_profile uses.
    Cached, and built as one matmul rather than a per-k Gauss sum: the
    naive loop dominated this file's runtime."""
    key = (npol,ngauss)
    if key not in _MOMENT_CACHE:
        x,w = _gauss(ngauss)
        th = np.arccos(x)
        _MOMENT_CACHE[key] = np.cos(np.outer(np.arange(npol),th))@(w*_rho(x))
    return _MOMENT_CACHE[key]

def _kernel(E,x,order,eta):
    """K_eta(E,x) straight from Eq. (15), for cross-checking."""
    z,w = kpm.hodc_poles_weights(order)
    return -np.sum((w[:,None]/(E-np.asarray(x)[None,:]+eta*z[:,None])).imag,
                   axis=0)/np.pi


def test_poles_satisfy_the_moment_matching_conditions():
    """sum_l w_l z_l^j = delta_{j0} for j < m -- the condition that makes
    the regularized delta O(eta^m) accurate."""
    for order in range(1,kpm.HODC_MAX_ORDER+1):
        z,w = kpm.hodc_poles_weights(order)
        assert np.all(z.imag>0) # every pole in the upper half plane
        got = np.array([np.sum(w*z**j) for j in range(order)])
        want = np.zeros(order,dtype=complex) ; want[0] = 1.0
        assert np.allclose(got,want,atol=1e-10)
        # symmetric pole set -> conjugate-symmetric weights -> even kernel
        assert np.allclose(w,np.conjugate(w[::-1]),atol=1e-10)


def test_rejects_out_of_range_order():
    with pytest.raises(ValueError): kpm.hodc_poles_weights(0)
    with pytest.raises(ValueError): kpm.hodc_poles_weights(kpm.HODC_MAX_ORDER+1)
    with pytest.raises(ValueError):
        kpm.generate_profile_hodc(np.ones(10),np.array([0.0]),eta=-1.0)


def test_coefficients_match_chebyshev_quadrature():
    """The closed form for nu_k replaces the paper's fast cosine
    transform; check it against brute-force Chebyshev-Gauss quadrature of
    K_eta(E,.) itself."""
    npol,ngrid = 30,4000
    th = (np.arange(ngrid)+0.5)*np.pi/ngrid
    xq = np.cos(th)
    for order in (1,2,4,6):
        eta = 0.13
        for E in (-0.6,0.0,0.37):
            f = _kernel(E,xq,order,eta)
            want = np.array([(2-(k==0))/ngrid*np.sum(f*np.cos(k*th))
                             for k in range(npol)])
            got = kpm.hodc_coefficients(np.array([E]),npol,order=order,
                                        eta=eta)[:,0]
            assert np.allclose(got,want,atol=1e-12)


def test_order_one_is_exactly_lorentzian_broadening():
    """m=1 is w=1, z=i: the regularized delta collapses to a Lorentzian,
    so the reconstruction must reproduce an explicit convolution. This
    pins the sign, the branch of sqrt(zeta^2-1) and the eta convention
    all at once."""
    x,w = _gauss(2000)
    for eta in (0.2,0.1):
        npol = int(40/eta)
        mus = _exact_moments(npol)
        for E in (-0.3,0.0,0.2):
            got = kpm.generate_profile_hodc(mus,np.array([E]),order=1,eta=eta)
            want = np.sum(w*_rho(x)*(eta/np.pi)/((E-x)**2+eta**2))
            assert got[0].real == pytest.approx(want,abs=1e-10)


@pytest.mark.parametrize("order,expected_rate",[(1,1.0),(2,2.0),(4,4.0)])
def test_convergence_is_order_m_in_eta(order,expected_rate):
    """Halving eta must cut the error at a smooth point by 2^m. npol is
    kept at 45/eta so the Chebyshev truncation (which costs ~exp(-npol*
    eta)) stays far below the O(eta^m) regularization error."""
    E = 0.2
    errs = []
    for eta in (0.32,0.16,0.08):
        npol = int(45/eta)
        got = kpm.generate_profile_hodc(_exact_moments(npol),np.array([E]),
                                        order=order,eta=eta)
        errs.append(abs(got[0].real-_rho(E)))
    rates = [np.log2(errs[i]/errs[i+1]) for i in range(len(errs)-1)]
    assert min(rates) > expected_rate-0.5


def test_beats_jackson_on_a_smooth_density_at_equal_moment_count():
    """The whole point: same moments, same cost, smaller error. Both
    kernels are given their own natural width -- pi/npol for Jackson,
    hodc_default_eta(npol) for HODC."""
    npol = 400
    xs = np.linspace(-0.75,0.75,300)
    mus = _exact_moments(npol)
    ref = _rho(xs)
    ej = np.linalg.norm(kpm.generate_profile(mus,xs).real-ref)
    eh = np.linalg.norm(kpm.generate_profile(mus,xs,kernel="hodc").real-ref)
    assert eh < ej/10.0


def test_complex_moments_decouple():
    """nu_k(E,eta) is real, so a complex moment array must reconstruct as
    the real reconstruction of its real part plus i times that of its
    imaginary part -- the property that lets A != B correlators use the
    same code path as jackson."""
    rng = np.random.default_rng(0)
    a,b = rng.normal(size=64),rng.normal(size=64)
    xs = np.linspace(-0.8,0.8,17)
    got = kpm.generate_profile(a+1j*b,xs,kernel="hodc")
    ra = kpm.generate_profile(a,xs,kernel="hodc")
    rb = kpm.generate_profile(b,xs,kernel="hodc")
    assert np.allclose(got,ra.real+1j*rb.real,atol=1e-12)
    # generate_profile_pair must agree with two scalar calls
    pr,pi = kpm.generate_profile_pair(a,b,xs,kernel="hodc")
    assert np.allclose(pr,ra,atol=1e-12) and np.allclose(pi,rb,atol=1e-12)


def test_kpmdmrg_reconstruction_units():
    """End-to-end through kpmdmrg.dynamical_correlator_from_moments: a
    synthetic stick spectrum, rescaled exactly the way the backends
    rescale a Hamiltonian, must come back broadened by K_eta at the
    requested delta in *physical* units."""
    from dmrgpy import kpmdmrg
    emin,emax,kpm_scale = -8.0,8.0,0.7
    scale = 1.0/((emax-emin)*kpm_scale)
    om = np.array([0.7,1.3,2.6,4.1])       # excitation energies above emin
    wt = np.array([0.3,0.25,0.3,0.15])
    xj = (emin+om)*scale                   # where they sit in [-1,1]
    npol,delta,order = 600,0.25,6
    mus = np.array([np.sum(wt*np.cos(k*np.arccos(xj))) for k in range(npol)])
    es = np.linspace(0.3,4.5,120)
    want = np.sum(wt[None,:]*np.array([_kernel(e,om,order,delta)
                                       for e in es]),axis=1)
    # via the delta default...
    _,ys = kpmdmrg.dynamical_correlator_from_moments(mus,emin,emax,scale,
            npol,es,kernel="hodc",delta=delta,hodc_order=order)
    assert np.allclose(ys.real,want,atol=2e-3)
    assert np.max(np.abs(ys.imag)) < 1e-9
    # ...and stating eta explicitly, in the same physical units as delta
    # and es (not the rescaled ones generate_profile works in)
    _,ys2 = kpmdmrg.dynamical_correlator_from_moments(mus,emin,emax,scale,
            npol,es,kernel="hodc",delta=0.9,hodc_eta=delta,hodc_order=order)
    assert np.allclose(ys2.real,want,atol=2e-3)


def test_dmrg_spectrum_carries_the_full_sum_rule():
    """A real (small) KPM-DMRG run on the pure-Python backend: whatever
    the kernel, integrating <Sz delta(w-H+E0) Sz> over w must return
    <Sz^2> = 1/4, and HODC must land on the same spectrum as Jackson to
    within the two kernels' own widths."""
    from dmrgpy import spinchain
    ns = 8
    sc = spinchain.Spin_Chain(["S=1/2"]*ns)
    sc.itensor_version = "python" ; sc.setup_python()
    h = 0
    for i in range(ns-1): h = h + sc.SS(i,i+1)
    sc.set_hamiltonian(h)
    sc.maxm,sc.kpmmaxm = 30,30
    es = np.linspace(-0.5,6.0,400)
    name = (sc.Sz[ns//2],sc.Sz[ns//2])
    out = {}
    for kern in ("jackson","hodc"):
        _,y = sc.get_dynamical_correlator(name=name,es=es,delta=0.2,
                                          kernel=kern)
        out[kern] = y.real
        assert np.trapezoid(y.real,es) == pytest.approx(0.25,abs=5e-3)
    # same physics, different smoothing: compare after a common blur
    blur = lambda y: np.convolve(y,np.ones(25)/25.,mode="same")
    assert np.max(np.abs(blur(out["jackson"])-blur(out["hodc"]))) < 0.05


@pytest.mark.parametrize("version",[
    pytest.param(2,marks=pytest.mark.skipif(
        not cppext.available(2),reason="mpscpp2 not compiled")),
    pytest.param(3,marks=pytest.mark.skipif(
        not cppext.available(3),reason="mpscpp3 not compiled")),
])
def test_same_hodc_spectrum_on_every_backend(version):
    """HODC only ever touches the moments, so the compiled backends must
    return the same spectrum as the pure-Python one -- this is what makes
    "kernel= is available on every backend" in the user guide a fact
    rather than an intention."""
    from dmrgpy import spinchain
    ns = 10
    es = np.linspace(0.0,5.0,300)
    def spectrum(itensor_version):
        sc = spinchain.Spin_Chain(["S=1/2"]*ns)
        sc.itensor_version = itensor_version
        if itensor_version == "python": sc.setup_python()
        else: sc.setup_cpp(version=itensor_version)
        h = 0
        for i in range(ns-1): h = h + sc.SS(i,i+1)
        sc.set_hamiltonian(h)
        sc.maxm,sc.kpmmaxm = 40,40
        _,y = sc.get_dynamical_correlator(name=(sc.Sz[ns//2],sc.Sz[ns//2]),
                                          es=es,delta=0.3,kernel="hodc")
        return y.real
    assert np.max(np.abs(spectrum(version)-spectrum("python"))) < 1e-3
