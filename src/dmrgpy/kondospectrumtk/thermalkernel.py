"""The cutoff kernel of the leading-log flow at a finite temperature.

At T=0 the flow of leadinglog.py integrates the kernel
K(x,D) = 1/(D+sqrt(x^2+Gamma0^2)) - 1/(D+omega0+|x|), whose integral over
the cutoff D is F0(x) of stepfunctions.py. At a temperature T the third
order of conductance.py uses the thermal log F(x,T), the convolution of F0
with the weight w(u) = Theta'(u/kT)/kT of stepfunctions.FBuilder. So the
thermal flow uses the thermal kernel

    K_T(x,D) = int du w(u) K(x+u,D),   int_0^inf dD K_T(x,D) = F(x,T),

and its first order in Jrho_s is exactly the thermal third order. The
derivation and the checks are in docs/new_methods/kondo_leading_log/
README.md, section "Finite temperature". The construction was developed
in a project on 1T-TaSe2 (2026-11-11) on top of dmrgpy, and moved here.

How it is computed. At each cutoff D the even function x -> K_T(x,D) is
tabulated once and interpolated with cubic splines:
  - near region, |x| < xnear kT: K_T itself on a uniform grid of step
    kT/hdiv. It is smooth on the scale kT there, while K has a peak of
    width max(D, Gamma0) at x=0. The integral is done in y = x+u on
    Gauss-Legendre panels that halve towards y=0 down to Gamma0/2, so the
    peak is resolved for every x and every D. The quadrature is a fixed
    matrix, so a new cutoff costs one matrix-vector product;
  - far region, up to xmax: the relative correction K_T/K - 1, of order
    3.3 (kT/x)^2, on a grid uniform in ln|x|, from a Gauss-Legendre rule
    in u. Beyond xmax the plain K is returned.
A fixed short quadrature in u would not do: for D and Gamma0 below kT it
misses the peak at u=-x. One table takes about 2 ms.
"""
import numpy as np
from numba import njit
from scipy.interpolate import CubicSpline

from .stepfunctions import Theta_prime

# the kernel rows of the exchange sum (hole_sum) are filled in blocks of
# final states that take at most this many bytes
ROWS_BYTES = 256e6


def _k0(ax, D, Dlow, omega0, Gamma0):
    """K(x,D) on one denominator (no cancellation at large D), or its
    integral over the cutoff from Dlow to D when Dlow >= 0; ax = |x|. The
    same functions as leadinglog._kernel and leadinglog._kernel_integral."""
    r = np.sqrt(ax*ax + Gamma0*Gamma0)
    if Dlow < 0.:
        return (omega0 + ax - r)/((D + r)*(D + omega0 + ax))
    return np.log((D + r)/(Dlow + r)) - np.log((D + omega0 + ax)/(Dlow + omega0 + ax))


@njit(cache=True, fastmath=True, inline="always")
def _kval(ax, D, Dlow, omega0, G2, h, nn, cn, t0, iht, nf, cf):
    """K_T at |x| = ax from the tables of one cutoff: a cubic spline of K_T
    itself below nn*h (the near region), and K times (1 + a spline in
    ln|x|) above."""
    if ax < h*nn:
        s = ax/h
        j = int(s)
        if j >= nn:
            j = nn - 1
        ds = s - j
        return ((cn[j, 0]*ds + cn[j, 1])*ds + cn[j, 2])*ds + cn[j, 3]
    r = np.sqrt(ax*ax + G2)
    if Dlow < 0.:
        k = (omega0 + ax - r)/((D + r)*(D + omega0 + ax))
    else:
        k = np.log((D + r)/(Dlow + r)) - np.log((D + omega0 + ax)/(Dlow + omega0 + ax))
    s = (np.log(ax) - t0)*iht
    j = int(s)
    if j >= nf:
        return k # beyond the table the thermal correction is below 1e-8
    ds = s - j
    return k*(1. + ((cf[j, 0]*ds + cf[j, 1])*ds + cf[j, 2])*ds + cf[j, 3])


@njit(cache=True, fastmath=True)
def _eval_flat(ax, D, Dlow, omega0, G2, h, nn, cn, t0, iht, nf, cf, out):
    for i in range(len(ax)):
        out[i] = _kval(ax[i], D, Dlow, omega0, G2, h, nn, cn, t0, iht, nf, cf)


@njit(cache=True, fastmath=True)
def _fill_rows(w, e, a0, a1, D, Dlow, omega0, G2, h, nn, cn, t0, iht, nf, cf, Kr):
    """Kr[a-a0, k, c] = K_T(w[k] - (e[a] - e[c])) for a in [a0, a1). One
    core, like the rest of the flow."""
    nk, d = len(w), len(e)
    for ia in range(a1 - a0):
        ea = e[a0 + ia]
        for k in range(nk):
            u = w[k] - ea
            for c in range(d):
                Kr[ia, k, c] = _kval(abs(u + e[c]), D, Dlow, omega0, G2, h, nn, cn,
                                     t0, iht, nf, cf)


def _gl_panels(edges, ngl):
    """Composite Gauss-Legendre nodes and weights on the panels between
    consecutive edges."""
    t, wt = np.polynomial.legendre.leggauss(ngl)
    a, b = edges[:-1, None], edges[1:, None]
    return (0.5*(a + b) + 0.5*(b - a)*t[None]).ravel(), (0.5*(b - a)*wt[None]).ravel()


class ThermalKernel():
    """K_T(x,D) = int du w(u) K(x+u,D), w(u) = Theta'(u/kT)/kT, the weight
    of stepfunctions.FBuilder, tabulated per cutoff D (see the module
    docstring). kT, omega0, Gamma0 and xmax in the units of the spectrum;
    xmax is the largest |x| the flow asks for. The keywords set the tables:
    the near region is xnear*kT wide with hdiv points per kT, the far region
    has per_efold points per e-fold, the quadratures ngl Gauss-Legendre
    points per panel and panels of upanel*kT in u, over +-width*kT. The
    defaults put the integral over D within 3e-10 of F(x,T)."""

    def __init__(self, kT, omega0, Gamma0, xmax, width=60., hdiv=32, xnear=48.,
                 per_efold=32, ngl=16, upanel=2.):
        if not kT > 0.:
            raise ValueError("ThermalKernel needs kT > 0")
        if not Gamma0 > 0.:
            raise ValueError("Gamma0 must be positive, it cuts the logs")
        self.kT, self.omega0, self.Gamma0 = float(kT), float(omega0), float(Gamma0)
        self.h = self.kT/hdiv
        self.nn = int(round(xnear*hdiv))
        xn = self.nn*self.h
        xj = self.h*np.arange(self.nn + 1)
        # near region: nodes y >= 0 (K is even), panels [0,G/2], [G/2,G],
        # doubling up to kT, then of width kT
        edges, ytop = [0.], xn + width*self.kT
        while edges[-1] < ytop:
            edges.append(edges[-1] + min(max(edges[-1], 0.5*self.Gamma0), self.kT))
        self.y, gw = _gl_panels(np.array(edges), ngl)
        wfun = lambda u: Theta_prime(u/self.kT)/self.kT
        self.M = gw[None, :]*(wfun(xj[:, None] - self.y[None, :])
                              + wfun(xj[:, None] + self.y[None, :]))
        # far region: uniform in ln|x| from xn to xmax, Gauss-Legendre in u
        self.t0 = np.log(xn)
        self.nf = max(4, int(np.ceil(per_efold*(np.log(max(xmax, 2*xn)) - self.t0))) + 1)
        self.iht = float(per_efold)
        self.xf = np.exp(self.t0 + np.arange(self.nf + 1)/per_efold)
        self.u, uw = _gl_panels(self.kT*np.arange(-width, width + 0.5*upanel, upanel), 12)
        self.uw = uw*wfun(self.u)
        self._cache = {}
        self.nfev = 0 # tables built

    def _tables(self, f):
        vn = self.M @ f(self.y)
        cn = CubicSpline(np.arange(self.nn + 1), vn, bc_type=((1, 0.), "not-a-knot")).c
        fx = f(self.xf)
        rel = ((f(np.abs(self.xf[:, None] + self.u[None, :])) - fx[:, None]) @ self.uw)/fx
        cf = CubicSpline(np.arange(self.nf + 1), rel).c
        return np.ascontiguousarray(cn.T), np.ascontiguousarray(cf.T)

    def table(self, D, Dlow=-1.):
        """The spline coefficients of K_T at the cutoff D, or of its
        integral over the cutoff from Dlow to D when Dlow >= 0. The flow
        asks for the same D several times in one step, so the last few are
        kept."""
        key = (float(D), float(Dlow))
        if key not in self._cache:
            if len(self._cache) > 4:
                self._cache.clear()
            self._cache[key] = self._tables(
                lambda y: _k0(y, D, Dlow, self.omega0, self.Gamma0))
            self.nfev += 1
        return self._cache[key]

    def _args(self, D, Dlow):
        cn, cf = self.table(D, Dlow)
        return (float(D), float(Dlow), self.omega0, self.Gamma0**2, self.h, self.nn,
                cn, self.t0, self.iht, self.nf, cf)

    def __call__(self, x, D, Dlow=-1.):
        """K_T(x,D), or its integral over the cutoff from Dlow to D when
        Dlow >= 0, at every element of x."""
        x = np.asarray(x, dtype=float)
        out = np.empty(x.size)
        _eval_flat(np.ascontiguousarray(np.abs(x).ravel()), *self._args(D, Dlow), out)
        return out.reshape(x.shape)

    def hole_sum(self, Q, w, e, D, Dlow, out):
        """out[k,a,l] = sum_c Q[a,l,c] K_T(w[k] - (e[a] - e[c])), the
        exchange sum of leadinglog._hole_sum with the thermal kernel, for a
        real Q. The kernel rows are filled per block of final states and
        contracted with a matrix product."""
        n, d = len(w), len(e)
        args = self._args(D, Dlow)
        na = max(1, min(d, int(ROWS_BYTES//(8*n*d))))
        Kr = np.empty((na, n, d))
        for a0 in range(0, d, na):
            a1 = min(d, a0 + na)
            _fill_rows(w, e, a0, a1, *args, Kr)
            out[:, a0:a1, :] = np.matmul(Kr[:a1 - a0],
                                         Q[a0:a1].transpose(0, 2, 1)).transpose(1, 0, 2)

    def average(self, f, x):
        """The thermal average int du w(u) f(x+u) of an even function f
        with the structure of K (a peak of width >= Gamma0 at 0, smooth on
        the scale |y| away from it), at the points x, through the same
        tables. Slow; for checks, such as F(x,T) from f = F0."""
        cn, cf = self._tables(f)
        ax = np.abs(np.atleast_1d(np.asarray(x, dtype=float)))
        out = np.array(f(ax), dtype=float)
        near = ax < self.h*self.nn
        s = ax[near]/self.h
        j = np.minimum(s.astype(int), self.nn - 1)
        ds = s - j
        out[near] = ((cn[j, 0]*ds + cn[j, 1])*ds + cn[j, 2])*ds + cn[j, 3]
        far = ~near
        s = (np.log(ax[far]) - self.t0)*self.iht
        j = s.astype(int)
        ok = j < self.nf
        j, ds = j[ok], (s - j)[ok]
        idx = np.nonzero(far)[0][ok]
        out[idx] *= 1. + ((cf[j, 0]*ds + cf[j, 1])*ds + cf[j, 2])*ds + cf[j, 3]
        return out
