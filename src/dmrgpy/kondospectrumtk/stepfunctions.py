import functools
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline

# Numerical building blocks for the third-order STM/Kondo perturbation
# theory of Ternes, New J. Phys. 17 063016 (2015), arXiv:1505.04430.
#
# Theta(x) and F(eps,T) below are NOT literal transcriptions of the
# closed-form equations printed in the arXiv source for them (eq.
# "step-fkt"/eq. 11 and equ. "F_2"/eq. 22 in the arXiv v1 PDF), because
# those printed closed forms fail basic requirements the paper itself
# states:
#   - eq. 11, Theta(x) = (1+(x-1)e^x)/e^(2x), diverges as x->-inf instead
#     of saturating at 0, so it cannot be the bounded, symmetric,
#     temperature-broadened step of the paper's own Fig. 4. Theta(x) here
#     was instead re-derived directly from eq. 9 (the paper's own
#     prescription: differentiate the current w.r.t. eV) and verified
#     against direct numerical integration to machine precision. (With
#     (e^x-1)^2 for the denominator the printed form is 1-Theta(x), the
#     textbook Lambe-Jaklevic step with the opposite sign convention for
#     x.)
#   - eq. 22 as printed, -int deps' ln[(w0+|eV-eps_m|)/(eV-eps_m+i*G0)]
#     Theta'(eV-eps_m+eps',T), has its log independent of the integration
#     variable eps' (which makes the integral trivial and drops exactly
#     the thermal broadening of the log singularity the surrounding text
#     and Fig. 5 describe) and an overall sign that makes it negative
#     where Fig. 5 is positive. The evident intent is the thermal
#     convolution of the closed-form log with the Theta' kernel,
#
#         F(x,T) = int deps' L0(x-eps') Theta'(eps'/kT)/kT,
#         L0(y)  = ln(w0+|y|) - 1/2 ln(y^2+G0^2),
#
#     which is what FBuilder/F below evaluate, and F0 = L0 is its exact
#     T=0 limit. Verified against the paper's own figures: the six peak
#     values of Fig. 5(b,c) (w0=200 meV, T=0.5..20 K) come out 8.12, 7.46,
#     6.78, 5.88, 5.19, 4.51 against 8.1, 7.5, 6.8, 5.9, 5.2, 4.5 read
#     off the plot, and the +-4 mV tails of the Fig. 7b spectra (w0=20
#     meV, T=1 K) come out 0.884 against 0.886 digitized from the figure
#     (tests/test_kondo_spectrum.py pins both).
#
# Until 2026-09-12 F was instead the paper's defining double integral for
# the ELECTRON-like process, eq. 20, evaluated exactly -- and used for the
# exchange diagram too, which the paper says is the hole-like eq. 21.
# The two are F_e(x) = ln|w0-x| - ln|x| and -F_h(x) = F_e(-x) at T=0, and
# eq. 22 is what both reduce to for |x| << w0; eq. 20's own thermal
# broadening turns out to be the same Theta' kernel (its outer f'
# convolution over its inner (1-f) cutoff), so at the Kondo peak nothing
# changes (7.459 -> 7.460 at 1 K/200 meV). What changes is the
# O(|x|/w0) band-edge term: F_e has the NEAR band edge for x>0, which put
# the Fig. 7b tails at 0.854 instead of 0.886, and a sharp-cutoff log
# singularity at x=w0 (a dip to ~-5, smeared only over kT) that eq. 22
# does not have -- a real hazard for a chain whose excitation energies
# approach the default w0=20 meV. Third-order numbers from before that
# date therefore differ from the current ones by up to ~3.5% of the total
# at |eV|=w0/5 (nothing at eV=0; the second-order term is untouched).


def Theta(x):
    """Temperature-broadened step function.

    Derived by differentiating eq. "current" w.r.t. eV (the paper's own
    prescription) and evaluating the resulting Fermi-function convolution
    in closed form:

        Theta(x) = 1/2 + [sinh(x) - x] / [2(cosh(x) - 1)]

    with x = eps/(kB T). Theta(x) -> 0 as x -> -inf, -> 1 as x -> +inf,
    Theta(0) = 1/2, and Theta(x) + Theta(-x) = 1 identically.

    Below |x| = 1e-4 it is the series 1/2 + x/6 - x^3/180 (next term
    x^5/5040), and above |x| = 40 the tail (|x|-1) e^-|x| of Theta(-|x|),
    both exact to double precision. The closed form between loses about
    1e-16/|x| to the cancellation in sinh(x) - x and cosh(x) - 1 (2e-13 at
    |x| = 2e-4), which is left as it is.
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    large = np.abs(x) > 40
    mid = ~small & ~large
    out[small] = 0.5 + x[small]/6.0 - x[small]**3/180.0
    xm = x[mid]
    out[mid] = 0.5 + (np.sinh(xm) - xm)/(2*(np.cosh(xm) - 1))
    xl = x[large]
    pos = xl > 0
    out[large] = np.where(pos,
                           1.0 - (np.abs(xl)-1)*np.exp(-np.abs(xl)),
                           (np.abs(xl)-1)*np.exp(-np.abs(xl)))
    return out


def Theta_prime(x):
    """d Theta / dx, the even, unit-normalized (integral over x is 1)
    thermal broadening kernel used inside F(eps,T).

    Closed form (x cosh(x/2) - 2 sinh(x/2)) / (4 sinh(x/2)^3), whose
    numerator is a difference of two nearly equal numbers for small x
    (relative error about 12 eps/x^2), so below |x|=0.05 the Taylor series
    1/6 - x^2/60 + x^4/1008 - x^6/21600 + x^8/532224 is used instead (the
    next term is 1e-13 at the switch). Both pieces keep 1e-13 or better."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 0.05
    large = np.abs(x) > 40
    mid = ~small & ~large
    xs = x[small]**2
    out[small] = 1.0/6 + xs*(-1.0/60 + xs*(1.0/1008 + xs*(-1.0/21600
                                                       + xs/532224.)))
    xm = x[mid]
    h = 0.5*xm
    out[mid] = (xm*np.cosh(h) - 2*np.sinh(h))/(4*np.sinh(h)**3)
    xl = x[large]
    out[large] = np.abs(xl)*np.exp(-np.abs(xl))
    return out


def Theta0(x):
    """T=0 limit of Theta(x): a Heaviside step, Theta0(0)=1/2. Both the
    spin system's Boltzmann occupation and the tunneling-electron thermal
    smearing collapse at T=0, so this is not merely "Theta at very small
    T" evaluated numerically -- it is the exact, closed-form limit."""
    x = np.asarray(x, dtype=float)
    return np.where(x > 0, 1., np.where(x < 0, 0., 0.5))


def F0(x, omega0=20e-3, Gamma0=5e-6):
    """T=0 limit of F(eps,T): the paper's closed-form log, eq. 22,

        F0(eps) = ln(omega0+|eps|) - 1/2 ln(eps^2+Gamma0^2),

    i.e. Re ln[(omega0+|eps|)/(eps+i*Gamma0)] -- even in eps, a
    Gamma0-regularized log peak at eps=0 decaying as ln(omega0/|eps|)
    (~1/|eps| beyond the band) with no feature at the band edge. See the
    module docstring for how this relates to the electron/hole-like
    defining integrals eqs. 20/21, and for what it replaced."""
    x = np.asarray(x, dtype=float)
    return np.log(omega0 + np.abs(x)) - 0.5*np.log(x**2 + Gamma0**2)


def _F_convolve(x, kT, omega0, Gamma0, npts=8001, width=60.):
    """F(x,T) = int du Theta'(u/kT)/kT * F0(x+u), vectorized over x.

    The -1/2 ln((x+u)^2+G0^2) part of F0 is log-singular at u=-x on the
    Gamma0 scale, far below any grid a kT-wide kernel wants, so when that
    point lies inside the kernel window it is subtracted exactly (the
    kernel value there times the closed-form integral of the log) and a
    plain Simpson rule only ever sees the bounded remainder
    (K(u)-K(-x))*L(x+u). The ln(w0+|x+u|) part has a mere kink there."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    W = width*kT
    u = np.linspace(-W, W, npts)
    K = Theta_prime(u/kT)/kT
    def L(t): return 0.5*np.log(t**2 + Gamma0**2)
    def IL(t): # antiderivative of L
        return 0.5*(t*np.log(t**2 + Gamma0**2) - 2*t + 2*Gamma0*np.arctan(t/Gamma0))
    inside = np.abs(x) < W
    K_at = np.where(inside, Theta_prime(-x/kT)/kT, 0.) # kernel at u=-x
    out = np.empty_like(x)
    chunk = max(1, int(2e7)//npts)
    for s0 in range(0, len(x), chunk):
        xc = x[s0:s0+chunk]
        y = xc[:, None] + u[None, :]
        integrand = (K[None, :]*np.log(omega0 + np.abs(y))
                     - (K[None, :] - K_at[s0:s0+chunk, None])*L(y))
        out[s0:s0+chunk] = simpson(integrand, x=u, axis=1)
    out -= K_at*(IL(x + W) - IL(x - W)) # int_{-W}^{W} L(x+u) du
    return out


def _merge_grid(points, min_gap):
    """Sorted union of grid points with near-duplicates removed: two
    nodes closer than min_gap make a cubic spline ill-conditioned (it
    oscillates between them), and np.unique only drops exact repeats."""
    pts = np.sort(np.asarray(points, dtype=float))
    keep = np.ones(len(pts), dtype=bool)
    last = pts[0]
    for i in range(1, len(pts)):
        if pts[i] - last < min_gap: keep[i] = False
        else: last = pts[i]
    return pts[keep]


class FBuilder():
    """Evaluator of the temperature-broadened Kondo log F(eps,T), the
    paper's eq. 22 read as the Theta'-convolution of its closed-form log
    (see the module docstring): F = F0 (*) Theta'(./kT)/kT.

    F is tabulated once at construction (dense at kT/20 over +-40kT,
    where it has its kT-wide peak, geometric out to +-grid_span*omega0
    where it is a slow log) and spline-interpolated by __call__; beyond
    the table the thermal broadening is negligible and F0 itself is
    returned. Building one costs ~0.4 s and evaluating it is O(1) per
    point (1e6 points in 0.15 s), so build it once per (T, omega0,
    Gamma0) and share it across the third-order terms
    (Spin_Chain.get_kondo_spectrum does; cached_fbuilder shares one
    across calls). Accuracy, for omega0 >~ kT (the regime of the paper,
    1 to 20 K with omega0 of 20 to 200 meV): the table reproduces the
    direct convolution to ~1e-8, and the convolution itself is converged
    in npts to ~1e-6 (npts=4001; 8001 buys 1e-6 at 20 K and nothing at
    1 K). For omega0 << kT the u step of 0.03 kT does not resolve the kink
    of ln(omega0+|y|), and both the table and the direct convolution are
    off by about 2e-4 of F (omega0=2 meV at 300 K). Until
    2026-09-12 __call__ ran a kT-wide quadrature per requested point on
    an unchunked (npoints x 2001) array -- ~5 GB and minutes for 21 bias
    points on 128 states; this is what made the T>0 third-order terms
    F-bound.

    tabulate=False skips the table and has __call__ run the convolution
    directly (_direct), for validating the spline against it.

    odd_part is max |F(x) - F(-x)| over the nodes of the table and the
    midpoints between them, relative to max |F| there: F is even, and the
    table is even to roundoff (about 1e-14), which conductance.py checks
    before it reads an s->t table off a t->s one (EVEN_TOL). It is None
    without a table, and the exact form is kept then."""
    def __init__(self, T, omega0=20e-3, Gamma0=5e-6, kB=8.617333262e-5,
                 npts=4001, width=60., grid_span=10., tabulate=True):
        if T <= 0.: raise ValueError("F(eps,T) requires T>0")
        self.kT = kB*T
        self.omega0 = omega0
        self.Gamma0 = Gamma0
        self._npts, self._width = npts, width
        self._F_spline = None
        self.odd_part = None
        if tabulate:
            kT = self.kT
            dense = np.linspace(-40*kT, 40*kT, 1601)
            top = grid_span*omega0
            parts = [dense]
            if top > 40*kT:
                outer = np.geomspace(40*kT, top, 400)
                parts += [outer, -outer]
            grid = _merge_grid(np.concatenate(parts), 0.01*kT)
            self._F_spline = CubicSpline(grid, self._direct(grid),
                                         extrapolate=False)
            self._F_min, self._F_max = grid[0], grid[-1]
            x = np.concatenate([grid, 0.5*(grid[1:] + grid[:-1])])
            Fx = self(x)
            self.odd_part = float(np.max(np.abs(Fx - self(-x)))/np.max(np.abs(Fx)))
    def _direct(self, x):
        return _F_convolve(x, self.kT, self.omega0, self.Gamma0,
                           npts=self._npts, width=self._width)
    def __call__(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        if self._F_spline is None: return self._direct(x)
        inside = (x >= self._F_min) & (x <= self._F_max)
        if inside.all(): return self._F_spline(x)
        out = np.empty_like(x)
        out[inside] = self._F_spline(x[inside])
        xo = x[~inside]
        out[~inside] = F0(xo, omega0=self.omega0, Gamma0=self.Gamma0)
        return out


@functools.lru_cache(maxsize=4)
def _cached_fbuilder(T, omega0, Gamma0, kB, npts, width, grid_span, tabulate):
    return FBuilder(T, omega0=omega0, Gamma0=Gamma0, kB=kB, npts=npts,
                    width=width, grid_span=grid_span, tabulate=tabulate)


def cached_fbuilder(T, omega0=20e-3, Gamma0=5e-6, kB=8.617333262e-5,
                    npts=4001, width=60., grid_span=10., tabulate=True):
    """An FBuilder shared between calls: the last four distinct argument
    sets are kept, keyed on the float values, so a sweep at fixed T,
    omega0, Gamma0 and kB builds the table once (about 0.2 to 0.6 s).
    The returned object is shared, so it must not be mutated."""
    return _cached_fbuilder(float(T), float(omega0), float(Gamma0),
                            float(kB), int(npts), float(width),
                            float(grid_span), bool(tabulate))


def F(x, T, omega0=20e-3, Gamma0=5e-6, kB=8.617333262e-5):
    """Convenience one-shot evaluation of F(eps,T) at eps=x (see FBuilder
    for the tabulated version to use when sweeping many eV points
    against a fixed set of intermediate-state energies)."""
    return FBuilder(T, omega0=omega0, Gamma0=Gamma0, kB=kB)(x)
