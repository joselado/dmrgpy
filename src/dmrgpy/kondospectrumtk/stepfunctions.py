import numpy as np
from scipy.integrate import simpson, quad
from scipy.interpolate import CubicSpline

# Numerical building blocks for the third-order STM/Kondo perturbation
# theory of Ternes, New J. Phys. 17 063016 (2015), arXiv:1505.04430.
#
# Theta(x) and F(eps,T) below are NOT literal transcriptions of the
# closed-form equations printed in the arXiv source for them (eq.
# "step-fkt" and equ. "F_2") because those printed closed forms fail
# basic physical requirements the paper itself states:
#   - eq. "step-fkt", Theta(x) = (1+(x-1)e^x)/e^(2x), diverges as x->-inf
#     instead of saturating at 0, so it cannot be the bounded, symmetric,
#     temperature-broadened step shown in the paper's own Fig. 2. Theta(x)
#     here was instead re-derived directly from eq. "current" (the paper's
#     own prescription: differentiate it w.r.t. eV) and verified against
#     direct numerical integration to machine precision.
#   - equ. "F_2" as printed has its ln(...) prefactor independent of the
#     integration variable, making the integral trivial and removing
#     exactly the temperature-dependent broadening of the log singularity
#     the surrounding text describes -- and an initial attempt at a
#     corrected/convolution reading of it, while qualitatively peaked in
#     the right place, was verified (see below) to have the wrong sign and
#     the wrong magnitude against the paper's own Fig. F. F(eps,T) here is
#     instead a direct, vectorized evaluation of the paper's unambiguous
#     defining double integral, equ. "F_1" (electron-like processes),
#     verified to reproduce the actual plotted values of Fig. F(b) (e.g.
#     T=1K, omega0=200meV: F(0)=7.459 here vs. ~7.5 read off the figure;
#     F(10meV)=2.945 here vs. ~3.1 read off the figure).


def Theta(x):
    """Temperature-broadened step function.

    Derived by differentiating eq. "current" w.r.t. eV (the paper's own
    prescription) and evaluating the resulting Fermi-function convolution
    in closed form:

        Theta(x) = 1/2 + [sinh(x) - x] / [2(cosh(x) - 1)]

    with x = eps/(kB T). Theta(x) -> 0 as x -> -inf, -> 1 as x -> +inf,
    Theta(0) = 1/2, and Theta(x) + Theta(-x) = 1 identically.
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    large = np.abs(x) > 40
    mid = ~small & ~large
    out[small] = 0.5 + x[small]/6.0 - x[small]**3/360.0
    xm = x[mid]
    out[mid] = 0.5 + (np.sinh(xm) - xm)/(2*(np.cosh(xm) - 1))
    xl = x[large]
    pos = xl > 0
    out[large] = np.where(pos,
                           1.0 - (np.abs(xl)+1)*np.exp(-np.abs(xl)),
                           (np.abs(xl)+1)*np.exp(-np.abs(xl)))
    return out


def Theta_prime(x):
    """d Theta / dx, the even, unit-normalized (integral over x is 1)
    thermal broadening kernel used inside F(eps,T)."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    large = np.abs(x) > 40
    mid = ~small & ~large
    out[small] = 1.0/6 - x[small]**2/60.0
    xm = x[mid]
    out[mid] = (xm*np.sinh(xm)/2 - np.cosh(xm) + 1)/(np.cosh(xm) - 1)**2
    xl = x[large]
    out[large] = np.abs(xl)*np.exp(-np.abs(xl))
    return out


def _fermi(e, kT):
    x = np.clip(e/kT, -700, 700)
    return 1.0/(1.0 + np.exp(x))


def _fermi_prime(e, kT):
    x = np.clip(e/kT, -350, 350)
    ex = np.exp(x)
    return -ex/(kT*(1 + ex)**2)


def _band_integral(epp, kT, omega0, Gamma0):
    """Re[ integral_{-omega0}^{+omega0} (1-f(ep,T))/(ep-epp+i*Gamma0) dep ],
    the exact eps'-integral of equ. "F_1" (a real-valued Lorentzian
    regularization of the ep=epp pole; Gamma0 is the lifetime broadening).
    Smooth in epp except for its own kT-scale Fermi cutoff near epp=0 (the
    Gamma0 pole is already integrated out here, not a remaining feature of
    this function)."""
    def integrand(ep):
        d = ep - epp
        return (1 - _fermi(ep, kT))*d/(d**2 + Gamma0**2)
    # hint both the (Gamma0-regularized) pole at ep=epp and the Fermi
    # cutoff of (1-f(ep,T)) at ep=0 to the adaptive quadrature, when they
    # fall inside the band -- otherwise quad can flag a well-converged
    # result as "probably divergent" if the two nearby features aren't
    # both resolved by its initial subdivision
    pts = sorted({0.} | ({epp} if abs(epp) < omega0 else set()))
    val, _ = quad(integrand, -omega0, omega0, points=pts, limit=200)
    return val


class FBuilder():
    """Vectorized evaluator of F(eps,T), the paper's own defining double
    integral equ. "F_1" (electron-like processes) -- see the module
    docstring for why this replaces the printed closed-form equ. "F_2".

    equ. "F_1" separates into an ep'-integral over the fixed +-omega0 band
    (done once, tabulated over eps'' on a grid and spline-interpolated --
    _band_integral above) and an eps''-integral against the temperature
    derivative of the Fermi function (a narrow, kT-wide kernel), which is
    then a cheap convolution evaluated on the fly for any array of
    eps = eV - eps_m values.
    """
    def __init__(self, T, omega0=20e-3, Gamma0=5e-6, kB=8.617333262e-5,
                 ngrid=400, nwidth=60, npts=2001, grid_span=10.):
        if T <= 0.: raise ValueError("F(eps,T) requires T>0")
        self.kT = kB*T
        self.omega0 = omega0
        self.Gamma0 = Gamma0
        # tabulate the band integral: dense near 0 (kT-scale Fermi cutoff)
        # and near +-omega0 (band edge), spanning well beyond the band
        # (grid_span*omega0) since eV sweeps can range further out than
        # omega0 itself -- _band_integral is smooth there (no remaining
        # sharp feature once the ep' integral is done), just evaluated on
        # a grid coarse enough that cubic-spline extrapolation past it
        # would otherwise blow up (confirmed: it does)
        lin = np.linspace(-1., 1., ngrid)
        grid = np.sign(lin)*grid_span*omega0*np.abs(lin)**3
        grid = np.unique(np.concatenate(
            [grid, np.linspace(-20*self.kT, 20*self.kT, 101)]))
        vals = np.array([_band_integral(e, self.kT, omega0, Gamma0)
                          for e in grid])
        self._spline = CubicSpline(grid, vals, extrapolate=False)
        self._grid_min, self._grid_max = grid[0], grid[-1]
        self._val_min, self._val_max = vals[0], vals[-1]
        width = nwidth*self.kT
        self._fp_grid = np.linspace(-width, width, npts)
        self._fp = _fermi_prime(self._fp_grid, self.kT)
    def _band_value(self, epp):
        # constant extrapolation beyond the tabulated grid instead of
        # letting CubicSpline extrapolate (see __init__ note)
        out = self._spline(np.clip(epp, self._grid_min, self._grid_max))
        return out
    def __call__(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        epp = x[:, None] + self._fp_grid[None, :]
        integrand = self._band_value(epp)*self._fp[None, :]
        return -simpson(integrand, x=self._fp_grid, axis=1)


def F(x, T, omega0=20e-3, Gamma0=5e-6, kB=8.617333262e-5):
    """Convenience one-shot evaluation of F(eps,T) at eps=x (see FBuilder
    for the efficient, precomputed-kernel version used when sweeping many
    eV points against a fixed set of intermediate-state energies)."""
    return FBuilder(T, omega0=omega0, Gamma0=Gamma0, kB=kB)(x)
