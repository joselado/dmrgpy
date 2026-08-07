"""
Linear-prediction extrapolation of a real-time correlator C(t), following
White & Affleck (arXiv:0801.4930) and Barthel, White & Schollwoeck,
"Spectral functions in one-dimensional quantum systems at T>0",
PRB 79, 245101 (2009) (arXiv:0901.2342): fit an autoregressive (AR) model
to the tail of the directly-simulated C(t), then run that model's own
recursion forward to synthesize samples far beyond the last directly
simulated time. This is the standard fix, in the DMRG/TDVP dynamical-
correlator literature, for the resolution loss caused by Fourier-
transforming only a finite simulated time window -- it sharpens the
resulting spectral function without any additional real-time evolution
(and thus no additional bond-dimension growth), unlike simply damping/
windowing a fixed-length series harder. See
docs/td_dynamical_correlator_sharpening_plan.md for the full literature
list and how this module composes with `timedependent._damping_window`.

TeNPy's `SpectralSimulation` (GPL-3.0, same license as this project) has
a built-in "linear_prediction" post-processing option for the same
purpose; consulted only for confirming this is established practice, not
copied from -- this module's AR fit/stabilization is its own
implementation.
"""
import numpy as np


def linear_predict_extend(ts, cs, order=20, extend_factor=10,
        fit_start_fraction=0.5, max_pole_radius=1.0, svd_rcond=1e-8):
    """
    Extrapolate a uniformly-sampled complex time series `cs(ts)` using
    linear prediction.

    Fits `cs[n] ~= sum_{k=1}^{order} a_k*cs[n-k]` on the tail of the
    series (the last `1-fit_start_fraction` fraction of samples, skipping
    the early-time transient/short-range physics by default) via an
    SVD-regularized least-squares solve of the linear-prediction normal
    equations (plain Yule-Walker can be ill-conditioned for the near-
    degenerate exponentials typical of these correlators). The fitted
    coefficients are then stabilized (`_stabilize_ar_coefficients`) before
    being run forward autoregressively to generate `extend_factor*len(cs)`
    additional samples at the same spacing `dt=ts[1]-ts[0]`.

    Returns `(ts_full, cs_full)`: the original series with the
    extrapolated tail appended, uniformly spaced, length
    `len(cs)*(1+extend_factor)`.

    Intended to run on the raw, undamped correlator, before any damping/
    windowing is applied (see `timedependent._fourier_transform_correlator`'s
    `predict=` kwarg) -- damping should act on the now-longer effective
    series, not the other way around.
    """
    ts = np.asarray(ts, dtype=float)
    cs = np.asarray(cs, dtype=complex)
    n = len(cs)
    if order < 1 or order >= n // 2:
        raise ValueError(
            "linear_predict_extend: order=%d must be >=1 and < len(cs)//2=%d"
            % (order, n // 2))
    dt = ts[1] - ts[0]

    start = int(fit_start_fraction * n)
    fit = cs[start:]
    m = len(fit)
    if m <= order:
        raise ValueError(
            "linear_predict_extend: not enough samples after "
            "fit_start_fraction=%g to fit order=%d (only %d left)"
            % (fit_start_fraction, order, m))

    # linear-prediction data matrix: D[r,k] = fit[order+r-1-k], predicting
    # b[r] = fit[order+r] as sum_k a[k]*D[r,k], i.e. a[k] is the weight of
    # lag (k+1)
    D = np.empty((m - order, order), dtype=complex)
    for k in range(order):
        D[:, k] = fit[order - 1 - k:m - 1 - k]
    b = fit[order:]
    a, _res, _rank, _sv = np.linalg.lstsq(D, b, rcond=svd_rcond)
    a = _stabilize_ar_coefficients(a, max_pole_radius=max_pole_radius)

    n_extend = int(extend_factor * n)
    ext = np.concatenate([cs[-order:], np.zeros(n_extend, dtype=complex)])
    for i in range(n_extend):
        lags = ext[i:i + order][::-1]  # lags[k] = c[n-1-k], n=i+order
        ext[i + order] = np.dot(a, lags)

    cs_full = np.concatenate([cs, ext[order:]])
    ts_full = ts[0] + dt * np.arange(len(cs_full))
    return ts_full, cs_full


def _stabilize_ar_coefficients(a, max_pole_radius=1.0):
    """
    Reflect any root of the AR characteristic polynomial that lies
    outside `max_pole_radius` back onto that radius, preserving its
    phase, then rebuilds the AR coefficients from the corrected roots.

    A least-squares AR fit has no built-in guarantee that every pole of
    `c[n]=sum_k a_k*c[n-k]` lies within the unit circle; a pole fit
    slightly outside it (from statistical noise in a finite, undamped
    signal that itself has no intrinsic decay) makes the forward
    recursion diverge instead of extrapolating a bounded oscillation.
    Since the physical correlator only has weight for poles on or inside
    the unit circle (unitary/marginally-stable evolution, aside from
    numerical noise), clipping is a correction of a fitting artifact, not
    an approximation of new physics -- the standard "root reflection"
    fix used in linear-prediction spectral estimation (e.g. LPSVD).
    """
    p = len(a)
    char_poly = np.concatenate(([1. + 0j], -np.asarray(a, dtype=complex)))
    roots = np.roots(char_poly)
    mag = np.abs(roots)
    too_big = mag > max_pole_radius
    if np.any(too_big):
        roots[too_big] = roots[too_big] / mag[too_big] * max_pole_radius
    new_poly = np.poly(roots)  # leading coefficient 1, length p+1
    return -new_poly[1:]
