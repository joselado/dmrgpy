import warnings
import numpy as np
from .stepfunctions import F0, Theta0

# T=0 third-order potential-interference dI/dV (conductance.py's
# third_order_potential_dIdV) via the dynamical correlator, the DMRG-side
# counterpart of secondorder_dc.py's second_order_dIdV_dc -- same idea,
# different kernel. conductance.third_order_potential_dIdV's T=0 limit
# (ks.p one-hot on the ground state) collapses its i,m sum, per tunneling
# direction, to
#   h(eV) = Theta0(eV) * sum_m sum_k |<m|Sk|GS>|^2 [F0(eV-eps_m0)-F0(eV+eps_m0)]
# (k running over Sx,Sy,Sz, matching that function's Xi[a,b,alpha] stack;
# direct MINUS exchange diagram -- see conductance.py's module docstring
# for why this term, unlike the Kondo one, has them with opposite signs,
# and third_order_potential_dIdV's docstring for the 2026-09-12 change
# from the summed form, which this module followed at the same time),
# with the measured term the odd combination h(eV)-h(-eV) over the two
# tunneling directions (eq. "asym_U"; see conductance.py's module
# docstring). Because the bracket vanishes at eps_m0=0, the m=GS
# (elastic) part of S(w) at w=0 drops out of the convolution.
# sum_m |<m|Sk|GS>|^2 delta(w-eps_m0) is exactly the T=0 dynamical
# structure factor S_kk(w) that get_dynamical_correlator already computes
# -- so the m-sum becomes an F0-weighted convolution of S_kk against the
# es frequency grid, in the same spirit as second_order_dIdV_dc's
# Theta0-weighted cumulative integral (a cumulative sum there vs. a
# genuine convolution here, since F0 -- unlike Theta0 -- is not a step).
#
# What es has to cover is NOT the same as for the second-order term. The
# Theta0 kernel there vanishes for every transition above max|eV|, so
# only those below it matter. The kernel here does not: for w above |eV|
# F0(eV-w)-F0(eV+w) is 2*eV*omega0/(w*(omega0+w)) to leading order in
# eV/w, i.e. about 2*eV/w inside the band, so every transition up to the
# top of the site's S_k spectrum contributes at every bias, and the ED
# reference sums all of them. Measured on a 3-site chain (Zeeman-split impurity, lowest
# transition 0.77 meV, lines at 10 to 16 meV carrying 0.40 of the weight
# S(S+1)=0.75): an es over +-3 meV dropped 0.0108 of a 0.1110 peak at
# eV=+-2 meV, all of it from those lines (2026-09-24 hole hunt, finding
# 14). The sum rule sum_k int S_kk(w) dw = S(S+1) is exact, so the
# weight the grid misses is known for free, and
# third_order_potential_dIdV_dc warns when it is more than 1e-2 of it.


def _convolved_F0_weight(chain, op, eVs, omega0, Gamma0, mode, submode,
                          delta, es, **kwargs):
    """(weighted, weight): weighted[e] = int dw S(w) [F0(eV_e-w) -
    F0(eV_e+w)] for every eV in eVs, and weight = int dw S(w), both by the
    trapezoid rule on the correlator's own grid, where S(w) = sum_m
    |<m|op|GS>|^2 delta(w-eps_m0) is the T=0 dynamical structure factor
    for A=op.get_dagger(), B=op (see secondorder_dc.py's module docstring
    for why the explicit get_dagger() is required).

    The quadrature weights are the trapezoid ones of the actual grid, so
    an es that is non-uniform and in any order (refined around a line,
    or a refined block appended after a coarse grid) is weighted
    correctly. It used to be the first spacing times a plain sum, which
    mis-weighted every point by the ratio of the first spacing to the
    local one: 101.7 times the exact value at the peak on a smooth
    sinh-mapped grid dense at the line (2026-09-24 hole hunt, finding 13).
    On a uniform grid the two differ only by the two endpoint half-bins.

    The grid is sorted first, since trapezoid weights from np.diff are
    signed: on a coarse grid followed by a refined block the join counted
    as a negative spacing, and the result was 0.654 off a 0.671 peak
    (2026-09-24b audit, finding 17). A stable sort leaves an increasing
    grid bit for bit as it was."""
    x, S = chain.get_dynamical_correlator(
            mode=mode, submode=submode, name=(op.get_dagger(), op),
            delta=delta, es=es, **kwargs)
    x = np.asarray(x, dtype=float)
    S = np.asarray(S)
    order = np.argsort(x, kind="stable")
    x, S = x[order], S[order]
    S = S.real
    dx = np.diff(x)
    wts = np.zeros_like(x) # trapezoid weights, built once for every eV
    wts[1:] += dx/2
    wts[:-1] += dx/2
    diff = eVs[:, None] - x[None, :]
    ssum = eVs[:, None] + x[None, :]
    kernel = (F0(diff.ravel(), omega0=omega0, Gamma0=Gamma0).reshape(diff.shape)
              - F0(ssum.ravel(), omega0=omega0, Gamma0=Gamma0).reshape(ssum.shape))
    return np.einsum('w,w,ew->e', wts, S, kernel), np.dot(wts, S)


def _spin_casimir(chain, site):
    """S(S+1) of a spin site, or None when the site is not a plain spin
    (the sum rule below then has no fixed value to check against)"""
    from ..spinchain import Spin_Chain
    if not isinstance(chain, Spin_Chain): return None
    try: d = int(chain.sites[site])
    except (TypeError, ValueError, IndexError): return None
    if d < 2: return None
    return (d*d - 1)/4. # S=(d-1)/2


def third_order_potential_dIdV_dc(chain, site, eVs, Jrho_s, U, T0=1.0,
                                   omega0=20e-3, Gamma0=5e-6, mode="DMRG",
                                   submode="KPM", delta=2e-6, es=None,
                                   **kwargs):
    """T=0 third-order potential-interference dI/dV (eq. "U-M"), computed
    via the dynamical correlator instead of the explicit excited-state
    sum conductance.third_order_potential_dIdV uses -- see module
    docstring. Matches that function's return convention (both tunneling
    directions included, so the result is odd in eV) and carries the same
    general-S extrapolation caveat -- see third_order_potential_dIdV's own
    docstring -- on top of the usual delta/es-resolution error already
    present in second_order_dIdV_dc.

    delta: see second_order_dIdV_dc's docstring.
    es: required, and with a STRICTER coverage than second_order_dIdV_dc
    asks for: it must reach past the top of the site's S_k spectrum, i.e.
    every transition from the ground state that carries weight in
    sum_k |<m|S_k|GS>|^2, plus several delta, and include w=0. The
    second-order kernel vanishes above max|eV|, so there only the
    transitions below it matter; this kernel falls off only as about
    2*eV/w inside the band, so every line contributes at every bias (see
    the module docstring for a measured example, where the lines above
    the grid were 10 per cent of the peak). Its spacing must resolve
    delta near the lines, and the grid may be non-uniform and in any
    order. Since sum_k int
    S_kk(w) dw = S(S+1) exactly, a spin site gets a RuntimeWarning when
    the grid misses more than 1e-2 of that weight -- a check of the
    coverage, not of the resolution, and one that counts the elastic
    weight at w=0 too.

    mode/submode: forwarded to chain.get_dynamical_correlator, as in
    second_order_dIdV_dc."""
    eVs = np.asarray(eVs, dtype=float)
    if es is None:
        raise ValueError(
            "es must be given explicitly: it needs to cover every "
            "eigenstate transition energy relevant to this system, which "
            "the eVs sweep range alone does not determine")
    ops = (chain.Sx[site], chain.Sy[site], chain.Sz[site])
    # both tunneling directions, as h(eV)-h(-eV) (eq. "asym_U", see
    # conductance.third_order_potential_dIdV). Sweeping [eVs, -eVs] in one
    # array keeps this to a single dynamical-correlator call per operator
    # -- S(w) does not depend on eV, only the F0 kernel does.
    nev = len(eVs)
    both = np.concatenate([eVs, -eVs])
    weighted, weight = 0., 0.
    for op in ops:
        wk, sk = _convolved_F0_weight(chain, op, both, omega0, Gamma0,
                                      mode, submode, delta, es, **kwargs)
        weighted, weight = weighted + wk, weight + sk
    casimir = _spin_casimir(chain, site)
    if casimir is not None and abs(casimir - weight) > 1e-2*casimir:
        warnings.warn(
            "third_order_potential_dIdV_dc: the es grid [%.4g, %.4g] holds "
            "%.4g of the spectral weight sum_k int S_kk = S(S+1) = %.4g. "
            "This term needs every transition up to the top of the site's "
            "S_k spectrum, not only those below max|eV|, so the missing "
            "weight is missing from the result; widen es past the highest "
            "transition by several delta (the sum rule also counts the "
            "elastic weight at w=0, so es must include w=0)."
            % (np.min(es), np.max(es), weight, casimir), RuntimeWarning,
            stacklevel=2)
    h = weighted*Theta0(both)
    total = h[:nev] - h[nev:]
    return 4*np.pi*T0**2*Jrho_s*U*total
