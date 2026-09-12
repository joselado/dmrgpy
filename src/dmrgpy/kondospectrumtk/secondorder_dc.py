import numpy as np

# T=0 second-order dI/dV via the dynamical correlator, avoiding explicit
# excited-state enumeration -- the natural DMRG-side counterpart of
# second_order_dIdV (conductance.py, which sums explicitly over ED
# eigenstates). At T=0 only the ground state is populated, so
#   sum_f |<f|S_alpha|GS>|^2 Theta0(eV-eps_f0)
# is exactly a Theta0-weighted cumulative integral of the T=0 dynamical
# structure factor S_aa(w) = sum_f |<f|S_alpha|GS>|^2 delta(w-eps_f0),
# which chain.get_dynamical_correlator already computes (as a function of
# w, broadened by its own `delta`/`eta` parameter) without ever
# diagonalizing beyond the ground state -- e.g. via submode="KPM" (a
# Chebyshev moment expansion) or "CVM" (correction-vector method), both
# already implemented for itensor_version=3 and itensor_version="python"
# (pyitensor) alike -- confirmed directly, this function works unmodified
# against a pyitensor chain (no compiled extension needed at all).
#
# Convention note: get_dynamical_correlator(name=(A,B)) computes
# G_AB(w) = <GS|A (w-H+E0+i*delta)^-1 B|GS>. To get the POSITIVE spectral
# weight sum_f |<f|B|GS>|^2 (not some other, possibly complex, bilinear
# combination), A must be B's own dagger: name=(B.get_dagger(), B) -- NOT
# name=(B,B) (confirmed directly: the "ED" submode does not
# auto-conjugate its first operator the way submode="EX" does).
#
# Verified (mode="ED", submode="ED", a small S=1/2 Zeeman-split system)
# against the already-validated excited-state-sum second_order_dIdV: sum
# rule (total spectral weight) matches to ~0.06%, and the swept dI/dV
# matches to within ~0.7% (the residual is the expected effect of the
# finite `delta` broadening the otherwise-sharp Theta0 threshold -- it
# shrinks as delta/es resolution are tightened). With mode="DMRG",
# submode="KPM" on a compiled itensor_version=3 backend (3-site chain,
# delta=2e-5, es of 800 points over +-3 meV) the swept dI/dV matches the
# exact sum to 0.2% of its maximum at every bias point, thresholds
# included. It used to be quoted as "a few tens of percent at the
# thresholds": that was the cumulative sum below, not KPM -- see
# _cumulative_theta0_weight (fixed 2026-09-12).


def _cumulative_theta0_weight(chain, op, eVs, mode, submode, delta, es,
                               **kwargs):
    """sum_f |<f|op|GS>|^2 [Theta0(eV-eps_f0) + Theta0(-eV-eps_f0)] for
    every eV in eVs (both tunneling directions, matching
    second_order_dIdV's convention), via a Theta0-weighted cumulative
    integral of the T=0 dynamical structure factor S(w) for
    A=op.get_dagger(), B=op."""
    x, S = chain.get_dynamical_correlator(
            mode=mode, submode=submode, name=(op.get_dagger(), op),
            delta=delta, es=es, **kwargs)
    S = np.asarray(S).real
    # cumulative integral from -inf up to each x point, by the trapezoid
    # rule -- NOT np.cumsum(S)*dw, which counts the whole bin around x_i
    # as lying below x_i and so, exactly at a threshold (where S is a
    # delta-like peak a few bins wide), attributes far more than half of
    # that peak's weight to eV=eps_f0: measured on a 3-site chain with
    # submode="KPM", delta=2e-5, the value at eV=0 came out 1.033
    # against an exact 0.808 while every other point agreed to 1e-4.
    from scipy.integrate import cumulative_trapezoid
    cum = cumulative_trapezoid(S, x, initial=0.)
    return np.interp(eVs, x, cum) + np.interp(-eVs, x, cum)


def second_order_dIdV_dc(chain, site, eVs, T0=1.0, U=0.0, mode="DMRG",
                          submode="KPM", delta=2e-6, es=None, **kwargs):
    """T=0 second-order dI/dV (eq. "conductance"), computed via the
    dynamical correlator instead of the explicit excited-state sum
    conductance.second_order_dIdV uses -- see module docstring. Matches
    that function's return convention (units of 2*pi*e^2*T0^2/hbar,
    U^2 elastic term added directly since it needs no dynamical
    information).

    delta: broadening of the dynamical correlator; controls how sharply
    the Theta0 threshold is resolved (smaller delta -> sharper, but needs
    a correspondingly fine `es` grid to stay well resolved).
    es: frequency grid passed to get_dynamical_correlator. Required (no
    default): it must cover every eigenstate transition energy from the
    ground state that matters for the eVs sweep, at several times finer
    spacing than delta -- these are properties of the chain's own
    spectrum, not of the eVs sweep range, so guessing a default from eVs
    alone is unsafe (confirmed directly: doing so silently misses real
    transitions whenever the eVs sweep happens to be narrower than, or
    offset from, the system's actual energy scale).

    mode/submode: forwarded to chain.get_dynamical_correlator (e.g.
    mode="ED", submode="ED" for an exact small-system check; mode="DMRG",
    submode="KPM"|"CVM" for itensor_version=3)."""
    eVs = np.asarray(eVs, dtype=float)
    if es is None:
        raise ValueError(
            "es must be given explicitly: it needs to cover every "
            "eigenstate transition energy relevant to this system, which "
            "the eVs sweep range alone does not determine")

    Sminus = chain.Sx[site] + (-1j)*chain.Sy[site]
    Splus = chain.Sx[site] + 1j*chain.Sy[site]
    Sz = chain.Sz[site]

    def weight(op):
        return _cumulative_theta0_weight(chain, op, eVs, mode, submode,
                                          delta, es, **kwargs)

    spin_term = 0.5*weight(Sminus) + 0.5*weight(Splus) + weight(Sz)
    # 4*|U|^2, not the bare |U|^2 of eq. "Matrix1sq" -- see the "SA factor"
    # note in conductance.py's module docstring
    U_term = 4*(U**2)*np.ones_like(eVs)
    return 2*np.pi*T0**2*(spin_term + U_term)
