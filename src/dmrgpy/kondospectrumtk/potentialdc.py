import numpy as np
from .stepfunctions import F0, Theta0

# T=0 third-order potential-interference dI/dV (conductance.py's
# third_order_potential_dIdV) via the dynamical correlator, the DMRG-side
# counterpart of secondorder_dc.py's second_order_dIdV_dc -- same idea,
# different kernel. conductance.third_order_potential_dIdV's T=0 limit
# (ks.p one-hot on the ground state) collapses its i,m sum to
#   Theta0(eV) * sum_m sum_k |<m|Sk|GS>|^2 [F0(eV-eps_m0)+F0(eV+eps_m0)]
# (k running over Sx,Sy,Sz, matching that function's Xi[a,b,alpha] stack)
# and sum_m |<m|Sk|GS>|^2 delta(w-eps_m0) is exactly the T=0 dynamical
# structure factor S_kk(w) that get_dynamical_correlator already computes
# -- so the m-sum becomes an F0-weighted convolution of S_kk against the
# es frequency grid, in the same spirit as second_order_dIdV_dc's
# Theta0-weighted cumulative integral (a cumulative sum there vs. a
# genuine convolution here, since F0 -- unlike Theta0 -- is not a step).


def _convolved_F0_weight(chain, op, eVs, omega0, Gamma0, mode, submode,
                          delta, es, **kwargs):
    """dw * sum_i S(w_i) * [F0(eV-w_i) + F0(eV+w_i)] for every eV in eVs,
    where S(w) = sum_m |<m|op|GS>|^2 delta(w-eps_m0) is the T=0 dynamical
    structure factor for A=op.get_dagger(), B=op (see secondorder_dc.py's
    module docstring for why the explicit get_dagger() is required)."""
    x, S = chain.get_dynamical_correlator(
            mode=mode, submode=submode, name=(op.get_dagger(), op),
            delta=delta, es=es, **kwargs)
    x = np.asarray(x, dtype=float)
    S = np.asarray(S).real
    dw = x[1] - x[0]
    diff = eVs[:, None] - x[None, :]
    ssum = eVs[:, None] + x[None, :]
    kernel = (F0(diff.ravel(), omega0=omega0, Gamma0=Gamma0).reshape(diff.shape)
              + F0(ssum.ravel(), omega0=omega0, Gamma0=Gamma0).reshape(ssum.shape))
    return dw*np.einsum('w,ew->e', S, kernel)


def third_order_potential_dIdV_dc(chain, site, eVs, Jrho_s, U, T0=1.0,
                                   omega0=20e-3, Gamma0=5e-6, mode="DMRG",
                                   submode="KPM", delta=2e-6, es=None,
                                   **kwargs):
    """T=0 third-order potential-interference dI/dV (eq. "U-M"), computed
    via the dynamical correlator instead of the explicit excited-state
    sum conductance.third_order_potential_dIdV uses -- see module
    docstring. Matches that function's return convention and carries the
    same LOWER CONFIDENCE caveat (general-S extrapolation from the
    paper's worked S=1/2 example -- see third_order_potential_dIdV's own
    docstring) on top of the usual delta/es-resolution error already
    present in second_order_dIdV_dc.

    delta/es: see second_order_dIdV_dc's docstring -- es is required for
    the same reason (it must cover every eigenstate transition energy
    from the ground state, a property of the chain's own spectrum, not
    of the eVs sweep range).

    mode/submode: forwarded to chain.get_dynamical_correlator, as in
    second_order_dIdV_dc."""
    eVs = np.asarray(eVs, dtype=float)
    if es is None:
        raise ValueError(
            "es must be given explicitly: it needs to cover every "
            "eigenstate transition energy relevant to this system, which "
            "the eVs sweep range alone does not determine")
    ops = (chain.Sx[site], chain.Sy[site], chain.Sz[site])
    weighted = sum(_convolved_F0_weight(chain, op, eVs, omega0, Gamma0,
                                         mode, submode, delta, es, **kwargs)
                   for op in ops)
    total = weighted*Theta0(eVs)
    return 4*np.pi*T0**2*Jrho_s*U*total
