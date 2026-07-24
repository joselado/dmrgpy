import numpy as np
from .stepfunctions import Theta, FBuilder

# Second- and third-order dI/dV for the STM/Kondo perturbation theory of
# Ternes, New J. Phys. 17 063016 (2015), arXiv:1505.04430, assembled from
# a KondoSpectrum (full ED spectrum + eigenbasis spin matrix elements,
# see edkondo.py).
#
# Throughout, matrices Sx/Sy/Sz stored on a KondoSpectrum are indexed
# [f,i] = <f|S|i>, i.e. directly in the paper's M_if convention (row is
# the *second*, bra, subscript). eps[f,i] = e_f - e_i matches eps_if.
#
# Results are returned in units of 2*pi*e^2*T0^2/hbar (i.e. the physical
# prefactor is folded into T0, which is otherwise a free, experiment-
# dependent tunneling-strength scale -- see eq. "current"/"conductance").
#
# Both tunneling directions (t->s and s->t) are summed, since a measured
# dI/dV is for the net current I = I^{t->s} - I^{s->t} (eq. "current").
# For unpolarized tip and sample the spin-exchange matrix elements are
# identical in both directions (the polarization-dependent prefactors
# (1+-eta)/2 both reduce to 1/2 at eta=0), and the s->t contribution to
# d(-I^{s->t})/d(eV) works out to the same Theta(x) evaluated at
# x=-eV-eps_if (see the derivation note in kondospectrum.py's module
# docstring) -- so only the sign of eV (not of eps_if) flips between the
# two directions.


def _spin_matrix_elements_squared(ks):
    """|M_if|^2 = 1/2|<f|S-|i>|^2 + 1/2|<f|S+|i>|^2 + |<f|Sz|i>|^2,
    eq. "SA-transitionmatrix" (unpolarized tip and sample)."""
    Splus = ks.Sx + 1j*ks.Sy
    Sminus = ks.Sx - 1j*ks.Sy
    return 0.5*np.abs(Sminus)**2 + 0.5*np.abs(Splus)**2 + np.abs(ks.Sz)**2


def second_order_dIdV(ks, eVs, T0=1.0, U=0.0):
    """Second-order (Fermi golden rule) dI/dV, eq. "conductance", summed
    over both tunneling directions. U is the dimensionless potential-
    scattering ratio (eq. "Matrix1"); its interference with the spin-
    exchange term (eq. "Matrix1sq"'s third term, the origin of magneto-
    resistive tunneling) vanishes identically for unpolarized tip/sample
    and is not included here (a polarized-lead generalization would need
    the appendix's eq. "Matrix-appendix")."""
    eVs = np.asarray(eVs, dtype=float)
    kT = ks.kB*ks.T
    M2 = _spin_matrix_elements_squared(ks) # [f,i]
    eps = ks.e[:, None] - ks.e[None, :] # eps[f,i] = e_f - e_i
    weight = ks.p[None, :]*M2 # weight[f,i] = p_i*|M_if|^2
    x_ts = (eVs[:, None, None] - eps[None, :, :])/kT
    x_st = (-eVs[:, None, None] - eps[None, :, :])/kT
    spin_term = (np.einsum('fi,efi->e', weight, Theta(x_ts))
                 + np.einsum('fi,efi->e', weight, Theta(x_st)))
    # elastic potential scattering: eps_ii=0, Theta(eV/kT)+Theta(-eV/kT)=1
    U_term = (U**2)*np.ones_like(eVs)
    return 2*np.pi*T0**2*(spin_term + U_term)


_EPS3 = np.zeros((3, 3, 3))
for _i, _j, _k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    _EPS3[_i, _j, _k] = 1.
for _i, _j, _k in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
    _EPS3[_i, _j, _k] = -1.


def _triple_product_coefficients(ks):
    """coeff[i,f,m] = sum_jkl eps_jkl <i|S_l|f><f|S_k|m><m|S_j|i>, the
    Levi-Civita triple product appearing in eq. "3rd-normal"/"3rd-reversed"
    (both give the identical coefficient; only the F(...) argument
    differs between the direct and exchange diagrams)."""
    Xi = np.stack([ks.Sx, ks.Sy, ks.Sz], axis=-1) # Xi[a,b,alpha]=<a|S_alpha|b>
    return np.einsum('jkl,ifl,fmk,mij->ifm', _EPS3, Xi, Xi, Xi, optimize=True)


def third_order_kondo_dIdV(ks, eVs, Jrho_s, T0=1.0, omega0=20e-3, Gamma0=5e-6,
                            Fb=None):
    """Third-order Kondo term, eqs. "3rd-normal" (direct diagram) +
    "3rd-reversed" (exchange diagram -- despite the name this is NOT the
    s->t tunneling direction: the paper's own "t -R-> s" notation on eq.
    "3rd-reversed" marks it as the reversed *interaction order* diagram,
    still for t->s tunneling; both equations are introduced as replacing
    M_1 with M_1+M_2 inside eq. "conductance", which is explicitly the
    t->s formula), summed over all i,f,m (m unrestricted -- a virtual
    intermediate state).

    Unlike second_order_dIdV, this returns d(I^{t->s})/dV only: the paper
    never gives a general s->t formula for the third-order terms (its own
    worked S=1/2 example, eqs. "sym_z"/"asym_U", shows the two directions
    are related by more than a plain eV -> -eV mirror -- an overall sign
    flip "due to the rearrangement of the interaction order together with
    the switching from hole-like to electron-like scattering"), so no
    such generalization is attempted here. For the field-split zero-bias
    Kondo peak this reproduces (Fig. 3a/b), t->s alone already carries the
    symmetric-in-eV lineshape shown there.

    This is an O(dim^3 * len(eVs)) calculation, inherent to the triple sum
    over eigenstates -- expected to be the bottleneck for larger Hilbert
    spaces.

    Fb: an existing FBuilder(ks.T, omega0=omega0, Gamma0=Gamma0, kB=ks.kB)
    to reuse instead of building a new one -- building it tabulates
    _band_integral over hundreds of adaptive-quadrature points, so
    callers that also need third_order_potential_dIdV with the same
    T/omega0/Gamma0 (e.g. Spin_Chain.get_kondo_spectrum) should build it
    once and pass it to both."""
    eVs = np.asarray(eVs, dtype=float)
    kT = ks.kB*ks.T
    coeff = np.imag(_triple_product_coefficients(ks))/4. # Re[X/(4i)] = Im[X]/4
    eps_if = ks.e[None, :] - ks.e[:, None] # eps_if[i,f] = e_f - e_i
    eps_im = ks.e[None, :] - ks.e[:, None] # eps_im[i,m] = e_m - e_i
    if Fb is None: Fb = FBuilder(ks.T, omega0=omega0, Gamma0=Gamma0, kB=ks.kB)
    dim = ks.dim
    Fim = Fb((eVs[:, None, None] - eps_im[None, :, :]).ravel()).reshape(len(eVs), dim, dim)
    Fmi = Fb((eVs[:, None, None] + eps_im[None, :, :]).ravel()).reshape(len(eVs), dim, dim)
    Th = Theta((eVs[:, None, None] - eps_if[None, :, :])/kT)
    Fsum = Fim + Fmi # direct (eps_im) + exchange (eps_mi=-eps_im) diagrams
    total = np.einsum('i,ifm,eif,eim->e', ks.p, coeff, Th, Fsum, optimize=True)
    return 4*np.pi*T0**2*Jrho_s*total


def third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0, omega0=20e-3,
                                Gamma0=5e-6, Fb=None):
    """Third-order potential-scattering interference term, eq. "U-M"
    (the origin of the bias asymmetry in Fig. 3c/d).

    LOWER CONFIDENCE than second_order_dIdV/third_order_kondo_dIdV: the
    paper gives eq. "U-M" itself unambiguously, but only spells out its
    electron-spin-averaged (unpolarized-lead) closed form -- the one
    actually needed here -- via a worked S=1/2 example (eqs. "sym_z"/
    "asym_U"), not a general-S formula. This implementation extrapolates
    to general S by replacing eq. "3rd-normal"'s antisymmetric
    Levi-Civita triple product (which comes from an electron-spin trace
    of three Pauli matrices, Tr[sigma_l sigma_k sigma_j] ~ i*eps_jkl) with
    the symmetric dot product Tr[sigma_l * I] ~ delta_lk gives when the
    outer, U-carrying vertex is the electron identity instead of a third
    Pauli matrix -- i.e. sum_k <i|Sk|m><m|Sk|i> = sum_k |<i|Sk|m>|^2,
    manifestly real. eq. "U-M"'s I_fi=delta_fi forces the impurity to
    return to its initial state (elastic), collapsing the i,f,m sum to a
    two-state i,m loop. The overall numeric prefactor is not pinned down
    by the paper and is set here to match third_order_kondo_dIdV's
    normalization by analogy; treat this term as provisional until
    checked against Fig. 3c/d directly.

    Like third_order_kondo_dIdV, this returns d(I^{t->s})/dV only -- see
    that function's docstring for why the s->t direction is not attempted
    here (the paper's own worked example shows it is not a plain eV ->
    -eV mirror, which is exactly the bias asymmetry this term is meant to
    produce in the first place).

    Fb: see third_order_kondo_dIdV's docstring -- pass the same FBuilder
    to both to avoid rebuilding its expensive tabulation twice."""
    eVs = np.asarray(eVs, dtype=float)
    Xi = np.stack([ks.Sx, ks.Sy, ks.Sz], axis=-1) # Xi[a,b,alpha]=<a|S_alpha|b>
    loop = np.real(np.einsum('imk,mik->im', Xi, Xi)) # sum_k |<i|Sk|m>|^2
    eps_im = ks.e[None, :] - ks.e[:, None] # eps_im[i,m] = e_m - e_i
    if Fb is None: Fb = FBuilder(ks.T, omega0=omega0, Gamma0=Gamma0, kB=ks.kB)
    dim = ks.dim
    Fim = Fb((eVs[:, None, None] - eps_im[None, :, :]).ravel()).reshape(len(eVs), dim, dim)
    Fmi = Fb((eVs[:, None, None] + eps_im[None, :, :]).ravel()).reshape(len(eVs), dim, dim)
    Fsum = Fim + Fmi
    weighted = np.einsum('i,im,eim->e', ks.p, loop, Fsum, optimize=True)
    total = weighted*Theta(eVs/(ks.kB*ks.T))
    return 4*np.pi*T0**2*Jrho_s*U*total
