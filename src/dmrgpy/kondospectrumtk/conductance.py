import numpy as np
from .stepfunctions import Theta, Theta0, FBuilder, F0

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
# Both tunneling directions (t->s and s->t) are summed at EVERY order,
# since a measured dI/dV is for the net current I = I^{t->s} - I^{s->t}
# (eq. "current"). For unpolarized tip and sample the matrix elements are
# identical in both directions (the polarization-dependent prefactors
# (1+-eta)/2 both reduce to 1/2 at eta=0), and the s->t contribution
# works out to the same expression with x=-eV-eps_if in place of
# x=eV-eps_if -- i.e. only the sign of eV (not of eps_if) flips between
# the two directions, so each term below is assembled as
#
#     g(eV) +- g(-eV),
#
# with g the t->s expression. The relative sign is fixed by the paper's
# own worked (121)/(121u) example, eqs. "sym_z"/"asym_U": the purely
# Kondo-like (spin-exchange) processes "[contribute] positive for both
# tunneling directions" (-> PLUS, an even-in-eV term), while "the
# conductance for processes that include potential scattering changes its
# sign when inverting the tunneling direction" (-> MINUS, an odd-in-eV
# term, which is exactly the bias asymmetry of Fig. 3c/d).
#
# Normalization ("SA factor"): the paper's spin-averaged transition
# matrix element, eq. "SA-transitionmatrix", is 2x the plain sum over the
# four electron-spin channels of |<..|s.S|..>|^2 built from its own
# electron matrix elements eq. "arvgspinmatrix" (s = sigma/2, so
# <up|s_z|up>=1/2 etc.). Writing that as W = 2*Tr_electron[...] fixes the
# normalization of every order at once, and reproduces the paper's own
# absolute figures:
#   - 2 vertices, spin: 2*Tr[s_a s_b] = delta_ab -> exactly eq.
#     "SA-transitionmatrix" (checked below in _spin_matrix_elements_squared)
#   - 2 vertices, potential: 2*Tr[I I]*|U|^2 = 4|U|^2 -- NOT the bare
#     |U|^2 of eq. "Matrix1sq" (which is written for the unnormalized
#     M^(1)). Fig. 3d minus Fig. 3b, dashed second-order curves, is
#     exactly 0.25 = 4*(0.25)^2 for the paper's U=0.25.
#   - 3 vertices, spin: 2*Tr[s_l s_k s_j] = (i/2) eps_lkj, so the
#     Levi-Civita coefficient is Im[X]/2, not Im[X]/4 (the naive reading
#     of eq. "3rd-normal"'s printed Re[.]/(4i) prefactor). Checked
#     against Fig. 3b/3d directly: for the paper's S=1/2, B=0, T=1K,
#     Jrho_s=-0.05, omega0=20meV, that gives a zero-bias peak of 1.137
#     (Fig. 3b reads 1.13) and 1.387 with U=0.25 (Fig. 3d reads 1.39).
#   - 3 vertices, one potential: 2*Tr[I s_k s_j] = delta_kj, i.e. the
#     plain sum_k |<m|S_k|i>|^2 already used below.
# See examples/kondo/kondo_spectrum_VS_paper for these checks as runnable
# assertions.


def _theta_raw(ks, x):
    """Theta evaluated at a raw energy argument x (not yet divided by kT);
    ks.T==0 dispatches to the closed-form Heaviside limit Theta0(x)."""
    if ks.T == 0.: return Theta0(x)
    return Theta(x/(ks.kB*ks.T))


def _get_F(ks, omega0, Gamma0, Fb):
    """Return a callable F(x_array)->array: the closed-form F0 at ks.T==0,
    otherwise a (possibly caller-supplied, shared) FBuilder."""
    if ks.T == 0.:
        return lambda x: F0(x, omega0=omega0, Gamma0=Gamma0)
    if Fb is None:
        return FBuilder(ks.T, omega0=omega0, Gamma0=Gamma0, kB=ks.kB)
    return Fb


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
    (2*Tr[I s_a] = 0) and is not included here (a polarized-lead
    generalization would need the appendix's eq. "Matrix-appendix").

    The elastic potential-scattering channel enters as 4*|U|^2, not eq.
    "Matrix1sq"'s bare |U|^2 -- see the "SA factor" note in this module's
    docstring, and Fig. 3d minus Fig. 3b (0.25 for the paper's U=0.25)."""
    eVs = np.asarray(eVs, dtype=float)
    M2 = _spin_matrix_elements_squared(ks) # [f,i]
    eps = ks.e[:, None] - ks.e[None, :] # eps[f,i] = e_f - e_i
    weight = ks.p[None, :]*M2 # weight[f,i] = p_i*|M_if|^2
    x_ts = eVs[:, None, None] - eps[None, :, :]
    x_st = -eVs[:, None, None] - eps[None, :, :]
    spin_term = (np.einsum('fi,efi->e', weight, _theta_raw(ks, x_ts))
                 + np.einsum('fi,efi->e', weight, _theta_raw(ks, x_st)))
    # elastic potential scattering: eps_ii=0, Theta(eV)+Theta(-eV)=1
    # (Theta0(0)+Theta0(-0)=1 too, by the same convention)
    U_term = 4*(U**2)*np.ones_like(eVs)
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
    intermediate state), and over both tunneling directions.

    Both directions matter and are summed here as g(eV)+g(-eV) (see this
    module's docstring): eq. "sym_z" states the Kondo-like processes
    contribute with the *same* sign in both directions, so the result is
    even in eV. That is what makes the zero-field feature the symmetric
    zero-bias resonance of Fig. 3a/b -- each individual process curve
    plotted there (112/112R, 121, 121R, 122, 122R) is itself symmetric
    about eV=0, which the t->s term alone (a one-sided Theta step) cannot
    be. The overall normalization carries the "SA factor" 2 of this
    module's docstring, i.e. the Levi-Civita coefficient is Im[X]/2.

    This is an O(dim^3 * len(eVs)) calculation, inherent to the triple sum
    over eigenstates -- expected to be the bottleneck for larger Hilbert
    spaces.

    Fb: an existing FBuilder(ks.T, omega0=omega0, Gamma0=Gamma0, kB=ks.kB)
    to reuse instead of building a new one -- building it tabulates
    _band_integral over hundreds of adaptive-quadrature points, so
    callers that also need third_order_potential_dIdV with the same
    T/omega0/Gamma0 (e.g. Spin_Chain.get_kondo_spectrum) should build it
    once and pass it to both. Ignored at ks.T==0 (uses the closed-form F0
    instead)."""
    eVs = np.asarray(eVs, dtype=float)
    coeff = np.imag(_triple_product_coefficients(ks))/2. # SA factor: see above
    # eps_if[i,f] = e_f - e_i and eps_im[i,m] = e_m - e_i are the same
    # array (e[None,:]-e[:,None]) under two names for readability at the
    # call sites below
    eps_if = eps_im = ks.e[None, :] - ks.e[:, None]
    Fcall = _get_F(ks, omega0, Gamma0, Fb)
    dim = ks.dim

    def one_direction(v):
        """d(I^{t->s})/dV at bias v; the s->t direction is this same
        expression at -v (only eV flips, not eps_if/eps_im)."""
        Fim = Fcall((v[:, None, None] - eps_im[None, :, :]).ravel()).reshape(len(v), dim, dim)
        Fmi = Fcall((v[:, None, None] + eps_im[None, :, :]).ravel()).reshape(len(v), dim, dim)
        Th = _theta_raw(ks, v[:, None, None] - eps_if[None, :, :])
        Fsum = Fim + Fmi # direct (eps_im) + exchange (eps_mi=-eps_im) diagrams
        return np.einsum('i,ifm,eif,eim->e', ks.p, coeff, Th, Fsum, optimize=True)

    total = one_direction(eVs) + one_direction(-eVs) # eq. "sym_z": same sign
    return 4*np.pi*T0**2*Jrho_s*total


def third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0, omega0=20e-3,
                                Gamma0=5e-6, Fb=None):
    """Third-order potential-scattering interference term, eq. "U-M"
    (the origin of the bias asymmetry in Fig. 3c/d), summed over both
    tunneling directions.

    The paper only spells out this term's electron-spin-averaged
    (unpolarized-lead) form via a worked S=1/2 example, not a general-S
    formula. This implementation extrapolates to general S by replacing
    eq. "3rd-normal"'s antisymmetric Levi-Civita triple product (which
    comes from the electron-spin trace of three s=sigma/2 operators,
    2*Tr[s_l s_k s_j] ~ i*eps_lkj) with the symmetric contraction
    2*Tr[I s_k s_j] = delta_kj that the electron identity at the outer,
    U-carrying vertex gives instead -- i.e. sum_k <i|Sk|m><m|Sk|i> =
    sum_k |<i|Sk|m>|^2, manifestly real. eq. "U-M"'s I_fi=delta_fi forces
    the impurity to return to its initial state (elastic), collapsing the
    i,f,m sum to a two-state i,m loop.

    Unlike third_order_kondo_dIdV, the two tunneling directions enter with
    OPPOSITE signs here, as h(eV)-h(-eV) -- eq. "asym_U": "the conductance
    for processes that include potential scattering changes its sign when
    inverting the tunneling direction". The result is therefore odd in eV
    (in particular exactly 0 at eV=0), which is precisely the
    bias-asymmetric lineshape of Fig. 3c (its "sum" curve is antisymmetric
    about zero bias) and the source of the asymmetry in Fig. 3d.

    Fb: see third_order_kondo_dIdV's docstring -- pass the same FBuilder
    to both to avoid rebuilding its expensive tabulation twice. Ignored at
    ks.T==0 (uses the closed-form F0 instead)."""
    eVs = np.asarray(eVs, dtype=float)
    Xi = np.stack([ks.Sx, ks.Sy, ks.Sz], axis=-1) # Xi[a,b,alpha]=<a|S_alpha|b>
    loop = np.real(np.einsum('imk,mik->im', Xi, Xi)) # sum_k |<i|Sk|m>|^2
    eps_im = ks.e[None, :] - ks.e[:, None] # eps_im[i,m] = e_m - e_i
    Fcall = _get_F(ks, omega0, Gamma0, Fb)
    dim = ks.dim

    def one_direction(v):
        Fim = Fcall((v[:, None, None] - eps_im[None, :, :]).ravel()).reshape(len(v), dim, dim)
        Fmi = Fcall((v[:, None, None] + eps_im[None, :, :]).ravel()).reshape(len(v), dim, dim)
        weighted = np.einsum('i,im,eim->e', ks.p, loop, Fim + Fmi, optimize=True)
        return weighted*_theta_raw(ks, v)

    total = one_direction(eVs) - one_direction(-eVs) # eq. "asym_U": sign flips
    return 4*np.pi*T0**2*Jrho_s*U*total
