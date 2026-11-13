import numpy as np
from .stepfunctions import Theta, Theta0, FBuilder, F0, cached_fbuilder

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
#
# Direct vs exchange diagram, per term. Every third-order process comes
# in two interaction orders (the paper's "normal" and "reversed", the
# latter tagged R in its Fig. 6/7 -- arXiv v1 numbering, which is what
# this package's "Fig. 3"/"Fig. F" references below mean: its Fig. 3 is
# the arXiv PDF's Fig. 7 and its Fig. F is Fig. 5). Reversing the order
# does two things: it reverses the electron-spin trace, and it turns the
# intermediate electron from electron-like into hole-like, whose F enters
# with the opposite sign (the paper's eqs. 20/21 "with a change of sign").
# For the Kondo term the trace is the antisymmetric Levi-Civita one, so
# the two sign flips cancel and both orders enter with the SAME
# coefficient. For the potential-interference term the trace is the
# symmetric delta_kj, so only the hole-like sign flip survives and the two
# orders enter with OPPOSITE signs, F(eV-eps_im) - F(eV+eps_im). That
# difference is what the paper's Fig. 7c actually shows (121u and 121uR
# are mirror images with opposite sign, their sum is smooth through zero
# bias), and it has two consequences the summed form does not have: the
# elastic intermediate state m=i cancels identically (no sign(eV)*F(eV)
# spike at zero bias), and at B=0 the whole term vanishes -- which is why
# the paper's Fig. 7d (U=0.25) is exactly symmetric at 0 T with the same
# peak position as Fig. 7b. Until 2026-09-12 this term used the summed
# form; see third_order_potential_dIdV's docstring for what that changed.
#
# Where the exchange diagram's log sits. The direct diagram tunnels first:
# its intermediate state is the tunnelled electron in the sample plus the
# impurity in m, on shell at the sample Fermi edge when eV = e_m - e_i, so
# F(eV - eps_im). The exchange diagram scatters a sample electron first,
# into the outgoing state, leaving a hole; that intermediate state is on
# shell when the hole is at the Fermi edge, which fixes the OUTGOING
# electron's energy at e_i - e_m, and the outgoing electron has eV - eps_if
# at the tip's Fermi edge. So the exchange log is F(eV - (e_f - e_m)).
# This DEPARTS FROM THE PAPER: its eq. 25 (and its Fig. 6 caption,
# F(E_m - E_i + eV)) has F(eV - eps_mi) = F(eV + eps_im), which measures
# the log from the incoming energy eV, i.e. as if the impurity absorbed
# nothing. The two agree when f = i -- so the potential term, whose
# I_fi forces f = i, is the same either way -- and whenever every state
# the initial one connects to is degenerate with it (a free S=1/2 at
# B=0, a Kramers doublet). Any field, anisotropy or exchange splitting
# separates them. tests/test_kondo_spectrum_tmatrix.py builds the
# second-order T-matrix of the paper's own Hamiltonian as explicit fermion
# x impurity operators, with every intermediate Fock state's denominator
# read off H0, and this module agrees with it to 1e-15 on anisotropic
# S=1 and exchange-coupled impurities in tilted fields; the printed form
# missed it by 5-20% of the third-order term.
#
# NUMBERS CHANGED on 2026-09-26 for every third-order Kondo term with an
# inelastic transition out of an occupied state (any B != 0, any
# anisotropy, any coupled spin), on mode="ED" and mode="DMRG" alike. For
# the paper's own Fig. 7d at 10 T the step overshoots move from
# 1.226/1.183 to 1.248/1.203 and the +-4 mV tails from 1.143/1.130 to
# 1.149/1.137 (2pi*e^2T0^2/h units); the paper's figure, drawn from its
# eq. 25, reads 1.231/1.177 and 1.146/1.128. Zero-bias values, B=0,
# order=2 and the U term are untouched.


# Thermal occupations below this fraction of the largest one are treated
# as zero when choosing which initial states i to sum over. Every term in
# this module is linear in p_i, so a state dropped here contributes less
# than P_CUT times the largest single-state contribution -- far below the
# accuracy of anything else in the model (Theta/F tabulations are
# ~1e-6..1e-8). It is what makes the third-order terms O(n_occ*dim^2)
# instead of O(dim^3): at T=0 n_occ is the ground-state degeneracy, and at
# any T it is the number of states within ~28 kT of the ground state.
P_CUT = 1e-12


# Final states f whose step Theta(eV - (e_f - e_i)) is closed for every
# occupied i and every bias are left out of the second- and third-order sums
# (_final_states). At T=0 the step of such a state is exactly 0, so nothing
# changes; at T>0 the cut is THETA_CUT kT above the largest open energy, where
# Theta is below Theta(-40) = 1.6e-16, so every term that goes is below about
# 1e-16 of a term that stays. The intermediate state m is virtual and is never
# cut. On the 10-spin chain at T=0 this keeps 99 final states of 1024.
THETA_CUT = 40.


# Byte budgets of the sums. CHUNK_ENTRIES is the size of one block of the
# Theta and F arrays (entries, not bytes); COEFF_BYTES bounds the complex
# Levi-Civita coefficient block of the third-order Kondo term, and
# EXCHANGE_BYTES the two accumulators of its exchange diagram (one per
# tunneling direction), which sets how many bias points share one pass over
# the blocks. Peak memory of the third-order term is a few times the larger
# of the two and does not grow with the number of occupied states.
CHUNK_ENTRIES = 4e6
COEFF_BYTES = 200e6
EXCHANGE_BYTES = 600e6


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
        return cached_fbuilder(ks.T, omega0=omega0, Gamma0=Gamma0, kB=ks.kB)
    return Fb


def _spin_matrix_elements_squared(ks, rows=None, cols=None):
    """|M_if|^2 = 1/2|<f|S-|i>|^2 + 1/2|<f|S+|i>|^2 + |<f|Sz|i>|^2,
    eq. "SA-transitionmatrix" (unpolarized tip and sample). rows and cols
    (index arrays of f and i) restrict it to a block; None is every state."""
    def block(S):
        S = np.asarray(S)
        if rows is not None: S = S[rows]
        if cols is not None: S = S[:, cols]
        return S
    Sx, Sy, Sz = block(ks.Sx), block(ks.Sy), block(ks.Sz)
    Splus = Sx + 1j*Sy
    Sminus = Sx - 1j*Sy
    return 0.5*np.abs(Sminus)**2 + 0.5*np.abs(Splus)**2 + np.abs(Sz)**2


def _final_states(ks, occ, eVs):
    """Indices of the final states f whose step is open for at least one
    occupied initial state i and one bias of either tunneling direction:
    min_i (e_f - e_i) <= max|eV| at T=0 (Theta0(0) = 1/2 keeps the equality),
    and <= max|eV| + THETA_CUT kT at T>0. The difference is the same float
    the sums form, so at T=0 every state left out has a step of exactly 0.
    Every occupied state is kept (e_i - e_i = 0)."""
    eVs = np.asarray(eVs, dtype=float)
    vmax = np.max(np.abs(eVs)) if eVs.size else 0.
    if ks.T != 0.: vmax = vmax + THETA_CUT*ks.kB*ks.T
    lowest = np.min(ks.e[None, :] - ks.e[occ, None], axis=0)
    return np.nonzero(lowest <= vmax)[0]


def second_order_dIdV(ks, eVs, T0=1.0, U=0.0, p_cut=P_CUT):
    """Second-order (Fermi golden rule) dI/dV, eq. "conductance", summed
    over both tunneling directions. U is the dimensionless potential-
    scattering ratio (eq. "Matrix1"); its interference with the spin-
    exchange term (eq. "Matrix1sq"'s third term, the origin of magneto-
    resistive tunneling) vanishes identically for unpolarized tip/sample
    (2*Tr[I s_a] = 0) and is not included here (a polarized-lead
    generalization would need the appendix's eq. "Matrix-appendix").

    The elastic potential-scattering channel enters as 4*|U|^2, not eq.
    "Matrix1sq"'s bare |U|^2 -- see the "SA factor" note in this module's
    docstring, and Fig. 3d minus Fig. 3b (0.25 for the paper's U=0.25).

    The sum runs over the thermally occupied initial states only, the same
    ones as the third-order terms (p_cut, see P_CUT; p_cut=0 keeps every
    state), and over the final states whose step is open (_final_states),
    so it costs O(n_eV n_f n_occ) and not O(n_eV dim^2). The Theta
    arrays are made for blocks of bias points of about CHUNK_ENTRIES
    entries."""
    eVs = np.asarray(eVs, dtype=float)
    occ = _occupied_states(ks, p_cut)
    fin = _final_states(ks, occ, eVs)
    M2 = _spin_matrix_elements_squared(ks, fin, occ) # [f,i], f over fin, i over occ
    eps = ks.e[fin][:, None] - ks.e[occ][None, :] # eps[f,i] = e_f - e_i
    weight = ks.p[occ][None, :]*M2 # weight[f,i] = p_i*|M_if|^2
    chunk = max(1, int(CHUNK_ENTRIES)//weight.size)
    spin_term = np.empty(len(eVs))
    for s0 in range(0, len(eVs), chunk):
        v = eVs[s0:s0+chunk]
        x_ts = v[:, None, None] - eps[None, :, :]
        x_st = -v[:, None, None] - eps[None, :, :]
        spin_term[s0:s0+len(v)] = (
            np.einsum('fi,efi->e', weight, _theta_raw(ks, x_ts))
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


def _occupied_states(ks, p_cut=P_CUT):
    """Indices of the initial states whose Boltzmann weight is not
    negligible (p_i > p_cut*max(p)); p_cut=0 keeps every state."""
    return np.nonzero(ks.p > p_cut*np.max(ks.p))[0]


def _triple_product_coefficients(ks, occ=None, fin=None):
    """coeff[i,f,m] = sum_jkl eps_jkl <i|S_l|f><f|S_k|m><m|S_j|i>, the
    Levi-Civita triple product appearing in eq. "3rd-normal"/"3rd-reversed"
    (both give the identical coefficient; only the F(...) argument
    differs between the direct and exchange diagrams, see the module
    docstring for which).

    occ: indices of the initial states i to build it for (the first axis
    of the result then runs over occ, in that order); None means all of
    them, i.e. the full dim^3 tensor -- only worth paying for when every
    state is thermally occupied. fin: the same for the final states f (the
    second axis); None is every state."""
    Xi = np.stack([ks.Sx, ks.Sy, ks.Sz], axis=-1) # Xi[a,b,alpha]=<a|S_alpha|b>
    if occ is None: occ = np.arange(ks.dim)
    if fin is None: fin = np.arange(ks.dim)
    return np.einsum('jkl,ifl,fmk,mij->ifm', _EPS3, Xi[np.ix_(occ, fin)],
                     Xi[fin, :, :], Xi[:, occ, :], optimize=True)


def third_order_kondo_dIdV(ks, eVs, Jrho_s, T0=1.0, omega0=20e-3, Gamma0=5e-6,
                            Fb=None, p_cut=P_CUT):
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

    The direct diagram enters as F(eV - eps_im) and the exchange diagram
    as F(eV - (e_f - e_m)), NOT the paper's printed F(eV - eps_mi): see
    the module docstring for why, and for what that changed on
    2026-09-26.

    This is an O(n_occ * n_f * dim * len(eVs)) calculation, n_occ being the
    number of thermally occupied initial states (p_cut, see P_CUT: one
    state at T=0 for a non-degenerate ground state) and n_f the number of
    final states whose step is open (_final_states, THETA_CUT), inherent to
    the triple sum over eigenstates -- the m sum always runs over the full
    spectrum, since m is a virtual intermediate state. The exchange
    diagram's F depends on the (f, m) pair rather than on (i, m), so it
    costs n_f*dim F evaluations per bias point whatever n_occ is; they are
    made in chunks of bias points, which bounds the memory at a few
    million entries. It used to be O(dim^3) in both time and memory
    regardless of T (the full coefficient tensor was built for every i,
    3.4 GB and 13 s at dim=512 even at T=0, where all but one row of it is
    multiplied by p_i=0).

    The coefficient tensor is built for blocks of occupied initial states
    (COEFF_BYTES), with both tunneling directions inside each block, so the
    memory does not grow with n_occ; one block, the usual case, costs what a
    single pass costs. Several blocks accumulate the exchange diagram's
    sum over i into one array per direction (EXCHANGE_BYTES) and contract it
    with F afterwards.

    Fb: an existing FBuilder(ks.T, omega0=omega0, Gamma0=Gamma0, kB=ks.kB)
    to reuse instead of building a new one -- building it tabulates
    _band_integral over hundreds of adaptive-quadrature points, so
    callers that also need third_order_potential_dIdV with the same
    T/omega0/Gamma0 (e.g. Spin_Chain.get_kondo_spectrum) should build it
    once and pass it to both. Ignored at ks.T==0 (uses the closed-form F0
    instead)."""
    eVs = np.asarray(eVs, dtype=float)
    occ = _occupied_states(ks, p_cut) # initial states i that carry weight
    fin = _final_states(ks, occ, eVs) # final states f whose step is open
    p = ks.p[occ]
    # eps_if[i,f] = e_f - e_i over occ x fin, eps_im[i,m] = e_m - e_i over
    # occ x every state, e_fm[f,m] = e_f - e_m over fin x every state
    eps_if = ks.e[None, fin] - ks.e[occ, None]
    eps_im = ks.e[None, :] - ks.e[occ, None]
    e_fm = ks.e[fin, None] - ks.e[None, :]
    Fcall = _get_F(ks, omega0, Gamma0, Fb)
    dim, nocc, nf = ks.dim, len(occ), len(fin)
    # blocks of initial states for the coefficient tensor, and groups of
    # bias points sharing one pass over the blocks (see COEFF_BYTES)
    nblock = int(min(nocc, max(1, COEFF_BYTES//(16*nf*dim))))
    group = int(max(1, EXCHANGE_BYTES//(16*nf*dim)))
    sub = max(1, int(CHUNK_ENTRIES)//(dim*nblock)) # biases per Theta/Fd array
    chunk = max(1, int(CHUNK_ENTRIES)//(nf*dim)) # biases per Fx array
    out = np.zeros(len(eVs))
    if nblock >= nocc:
        # one block: the coefficients are built once and both diagrams are
        # contracted per block of bias points, with no exchange accumulator
        # (this is what the code did before the blocking, at the same cost)
        pc = p[:, None, None]*np.imag(_triple_product_coefficients(ks, occ, fin))/2.
        for d, sign in enumerate((1., -1.)):
            v = sign*eVs
            for s0 in range(0, len(v), chunk):
                vc = v[s0:s0+chunk]
                nv = len(vc)
                Th = _theta_raw(ks, vc[:, None, None] - eps_if[None, :, :])
                Fd = Fcall((vc[:, None, None] - eps_im[None, :, :]).ravel()
                           ).reshape(nv, nocc, dim)
                Fx = Fcall((vc[:, None, None] - e_fm[None, :, :]).ravel()
                           ).reshape(nv, nf, dim)
                out[s0:s0+nv] += (np.einsum('ifm,eif,eim->e', pc, Th, Fd, optimize=True)
                                  + np.einsum('ifm,eif,efm->e', pc, Th, Fx, optimize=True))
        return 4*np.pi*T0**2*Jrho_s*out
    for g0 in range(0, len(eVs), group):
        vg = eVs[g0:g0+group]
        ng = len(vg)
        # the s->t direction is the same expression at -v (only eV flips,
        # not eps_if/eps_im/e_fm), so both directions share each block of
        # coefficients; W[d][e,f,m] = sum_i p_i coeff[i,f,m] Th[e,i,f]
        dirs = (vg, -vg)
        W = [np.zeros((ng, nf, dim)) for _ in dirs]
        res = np.zeros(ng)
        for b0 in range(0, nocc, nblock):
            sl = slice(b0, b0+nblock)
            # SA factor: see above
            pc = p[sl, None, None]*np.imag(
                _triple_product_coefficients(ks, occ[sl], fin))/2.
            for d, v in enumerate(dirs):
                for s0 in range(0, ng, sub):
                    vc = v[s0:s0+sub]
                    nv = len(vc)
                    Th = _theta_raw(ks, vc[:, None, None] - eps_if[None, sl, :])
                    # direct diagram: log at eV = e_m - e_i
                    Fd = Fcall((vc[:, None, None] - eps_im[None, sl, :]).ravel()
                               ).reshape(nv, -1, dim)
                    res[s0:s0+nv] += np.einsum('ifm,eif,eim->e', pc, Th, Fd,
                                               optimize=True)
                    W[d][s0:s0+nv] += np.einsum('ifm,eif->efm', pc, Th,
                                                optimize=True)
        # exchange diagram: log at eV = e_f - e_m
        for d, v in enumerate(dirs):
            for s0 in range(0, ng, chunk):
                vc = v[s0:s0+chunk]
                Fx = Fcall((vc[:, None, None] - e_fm[None, :, :]).ravel()
                           ).reshape(len(vc), nf, dim)
                res[s0:s0+len(vc)] += np.einsum('efm,efm->e', W[d][s0:s0+len(vc)], Fx)
        out[g0:g0+ng] = res # eq. "sym_z": same sign in both directions
    return 4*np.pi*T0**2*Jrho_s*out


def third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0, omega0=20e-3,
                                Gamma0=5e-6, Fb=None, p_cut=P_CUT):
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

    The direct and exchange diagrams enter with OPPOSITE signs,
    F(eV-eps_im) - F(eV+eps_im), not the summed combination the Kondo term
    has (with f = i the exchange log F(eV - (e_f - e_m)) is F(eV+eps_im),
    the paper's eq. 25 argument, so this term is unaffected by the
    2026-09-26 change to the Kondo term) -- see the module docstring for
    why (a symmetric electron trace
    leaves the hole-like sign flip of the reversed order uncompensated),
    and the paper's Fig. 3c (arXiv v1: Fig. 7c), where the 121u and 121uR
    curves are mirror images of opposite sign and their sum is smooth
    through zero bias. Two consequences: the m=i term cancels identically
    (no zero-bias spike), and the whole term vanishes at B=0, so the
    zero-field Kondo peak keeps its position and symmetry at U!=0 (the
    paper's Fig. 3d/7d at 0 T is exactly symmetric, peak 1.39 at eV=0).

    NUMBERS CHANGED on 2026-09-12: this term used F(eV-eps_im)+F(eV+eps_im)
    before, which (i) kept a spurious m=i contribution, a sign(eV)*F(eV)
    spike at zero bias that shifted and tilted the zero-field peak
    (measured on the paper's own Fig. 7d parameters: peak at -0.2 mV
    instead of 0, 1.47 instead of 1.39, tails 1.156/1.052 instead of
    1.145/1.135), and (ii) gave the exchange diagram the wrong sign, so
    the field-split step asymmetry at 10 T came out ~0.27 in the paper's
    units against ~0.06 in its Fig. 7d. Every U!=0, order=3 result from
    before this date is not comparable (U=0 and order=2 are untouched).

    The two tunneling directions also enter with OPPOSITE signs, as
    h(eV)-h(-eV) -- eq. "asym_U": "the conductance for processes that
    include potential scattering changes its sign when inverting the
    tunneling direction". The result is therefore odd in eV (in particular
    exactly 0 at eV=0), which is precisely the bias-asymmetric lineshape
    of Fig. 3c (its "sum" curve is antisymmetric about zero bias) and the
    source of the asymmetry in Fig. 3d.

    Fb: see third_order_kondo_dIdV's docstring -- pass the same FBuilder
    to both to avoid rebuilding its expensive tabulation twice. Ignored at
    ks.T==0 (uses the closed-form F0 instead)."""
    eVs = np.asarray(eVs, dtype=float)
    occ = _occupied_states(ks, p_cut) # initial states i that carry weight
    p = ks.p[occ]
    Xi = np.stack([ks.Sx, ks.Sy, ks.Sz], axis=-1) # Xi[a,b,alpha]=<a|S_alpha|b>
    # loop[i,m] = sum_k |<i|Sk|m>|^2, first axis over occ
    loop = np.real(np.einsum('imk,mik->im', Xi[occ, :, :], Xi[:, occ, :]))
    eps_im = ks.e[None, :] - ks.e[occ, None] # eps_im[i,m] = e_m - e_i
    Fcall = _get_F(ks, omega0, Gamma0, Fb)
    dim, nocc = ks.dim, len(occ)

    def one_direction(v):
        Fim = Fcall((v[:, None, None] - eps_im[None, :, :]).ravel()).reshape(len(v), nocc, dim)
        Fmi = Fcall((v[:, None, None] + eps_im[None, :, :]).ravel()).reshape(len(v), nocc, dim)
        # direct MINUS exchange: the m=i (eps_im=0) terms cancel here
        weighted = np.einsum('i,im,eim->e', p, loop, Fim - Fmi, optimize=True)
        return weighted*_theta_raw(ks, v)

    total = one_direction(eVs) - one_direction(-eVs) # eq. "asym_U": sign flips
    return 4*np.pi*T0**2*Jrho_s*U*total
