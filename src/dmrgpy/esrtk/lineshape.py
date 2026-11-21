"""The steady state of the driven rate equation and its readout.

With the drive Omega V cos(w t) and every line resolved (the secular
approximation, a rotating-wave approximation per transition), the
populations obey

    (R + W(w)) p = 0,   sum p = 1,

R the relaxation of rates.py and W the driven rates, the golden rule with
a Lorentzian, for every pair n (upper), m (lower):

    W_nm(w) = (Omega^2/2) |V_nm|^2 gamma_nm/((w - w_nm)^2 + gamma_nm^2).

Writing W = -sum_p w_p a_p a_p^T with a_p = e_n - e_m, and y_p = p_n - p_m,
the steady state is y = y0 + K diag(w) y, K = A^T G A, G = rates.resolvent,
and the change of the populations is G A diag(w) y. A pair is solved
exactly (it saturates) when its saturation parameter s_p = -w_p K_pp, which
is Omega_R^2 T1 T2/(1 + (w - w_nm)^2 T2^2) for a two-level system, is above
SAT_TOL at that frequency; the others enter at first order in their w_p,
with the populations the saturated pairs leave. For one spin 1/2 this is
the Bloch steady state. The coherence of each pair follows its populations,

    rho_nm = -(Omega/2) V_nm y_p/(w - w_nm + i gamma_nm)   (times e^{-iwt}),

and the readout of the current, in units of G_j a V_DC, is

    dI = sum_k (p_k - p0_k) D_kk + (V_RF/V_DC) Re[e^{i phase} sum rho_nm D_mn].

Frequencies f are in GHz and angular frequencies w = 2 pi f in rad/ns.
"""
import warnings

import numpy as np

from . import rates

# A pair is driven only when |V_nm|^2 is above V_TOL of the largest one.
V_TOL = 1e-20
# A pair outside the window [f_min, f_max] is kept when it is closer to it
# than MARGIN of its full widths; the tail of a Lorentzian at MARGIN widths
# is 1/(4 MARGIN^2) of its peak.
MARGIN = 1000.
# The states kept are those up to WINDOW_MARGIN full widths of a free spin
# (1/(pi T2) GHz) above the top of the window, plus kT ln(1/p_cut) at T > 0
# (the populated states): a line further out than that is dropped, and its
# DC tail inside the window is below 1/(4 WINDOW_MARGIN^2) = 3e-4 of its own
# peak. Among the kept states, pairs are kept up to MARGIN widths.
WINDOW_MARGIN = 30.
# The homodyne part of a line is dispersive, with a tail 1/detuning that
# a sum over the kept states cuts off wherever emax happens to end: it
# moved by 1e-2 to 4e-2 of the DC line when the window was zoomed onto a
# line (6 and 12 Heisenberg spins). Each line's coherence is therefore
# multiplied by homodyne_taper, 1 up to HOMODYNE_INNER free-spin widths
# from its resonance and 0 beyond WINDOW_MARGIN, a smoothstep in between.
# The far dispersive background is left out: it would need every state up
# to the exchange scale and the counter-rotating terms, and experiments
# subtract a smooth background anyway.
HOMODYNE_INNER = 15.
# A line is read, for the position warning and for the energy errors the
# DMRG Krylov route keeps, when its position, give or take the error bound
# of its two energies, comes within READ_MARGIN half widths 1/(2 pi T2) of
# a free spin of the window; a line further out adds only a tail.
READ_MARGIN = 5.


def homodyne_taper(det, T2):
    """The factor of a coherence at a detuning |det| (rad/ns) from its
    resonance: 1 up to HOMODYNE_INNER full widths 2/T2 of a free spin, 0
    beyond WINDOW_MARGIN, and a smoothstep in between"""
    x1 = HOMODYNE_INNER*2./T2
    x2 = WINDOW_MARGIN*2./T2
    t = np.clip((np.abs(det) - x1)/(x2 - x1), 0., 1.)
    return 1. - t*t*(3. - 2.*t)
# Pairs with a saturation parameter above SAT_TOL at a frequency are solved
# exactly there; the error of the others is of order SAT_TOL of their own
# small contribution.
SAT_TOL = 1e-6
# Two states closer than DEG_TOL of the width of the kept spectrum (and
# 1e-15 in absolute) are degenerate.
DEG_TOL = 1e-10
# A matrix element below COUPLING_TOL of the largest one does not couple
COUPLING_TOL = 1e-10
# The position bound of a line is compared with its half width for the lines
# with a weight above WARN_WEIGHT of the strongest
WARN_WEIGHT = 1e-3


class Lines():
    """Everything about the kept spectrum that does not depend on the
    frequency: rates, populations, pairs, resolvent."""
    def __init__(self, spec, f, h, kT, T1, T2, rabi, p_cut, f_ref):
        f = np.asarray(f, dtype=float)
        if f.ndim != 1 or len(f) == 0:
            raise ValueError("f must be a one-dimensional array of frequencies")
        self.spec = spec
        e = spec.e
        self.kT = kT
        self.T2 = float(T2)
        self.fmin, self.fmax = float(np.min(f)), float(np.max(f))
        # f_ref: the frequency at which a free spin 1/2 relaxes in T1
        self.f_ref = float(f_ref)
        self.g = rates.korringa_scale(T1, h*self.f_ref, kT)
        self.gphi = rates.dephasing_rate(T1, T2)
        self.Gamma = rates.rate_matrix(e, spec.K, self.g, kT)
        self.R = rates.generator(self.Gamma)
        self.out = self.Gamma.sum(axis=0) # total rate out of each state
        self.p0 = rates.stationary(e, kT)
        self.w_n = 2.*np.pi*e/h # angular frequency of each level, rad/ns
        self.Omega = 4.*np.pi*rabi # Omega_R = Omega |V_nm|
        n = spec.dim
        up, lo = np.triu_indices(n, k=1)
        # sort each pair as (upper, lower)
        swap = e[up] < e[lo]
        u = np.where(swap, lo, up)
        l = np.where(swap, up, lo)
        width = max(np.max(e) - np.min(e), 0.)
        self.deg = DEG_TOL*width + 1e-15
        # on the MPS route two Ritz energies are only known to sigma_n +
        # sigma_m, so a pair closer than that counts as degenerate
        degp = np.full(len(u), self.deg)
        if spec.sigma is not None:
            degp = degp + spec.sigma[u] + spec.sigma[l]
        self.degp = degp
        v2 = np.abs(spec.V[u, l])**2
        vmax = np.max(v2) if len(v2) else 0.
        gam = 0.5*(self.out[u] + self.out[l]) + self.gphi
        fp = (e[u] - e[l])/h
        marg = MARGIN*gam/np.pi # MARGIN full widths, in GHz
        inwin = (fp >= np.min(f) - marg) & (fp <= np.max(f) + marg)
        cand = (e[u] - e[l] > degp) & (v2 > V_TOL*vmax) & inwin & (vmax > 0.)
        # the states that can carry population: the populated ones and every
        # state a driven pair reaches from them
        pop = self.p0 > p_cut*np.max(self.p0)
        while True:
            act = cand & (pop[u] | pop[l])
            new = pop.copy()
            new[u[act]] = True
            new[l[act]] = True
            if np.array_equal(new, pop): break
            pop = new
        self._check_degenerate(u, l, pop)
        self.u, self.l = u[act], l[act]
        self.v2 = v2[act]
        self.gam = gam[act]
        self.w_p = self.w_n[self.u] - self.w_n[self.l]
        self.f_p = fp[act]
        self.y0 = self.p0[self.u] - self.p0[self.l]
        self.Vul = spec.V[self.u, self.l]
        self.Dlu = spec.D[self.l, self.u]
        self.Ddiag = np.real(np.diag(spec.D))
        self.npairs = len(self.u)
        self.labels = rates.blocks(self.Gamma)
        if self.npairs:
            cross = self.labels[self.u] != self.labels[self.l]
            if np.any(cross):
                i = np.flatnonzero(cross)[0]
                raise ValueError(
                    "get_esr_spectrum: the drive connects the states %d and %d, "
                    "which the relaxation does not connect (different blocks "
                    "of the rate matrix, for instance two charge sectors), so "
                    "the steady state is not defined; give the relaxation a "
                    "path between them (site_weights)" % (self.l[i], self.u[i]))
            G = rates.resolvent(self.R, self.p0, self.labels, e)
            self.G = G
            gD = G.T @ self.Ddiag
            self.d = gD[self.u] - gD[self.l] # D^T G a_p
            self.Kdiag = (G[self.u, self.u] - G[self.u, self.l]
                          - G[self.l, self.u] + G[self.l, self.l])

    def _check_degenerate(self, u, l, pop):
        """Raise on two degenerate states that V, D or the relaxation couple,
        at least one of them populated or reached by a driven pair (pop):
        their populations and D_kk depend on the basis the eigensolver
        picked inside the level."""
        e = self.spec.e
        close = np.abs(e[u] - e[l]) <= self.degp
        if not np.any(close): return
        uu, ll = u[close], l[close]
        def coupled(M):
            M = np.abs(np.asarray(M))
            m = np.max(M) if M.size else 0.
            return M[uu, ll] > COUPLING_TOL*max(m, 1e-300)
        c = coupled(self.spec.V) | coupled(self.spec.D) | coupled(self.spec.K)
        c &= pop[uu] | pop[ll]
        if np.any(c):
            i = np.flatnonzero(c)[0]
            raise ValueError(
                "get_esr_spectrum: the states %d and %d are degenerate (energy "
                "%.6g, split by %.3g, below the tolerance %.3g, which on "
                "mode=\"DMRG\" includes the energy uncertainty of the two "
                "states) and coupled by the drive, the detection or the "
                "relaxation, so their populations depend on the basis inside "
                "the level. This is a zero field or a level crossing; ESR "
                "needs the levels split (method=\"floquet\" on a small system "
                "handles close levels, and on DMRG converge the states "
                "further)" % (ll[i], uu[i], e[uu[i]], abs(e[uu[i]] - e[ll[i]]),
                              self.degp[close][i]))

    def K(self, rows, cols):
        """K[p, q] = a_p^T G a_q for the pairs rows x cols"""
        G = self.G
        ur, lr = self.u[rows], self.l[rows]
        uc, lc = self.u[cols], self.l[cols]
        return (G[np.ix_(ur, uc)] - G[np.ix_(ur, lc)]
                - G[np.ix_(lr, uc)] + G[np.ix_(lr, lc)])

    def solve(self, w):
        """(y, wp): the population differences of the pairs and their driven
        rates at the angular frequency w"""
        if self.npairs == 0:
            return np.zeros(0), np.zeros(0)
        dw = w - self.w_p
        wp = 0.5*self.Omega**2*self.v2*self.gam/(dw**2 + self.gam**2)
        s = -wp*self.Kdiag
        near = np.flatnonzero(s > SAT_TOL)
        far = np.flatnonzero(s <= SAT_TOL)
        y = self.y0.copy()
        if len(near):
            # y_N = y0_N + K_NN w_N y_N + K_NF w_F y_F, with the far pairs at
            # first order in their w_F: y_F = y0_F + K_FN w_N y_N
            KNN = self.K(near, near)
            rhs = self.y0[near].copy()
            A = np.eye(len(near)) - KNN*wp[near][None, :]
            if len(far):
                KNF = self.K(near, far)*wp[far][None, :]
                rhs += KNF @ self.y0[far]
                A -= (KNF @ self.K(far, near))*wp[near][None, :]
            yN = np.linalg.solve(A, rhs)
            y[near] = yN
            if len(far):
                y[far] = self.y0[far] + self.K(far, near) @ (wp[near]*yN)
        return y, wp

    def signal(self, f, vrf_over_vdc=0., phase=0.):
        """(dc, homodyne) at the frequencies f (GHz)"""
        f = np.asarray(f, dtype=float)
        dc = np.zeros(len(f))
        hom = np.zeros(len(f))
        if self.npairs == 0:
            return dc, hom
        ph = np.exp(1j*phase)
        for i, fi in enumerate(f):
            w = 2.*np.pi*fi
            y, wp = self.solve(w)
            dc[i] = np.sum(wp*y*self.d)
            if vrf_over_vdc != 0.:
                rho = -0.5*self.Omega*self.Vul*y/(w - self.w_p + 1j*self.gam)
                rho = rho*homodyne_taper(w - self.w_p, self.T2)
                hom[i] = vrf_over_vdc*np.real(ph*np.sum(rho*self.Dlu))
        return dc, hom

    def poles(self, h):
        """The driven transitions, one entry per pair"""
        D = np.real(np.diag(self.spec.D))
        out = dict(frequency=self.f_p.copy(), upper=self.u.copy(),
                   lower=self.l.copy(), weight=self.v2.copy(),
                   form_factor=D[self.u] - D[self.l],
                   population_difference=-self.y0,
                   width=self.gam/np.pi) # FWHM in GHz at weak drive
        if self.spec.sigma is not None:
            s = self.spec.sigma
            out["position_bound"] = (s[self.u] + s[self.l])/h
        return out

    def warn_positions(self, h):
        """On the MPS route, warn when the bound on the position of a strong
        line (sigma_n + sigma_m, which has an eigenvalue within it) is
        larger than its half width. Only lines read in the window count
        (READ_MARGIN): one further out adds only a tail, and on
        submode="KRYLOV" the states such a line ends on can be relaxation
        partners whose energies are not converged."""
        if self.spec.sigma is None or self.npairs == 0: return
        s = self.spec.sigma
        bound = (s[self.u] + s[self.l])/h
        strength = np.abs(self.v2*self.y0)
        rm = READ_MARGIN/(2.*np.pi*self.T2)
        near = ((self.f_p + bound >= self.fmin - rm)
                & (self.f_p - bound <= self.fmax + rm))
        strength = strength*near
        if not np.any(strength > 0.): return
        strong = strength > WARN_WEIGHT*np.max(strength)
        hw = 0.5*self.gam/np.pi
        bad = strong & (bound > hw)
        if np.any(bad):
            i = np.flatnonzero(bad)[np.argmax((bound/hw)[bad])]
            warnings.warn(
                "get_esr_spectrum(%s): the position of the line at %.6g GHz is "
                "bounded only to %.3g GHz by the error of the energies of its "
                "two states, more than its half width %.3g GHz; converge the "
                "states further (maxm, nsweeps or nkry) before reading line "
                "positions at that level" % (self.spec.route, self.f_p[i],
                                             bound[i], hw[i]),
                RuntimeWarning, stacklevel=4)
