import numpy as np

# The T=0 Kondo spectrum from the tangent-space Krylov (TaSK) spectrum of
# arXiv:2510.07279, as a drop-in for edkondo.KondoSpectrum.
#
# conductance.py reads a spectrum object through six attributes only (e, p,
# Sx, Sy, Sz in the eigenbasis, T, kB, dim), so the exact second-order,
# third-order Kondo and potential-interference sums run unchanged once the
# full eigenbasis is replaced by a smaller one: the ground state plus the
# Ritz states of P(H-E0)P on the block Krylov space started from
# P S_a|GS>, a = x, y, z, with P the projector on the one-site tangent space
# of the ground-state MPS (pyitensor/task.py).
#
# Why this is the same approximation as the TaSK correlator and no worse:
# the three spin operators act on the one impurity site, so S_j|GS> and
# S_k S_j|GS> are both one-site changes of the ground state and lie in
# span{|GS>} + tangent space exactly. Every column <a|S|GS> is exact; the
# middle vertex <f|S_k|m> is taken between two tangent-space states, as
# <f|P S_k P|m>, which is exact for those states. What is approximate is
# only which states the sums run over: Ritz states of PHP instead of
# eigenstates of H, measured pole by pole by the two-site variance, and
# exact when the bond dimension saturates the Hilbert space and the Krylov
# space is exhausted. The occupied initial state is the ground state, so
# T>0 here is the temperature of the electrodes only (the thermal steps and
# the thermal Kondo log of conductance.py), which is the thermal state of
# the whole problem while kT is far below the lowest excitation; above
# that the route is edkondo.KondoSpectrum, and the constructor raises.
#
# Which Ritz states are eigenstates of PHP is measured, by the residual norm
# ||(PHP - e)|a>|| the block Krylov returns for each of them. It matters at
# T=0, where every pole is a sharp step of dI/dV: a pole with a residual is
# a quadrature node that stands for a group of eigenstates, its energy moves
# with nkry, and it moved between two solves of the same chain whose states
# had an overlap of 1 to 1e-15 (24-site S=1/2 chain, maxm=40, nkry=100: 21
# of the poles with weight converged, and dI/dV different by 3e-3 of its
# peak at the biases a pole had crossed, by 1.3e-2 against nkry=200). The
# converged poles agreed to 3e-13. convergence() reports it and
# get_kondo_spectrum warns. A larger nkry shrinks the open share (0.61 at
# 100, 0.09 at 400 on that chain) and does not close it on a long chain.

# A Ritz state counts as an eigenstate of PHP when its residual norm is
# below this fraction of the width of the spectrum; its energy is then
# right to the square of that over the gap to its neighbours.
CONVERGED_RESIDUAL = 1e-8
# get_kondo_spectrum warns when a pole that is not converged makes a step
# larger than this fraction of the largest second-order dI/dV.
UNCONVERGED_STEP = 1e-3
# The largest Boltzmann factor exp(-e1/kT) of the lowest pole above the
# ground state at which the sample still counts as being in its ground
# state; T>0 raises above it.
THERMAL_OCCUPATION = 1e-3
# The thermal step Theta has a derivative of full width 5.4 kT. Poles that
# are not converged do not matter once they are closer to each other than
# that: the sum over them is then a quadrature of a smooth function. On a
# 24-site chain at nkry=100 the steps are 9e-5 eV apart; at 0.3 K (width
# 1.4e-4 eV) the curve is within 1e-5 of the nkry=200 one and two solves
# agree to 9e-7, at 0.1 K (4.7e-5 eV) within 8e-4 and 2e-4, at T=0 within
# 1.3e-2 and 3e-3.
THERMAL_WIDTH = 5.4


class TaskKondoSpectrum():
    """The ground state and the Ritz states that the impurity spin reaches,
    with the impurity Sx/Sy/Sz in that basis, in the shape of
    edkondo.KondoSpectrum at T=0. p is one-hot on the converged ground
    state, the DMRG route's own T=0 (get_kondo_spectrum's n_gs averages a
    degenerate manifold on top of it)."""
    def __init__(self, chain, site, nkry=100, kB=8.617333262e-5,
                 degeneracy_tol=1e-9, T=0., n_gs=1):
        from .. import taskdmrg
        taskdmrg._check_backend(chain, "get_kondo_spectrum(mode=\"DMRG\")")
        if isinstance(nkry, bool) or not isinstance(nkry, (int, np.integer)) or nkry < 1:
            raise ValueError("nkry must be a positive integer, got %r" % (nkry,))
        if not T >= 0.:
            raise ValueError("T must be >= 0, got %r" % (T,))
        self.chain = chain
        self.site = site
        self.T = float(T)
        self.kB = kB
        taskdmrg.ground_state(chain) # no session on julia_live
        ops = []
        for name in ("Sx", "Sy", "Sz"):
            op = getattr(chain, name)[site]
            ops.append((op.to_terms(), op.get_dagger().to_terms()))
        self.nkry = int(nkry)
        e, mats, dim, res = taskdmrg.session_spectrum(chain, ops, int(nkry),
                                                      residual=True)
        e = np.asarray(e, dtype=float)
        self.residual = np.asarray(res, dtype=float)
        # A state degenerate with the ground state comes out of the Ritz
        # solve at +-1e-16 rather than at 0, and the zero-bias point then
        # reads the Theta0 step on one side or the other at random: on a
        # Kramers doublet the n_gs=2 average missed ED by 37% at eV=0, in 5
        # runs of 6, and agreed to 1e-15 everywhere else. Snapped to 0
        # within edkondo.KondoSpectrum's own degeneracy_tol, relative to
        # the bandwidth, so the step is read at Theta0(0) = 1/2 as on ED.
        width = max(np.max(np.abs(e)), 1e-300)
        e[np.abs(e) <= degeneracy_tol*width] = 0.
        # A Ritz energy of P(H-E0)P is the energy of an actual state, so one
        # below the snapping window means the state handed over is not the
        # ground state: on a near-degenerate doublet DMRG can stop at a
        # mixture of the two members, and the sums below then measured an
        # order-one deviation from ED (2026-10 audit, finding 3). The
        # exception is a state the impurity spin does not connect to the
        # ground state: on a state promoted from a conserved sector, P
        # commutes with the charge, rounding seeds a lower sector and the
        # Krylov space finds it at the bottom of the spectrum, with no
        # matrix element to anything the sums reach (Ne=21 next to an
        # Ne=20 sector 89 meV lower, on 7 orbitals of 20).
        if np.min(e) < 0.:
            import warnings
            below = np.flatnonzero(e < 0.)
            reached = connected_to_ground_state(mats)
            coupled = below[reached[below]]
            if len(coupled):
                warnings.warn(
                    "get_kondo_spectrum(submode=\"TASK\"): the tangent space "
                    "of the ground state reaches a state %.3g below it, so the "
                    "state DMRG converged to is not the ground state (an "
                    "accidental near-degeneracy, or an unconverged solve). "
                    "Converge further, or pass n_gs= for a degenerate manifold."
                    % (np.min(e[coupled]),), RuntimeWarning, stacklevel=4)
            else:
                warnings.warn(
                    "get_kondo_spectrum(submode=\"TASK\"): %d Ritz states "
                    "lie below the ground state (the lowest %.3g below it), "
                    "but the impurity spin does not connect them to it, so "
                    "they do not enter the spectrum: on a state promoted "
                    "from a conserved sector, rounding seeds a lower charge "
                    "sector. Row 0 is the ground state, not argmin(e)."
                    % (len(below), -np.min(e)), RuntimeWarning, stacklevel=4)
        self.e = e
        self.dim = len(self.e)
        self.tangent_dimension = dim
        self.p = np.zeros(self.dim)
        self.p[0] = 1.
        self.Sx, self.Sy, self.Sz = mats
        if self.T > 0.:
            self._check_cold(mats, n_gs)

    def _check_cold(self, mats, n_gs):
        """Raise unless the state the sums run from is the thermal state at
        self.T: no state the spin reaches may be occupied, beyond the
        degenerate manifold that n_gs declares.

        The members of a manifold are the poles within THERMAL_OCCUPATION
        kT of the ground state, whose Boltzmann factors are 1 to that
        accuracy, so the equal weights of n_gs are the thermal ones; at a
        truncated maxm the partner of a member is not at zero exactly
        (1.5e-9 eV on 9 sites at maxm=10). The lowest pole above them
        enters with its residual subtracted: a Ritz energy is an upper
        bound, and a short Krylov space puts it too high (0.47 meV at
        nkry=3 for 0.29 meV on a 24-site chain), which would let a warm
        sample through."""
        e, kT = self.e, self.kB * self.T
        near = THERMAL_OCCUPATION * kT
        rows = np.arange(1, self.dim)
        reached = connected_to_ground_state(mats)[1:]
        partners = rows[np.abs(e[1:]) <= near]
        if len(partners) > n_gs - 1:
            raise ValueError(
                "get_kondo_spectrum(mode=\"DMRG\", T=%g): %d more state(s) "
                "lie within %.3g of the ground state, so at this "
                "temperature the sample is in a manifold of %d states and "
                "not in the %d that n_gs=%d averages. Pass n_gs=%d, or use "
                "mode=\"ED\"." % (self.T, len(partners), near,
                                  len(partners) + 1, n_gs, n_gs,
                                  len(partners) + 1))
        below = rows[(e[1:] < -near) & reached]
        if len(below):
            raise ValueError(
                "get_kondo_spectrum(mode=\"DMRG\", T=%g): the spin of site "
                "%d reaches a state %.3g below the state DMRG converged to, "
                "which is then not the ground state, and no thermal "
                "spectrum follows from it. Converge the ground state "
                "further, or use mode=\"ED\"."
                % (self.T, self.site, -np.min(e[below])))
        above = rows[e[1:] > near]
        if not len(above):
            return
        k = above[np.argmin(e[above])]
        lowest, floor = e[k], max(e[k] - self.residual[k], 0.)
        if np.exp(-floor / kT) <= THERMAL_OCCUPATION:
            return
        bound = "T below %.3g" % (
            lowest / (self.kB * np.log(1. / THERMAL_OCCUPATION)))
        width = max(float(np.max(np.abs(e))), 1e-300)
        if self.residual[k] > CONVERGED_RESIDUAL * width:
            raise ValueError(
                "get_kondo_spectrum(mode=\"DMRG\", T=%g): the lowest pole, "
                "%.3g above the ground state, is not converged at nkry=%d "
                "(residual %.2g), so the lowest excitation is only known to "
                "be above %.3g, a Boltzmann factor of up to %.2g, above "
                "%.0e. Raise nkry, or use mode=\"ED\"."
                % (self.T, lowest, self.nkry, self.residual[k], floor,
                   np.exp(-floor / kT), THERMAL_OCCUPATION))
        raise ValueError(
            "get_kondo_spectrum(mode=\"DMRG\", T=%g): the lowest "
            "excitation the spin of site %d reaches is %.3g above "
            "the ground state, a Boltzmann factor of %.2g at this "
            "temperature, above %.0e. The DMRG route keeps the "
            "sample in its ground state and only the electrodes at "
            "T, so it is not the thermal spectrum here; use "
            "mode=\"ED\", or %s."
            % (self.T, self.site, lowest, np.exp(-lowest / kT),
               THERMAL_OCCUPATION, bound))

    def convergence(self, emax=np.inf):
        """How much of the spectrum below the excitation energy emax the
        Krylov space has converged. Returns a dict: "poles", the number of
        poles with weight there; "unconverged", how many of them have a
        residual norm above CONVERGED_RESIDUAL times the width of the
        spectrum (and, at T>0, lie further from their neighbours than the
        width of the thermal step, which does not resolve closer ones);
        "weight", the fraction of the spectral weight
        sum_k |<a|S_k|GS>|^2 in that window they carry; and "step", the
        largest of their weights over the total weight, which is the size
        of the largest step of the second-order dI/dV whose position is
        not converged, relative to its largest value. At T>0 a step is
        smooth, and a pole that moves by d changes the curve by about
        d/(THERMAL_WIDTH kT) of its weight, so there each weight in "step"
        is scaled by min(1, d/(THERMAL_WIDTH kT)), with d = min(r, r^2/g)
        the distance from a Ritz energy with residual r to an eigenvalue
        of PHP, g being the distance to the nearest other Ritz energy.
        Measured against nkry=300 on one state at nkry=100: on a gapped
        20-site chain the size went from 4.1e-3 to 3.0e-5 at 0.1 K (true
        change 6.5e-8) and to 1.1e-4 at 0.05 K (7.0e-7), so the warning no
        longer fires; on a uniform 24-site chain it went from 3.7e-2 to
        1.5e-2 at 0.05 K (8.7e-4) and still fires. In both it stays above
        the true change."""
        w = sum(np.abs(M[:, 0]) ** 2 for M in (self.Sx, self.Sy, self.Sz))
        total = max(float(np.sum(w)), 1e-300)
        width = max(float(np.max(np.abs(self.e))), 1e-300)
        window = (self.e > 0.) & (self.e <= emax) & (w > 1e-14 * total)
        loose = window & (self.residual > CONVERGED_RESIDUAL * width)
        if self.T > 0. and loose.any():
            # a thermal step does not resolve poles closer than its width.
            # The distance is the one to the third nearest pole with
            # weight: the three spin operators give up to three poles at
            # nearly the same energy, which are one step, not three
            ew = np.sort(self.e[w > 1e-14 * total])
            k = min(3, len(ew) - 1)
            apart = np.array([np.sort(np.abs(ew - x))[k] if k > 0 else np.inf
                              for x in self.e])
            loose &= apart > THERMAL_WIDTH * self.kB * self.T
        inside = max(float(np.sum(w[window])), 1e-300)
        step = w[loose]
        if self.T > 0. and loose.any():
            step = step * self._thermal_shift_factor(np.flatnonzero(loose))
        return {"poles": int(np.sum(window)),
                "unconverged": int(np.sum(loose)),
                "weight": float(np.sum(w[loose])) / inside,
                "step": float(np.max(step)) / total if loose.any() else 0.}

    def _thermal_shift_factor(self, rows):
        """For the Ritz states in rows, min(1, d/(THERMAL_WIDTH kT)): the
        part of its weight by which a thermal step moves when the pole
        moves by d. A Ritz energy with residual r is within r of an
        eigenvalue of PHP, and within r^2/g when the distance g to the
        nearest other Ritz energy is larger than r, so d = min(r, r^2/g).
        A Ritz energy degenerate with another one gets d = r."""
        es = np.sort(self.e)
        r = self.residual[rows]
        g = np.array([np.sort(np.abs(es - self.e[k]))[1] if len(es) > 1
                      else np.inf for k in rows])
        quad = np.full(len(rows), np.inf)
        np.divide(r ** 2, g, out=quad, where=g > 0.)
        d = np.minimum(r, quad)
        return np.minimum(1., d / (THERMAL_WIDTH * self.kB * self.T))

    def warn_unconverged(self, emax, stacklevel=5):
        """Warn when a pole below emax whose position is not converged
        makes a step of dI/dV above UNCONVERGED_STEP of the largest value."""
        c = self.convergence(emax)
        if c["step"] <= UNCONVERGED_STEP:
            return c
        import warnings
        warnings.warn(
            "get_kondo_spectrum(submode=\"TASK\"): %d of the %d poles below "
            "max|eV| = %.3g are not converged at nkry=%d, and they carry "
            "%.2g of the spectral weight there. %s, so at the biases next "
            "to these poles the curve changes with nkry, and from one solve "
            "of the ground state to the next, by up to %.2g of its largest "
            "value. A larger nkry makes this smaller; on a long chain at "
            "T=0 it does not remove it, and a temperature of the electrodes "
            "above the spacing of the poles does."
            % (c["unconverged"], c["poles"], emax, self.nkry, c["weight"],
               "At T=0 every pole is a sharp step of dI/dV" if self.T == 0.
               else "At T=%g the thermal step is narrower than their "
               "spacing and than how far they can still move" % self.T,
               c["step"]),
            RuntimeWarning, stacklevel=stacklevel)
        return c


def connected_to_ground_state(mats, tol=1e-8):
    """Which rows of the basis the operators in `mats` connect to row 0,
    the ground state, through any chain of matrix elements above tol times
    the largest one. A state outside that set has no amplitude in any of
    the conductance sums, whatever its energy."""
    A = sum(np.abs(M) for M in mats)
    A = (A + A.T) > tol*max(np.max(A), 1e-300)
    reached = np.zeros(len(A), dtype=bool)
    reached[0] = True
    front = [0]
    while front:
        new = np.flatnonzero(A[front].any(axis=0) & ~reached)
        reached[new] = True
        front = list(new)
    return reached
