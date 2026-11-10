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
# space is exhausted. T=0 only, since the occupied initial state is the
# ground state; the T>0 route stays edkondo.KondoSpectrum.
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


class TaskKondoSpectrum():
    """The ground state and the Ritz states that the impurity spin reaches,
    with the impurity Sx/Sy/Sz in that basis, in the shape of
    edkondo.KondoSpectrum at T=0. p is one-hot on the converged ground
    state, the DMRG route's own T=0 (get_kondo_spectrum's n_gs averages a
    degenerate manifold on top of it)."""
    def __init__(self, chain, site, nkry=100, kB=8.617333262e-5,
                 degeneracy_tol=1e-9):
        from .. import taskdmrg
        taskdmrg._check_backend(chain, "get_kondo_spectrum(mode=\"DMRG\")")
        if isinstance(nkry, bool) or not isinstance(nkry, (int, np.integer)) or nkry < 1:
            raise ValueError("nkry must be a positive integer, got %r" % (nkry,))
        self.chain = chain
        self.site = site
        self.T = 0.
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


    def convergence(self, emax=np.inf):
        """How much of the spectrum below the excitation energy emax the
        Krylov space has converged. Returns a dict: "poles", the number of
        poles with weight there; "unconverged", how many of them have a
        residual norm above CONVERGED_RESIDUAL times the width of the
        spectrum; "weight", the fraction of the spectral weight
        sum_k |<a|S_k|GS>|^2 in that window they carry; and "step", the
        largest of their weights over the total weight, which is the size
        of the largest step of the second-order dI/dV whose position is
        not converged, relative to its largest value."""
        w = sum(np.abs(M[:, 0]) ** 2 for M in (self.Sx, self.Sy, self.Sz))
        total = max(float(np.sum(w)), 1e-300)
        width = max(float(np.max(np.abs(self.e))), 1e-300)
        window = (self.e > 0.) & (self.e <= emax) & (w > 1e-14 * total)
        loose = window & (self.residual > CONVERGED_RESIDUAL * width)
        inside = max(float(np.sum(w[window])), 1e-300)
        return {"poles": int(np.sum(window)),
                "unconverged": int(np.sum(loose)),
                "weight": float(np.sum(w[loose])) / inside,
                "step": float(np.max(w[loose])) / total if loose.any() else 0.}

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
            "%.2g of the spectral weight there. At T=0 every pole is a sharp "
            "step of dI/dV, so at the biases next to these poles the curve "
            "changes with nkry, and from one solve of the ground state to "
            "the next, by up to %.2g of its largest value. A larger nkry "
            "makes this smaller; on a long chain it does not remove it, and "
            "curves are then comparable only at a resolution above the "
            "spacing of the poles."
            % (c["unconverged"], c["poles"], emax, self.nkry, c["weight"],
               c["step"]), RuntimeWarning, stacklevel=stacklevel)
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
