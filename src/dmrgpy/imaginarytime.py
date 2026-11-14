"""Imaginary-time evolution of a state by implicit steps on a logarithmic
time grid.

Reference: J. P. Zima, E. M. Stoudenmire, S. R. White, O. Parcollet and
J. Kaye, "Fast Tensor Network Imaginary Time Evolution by Implicit
Stepping on Logarithmic Grids", arXiv:2606.02930.

The state psi(tau) = e^{-tau H} psi(0) is a sum of decaying exponentials
e^{-E_n tau}, so the energy scale that matters falls as tau grows, and a
grid with steps that grow exponentially resolves it: P panels with
endpoints 0, tau/2^(P-1), tau/2^(P-2), ..., tau, n uniform steps on each,
O(log(tau E_max)) steps in total against O(tau E_max) for uniform steps.
Long steps need an A-stable stepper, one whose discrete propagator stays
below one in modulus on every positive eigenvalue whatever the step size,
and only implicit ones are: here the implicit trapezoid rule (order 2),

    (1 + dtau/2 (H-s)) psi_k = (1 - dtau/2 (H-s)) psi_{k-1},

or the s-stage Gauss-Legendre rule (order 2s), whose propagator on a
linear equation is the diagonal Pade approximant Q(-z)/Q(z) of e^{-z},
applied here as the product of its s factors (1 + z/q_j)/(1 - z/q_j) over
the roots q_j of Q, each one linear solve. The paper solves its spectral
(Legendre collocation) variant as one coupled system per panel, whose
details it leaves to a later publication; the product form gives the same
end-of-step propagator. The shift s is the ground-state energy by default,
so that H - s >= 0 and every factor is A-stable. It has to sit at E0 and
not merely below it: neither rule is L-stable, the propagator of a mode
tends to (-1)^s once the step is long against its energy, so a shift a
distance m below E0 stops the projection once the steps exceed about
1/m, the ground state being sent to -1 with everything else. The energy
of the current state, the paper's choice for the Anderson model, is the
other option; it approaches E0 from above, and while it is still far
above it a long step can meet the pole of the trapezoid factor on the
states below it.

Each factor is a two-site sweep with local GMRES (pyitensor/implicit.py,
mpscpp3's Chain::implicit_step), the state is normalized after every
step, and the cost of one step does not grow with its size, so at fixed
bond dimension the cost of reaching tau grows as log(tau). mode="ED" is
the exact e^{-tau H} on the same grid, the reference.
"""

from math import factorial

import numpy as np


def gauss_legendre_roots(stages):
    """The roots q_j of Q(z) = sum_k c_k z^k, the denominator of the
    (stages, stages) Pade approximant e^{-z} ~ Q(-z)/Q(z), with
    c_k = (2s-k)! s! / ((2s)! k! (s-k)!). They lie in the left half-plane;
    stages=1 gives q = -2, the trapezoid rule."""
    s = int(stages)
    if s < 1:
        raise ValueError("gauss_legendre_roots: stages must be >= 1")
    c = [factorial(2*s-k)*factorial(s)/(factorial(2*s)*factorial(k)*factorial(s-k))
         for k in range(s+1)]
    return np.roots(c[::-1])


def step_factors(dtau, order=2):
    """The betas of one step of size dtau: x -> (1 + b(H-s))^-1 (1 - b(H-s)) x
    for each b, applied in turn. order=2 is the trapezoid rule (b = dtau/2),
    order=2s the s-stage Gauss-Legendre rule (b = -dtau/q_j)."""
    if order % 2 or order < 2:
        raise ValueError("imaginary-time order must be an even integer >= 2 "
                         "(got %r)" % (order,))
    if order == 2:
        return [complex(dtau/2.)]
    return [complex(-dtau/q) for q in gauss_legendre_roots(order//2)]


def log_grid(tau, tau0, nsteps):
    """The step sizes of the logarithmic grid from 0 to tau: panels
    [0,t1], [t1,2t1], [2t1,4t1], ..., [tau/2,tau] with t1 = tau/2^(P-1) the
    largest such value at most tau0, and nsteps uniform steps on each.
    tau <= tau0 gives one panel."""
    tau, tau0 = float(tau), float(tau0)
    if tau <= 0.: raise ValueError("imaginary time tau must be > 0")
    if tau0 <= 0.: raise ValueError("tau0 must be > 0")
    if int(nsteps) < 1: raise ValueError("nsteps must be >= 1")
    P = 1 if tau <= tau0 else int(np.ceil(np.log2(tau/tau0) - 1e-12)) + 1
    t1 = tau/2.**(P-1)
    ends = [0.] + [t1*2.**j for j in range(P)]
    steps = []
    for a, b in zip(ends[:-1], ends[1:]):
        steps += [(b-a)/int(nsteps)]*int(nsteps)
    return np.array(steps)


def _unit(H):
    """The largest coefficient of a MultiOperator, its energy unit."""
    return max([abs(t[0]) for t in H.op]+[0.])


def _plan(self, tau, tau0, nsteps, order):
    H = self.hamiltonian
    if H is None:
        raise ValueError("evolve_imaginary_time: set a Hamiltonian first")
    if tau0 is None:
        u = _unit(H)
        tau0 = 0.01/u if u > 0. else 0.01
    steps = log_grid(tau, tau0, nsteps)
    step_factors(steps[0], order) # validates order
    return steps


def _check_dmrg(self):
    if self.itensor_version not in (3, "python"):
        raise NotImplementedError(
            "evolve_imaginary_time is implemented for itensor_version=3 and "
            "'python' (got itensor_version=%r); use mode='ED' for a small "
            "chain" % (self.itensor_version,))
    if getattr(self, "conserved_sector", None):
        raise ValueError("conserved-sector mode (%s): evolve_imaginary_time "
                "works on dense sites only; clear the sector with "
                "set_conserved_sector()" % (self.conserved_sector,))
    if self.itensor_version == 3 and self.ns < 3:
        raise NotImplementedError("evolve_imaginary_time with "
                "itensor_version=3 needs at least 3 sites; use "
                "itensor_version='python' or mode='ED'")
    if not hasattr(getattr(self, "_session", None), "implicit_step"):
        raise NotImplementedError("evolve_imaginary_time needs the "
                "itensor_version=%r session, which is not available here "
                "(is the extension compiled?); use itensor_version='python' "
                "or mode='ED'" % (self.itensor_version,))
    if not self.is_hermitian(self.hamiltonian):
        raise ValueError("evolve_imaginary_time needs a Hermitian "
                "Hamiltonian: A-stability of the implicit steps rests on a "
                "real spectrum bounded below")


def _ground_energy(self):
    """The lowest eigenvalue of H, the shift every factor needs. The
    chain's own e0 is that only when a plain ground-state solve produced
    it: after set_gs()/set_initial_wf() it is the energy of the state the
    caller set, and after gs_energy_generalized() that of the generalized
    eigenvector, both above E0, and a shift above E0 puts the pole of the
    trapezoid factor inside the spectrum once dtau (E_s - E0) reaches 2.
    In those cases the energy comes from a solve on a clone, which leaves
    the chain's state as it is."""
    from . import groundstate
    if groundstate.state_supplied(self) or \
            getattr(self, "lam_generalized", None) is not None:
        return float(np.real(self.lowest_eigenvalue(self.hamiltonian)))
    return float(np.real(self.gs_energy()))


def _shift(self, shift):
    """(value, from_energy) for the session: the ground-state energy for
    "ground", the running energy for "energy", or a number."""
    if isinstance(shift, str):
        if shift == "ground": return _ground_energy(self), False
        if shift == "energy": return 0., True
        raise ValueError("shift must be 'ground', 'energy' or a number "
                         "(got %r)" % (shift,))
    return float(shift), False


def _ed(self, wf, steps, operators):
    """The exact normalized e^{-tau H} wf at every grid time, by
    expm_multiply, and the operators measured there. H is shifted by its
    ground-state energy, so that e^{-dtau (H-E0)} <= 1 does not overflow
    on a long step."""
    import scipy.sparse as sparse
    from scipy.sparse.linalg import expm_multiply
    from .edtk.edchain import State
    ed = self.get_ED_obj()
    Hm = sparse.csr_matrix(ed.MO2matrix(self.hamiltonian))
    Hm = Hm - float(self.gs_energy(mode="ED"))*sparse.identity(Hm.shape[0],
                                                               format="csr")
    Om = [ed.MO2matrix(op) for op in operators]
    v = np.array(wf.v, dtype=complex)
    v = v/np.linalg.norm(v)
    vals = [[np.vdot(v, A@v) for A in Om]]
    for dt in steps:
        v = expm_multiply(-dt*Hm, v)
        v = v/np.linalg.norm(v)
        vals.append([np.vdot(v, A@v) for A in Om])
    return State(v, wf.MBO), np.array(vals)


def evolve(self, wf, tau, operators=None, nsteps=4, tau0=None, order=2,
           shift="ground", nsweeps=1, tol=1e-10, krylov=10, maxiter=100,
           mode="DMRG"):
    """Return (taus, values, state): the grid times, the operators measured
    on the normalized state at each of them (shape (len(taus), number of
    operators), empty without operators), and the normalized
    e^{-tau H} wf. Keywords as in Many_Body_Chain.evolve_imaginary_time."""
    from .mpsalgebra import wavefunction_mode
    from . import mps
    operators = list(operators) if operators is not None else []
    steps = _plan(self, tau, tau0, nsteps, order)
    taus = np.concatenate([[0.], np.cumsum(steps)])
    if mode == "ED" or wavefunction_mode(wf) == "ED":
        if wavefunction_mode(wf) != "ED":
            raise ValueError("evolve_imaginary_time(mode='ED') needs an ED "
                             "state (get_gs(mode='ED') or random_state(mode='ED'))")
        state, vals = _ed(self, wf, steps, operators)
        return taus, vals, state
    if mode != "DMRG":
        raise ValueError("mode must be 'DMRG' or 'ED' (got %r)" % (mode,))
    _check_dmrg(self)
    value, from_energy = _shift(self, shift)
    from .groundstate import send_hamiltonian
    session = self._session
    session.set_sweep_params(self.maxm, self.nsweeps, self.cutoff, self.noise)
    session.set_mpomaxm(max(self.maxm, self.mpomaxm))
    send_hamiltonian(self)
    norm = wf.dot(wf).real
    if norm <= 0.:
        raise ValueError("evolve_imaginary_time: the initial state is zero")
    psi = (wf*(1./np.sqrt(norm))).cpp_handle
    def measure(handle):
        state = mps.MPS(self, cpp_handle=handle)
        return [state.aMb(op, state) for op in operators]
    vals = [measure(psi)]
    for dt in steps:
        psi = session.implicit_step(psi, step_factors(dt, order), value,
                                    from_energy, int(nsweeps), float(tol),
                                    int(krylov), int(maxiter))
        if operators: vals.append(measure(psi))
    state = mps.MPS(self, cpp_handle=psi).copy()
    vals = np.array(vals) if operators else np.zeros((len(taus), 0))
    return taus, vals, state


_PLAN_KEYWORDS = ("nsteps", "tau0", "order", "shift", "nsweeps", "tol",
                  "krylov", "maxiter")


def _implicit_available(self):
    """True when the implicit steps can run on this chain: a v3 or
    "python" session, no conserved sector, at least 3 sites on v3, and a
    Hermitian Hamiltonian."""
    if self.itensor_version not in (3, "python"): return False
    if getattr(self, "conserved_sector", None): return False
    if self.itensor_version == 3 and self.ns < 3: return False
    if not hasattr(getattr(self, "_session", None), "implicit_step"):
        return False
    return bool(self.is_hermitian(self.hamiltonian))


def metts_plan(self, T, itime=None, itime_options=None, dbeta_half_step=0.05):
    """The imaginary-time evolution a METTS sample takes to beta/2: None
    for TDVP steps of dbeta_half_step, or the implicit plan the sessions
    run (a list of betas per step of the logarithmic grid from 0 to
    1/(2T), and the solver keywords), with itime_options the keywords of
    evolve_imaginary_time.

    The implicit steps default to nsteps=8 here, against 4 in
    evolve_imaginary_time, so that their time-step error is that of the
    TDVP default. Measured on a 12-site Heisenberg chain from the Neel
    state against ED at tau = beta/2 (itensor_version=3): TDVP with
    dbeta_half_step=0.05 is off by 2.1e-5, 1.7e-6 and 5.9e-8 in the energy
    at beta/2 = 5, 10, 25, and the trapezoid rule with nsteps=8 by 2.0e-5,
    1.4e-5 and 4.6e-8, in 1.0, 1.2 and 1.6 s against 1.4, 2.8 and 7.3 s;
    with nsteps=4 it is 5 to 30 times less accurate at intermediate beta.

    itime=None (the default) takes the implicit steps whenever they can
    run here and their grid needs fewer linear solves than TDVP needs
    steps, and TDVP otherwise. One solve costs about one TDVP step, so
    this is the cheaper of the two at matched accuracy: TDVP at high T,
    where its few steps beat the panels of the grid, the implicit steps
    below T of about 0.12 in units of the largest coupling, where their
    number grows as log(beta) against beta."""
    if itime not in (None, "tdvp", "implicit"):
        raise ValueError("itime must be None, 'tdvp' or 'implicit' (got %r)"
                         % (itime,))
    opts = dict(itime_options or {})
    if itime == "tdvp":
        if opts:
            raise ValueError("itime_options= is not used with itime='tdvp'")
        return None
    unknown = sorted(set(opts) - set(_PLAN_KEYWORDS))
    if unknown:
        raise TypeError("itime_options: unknown keyword(s) %s; the known "
                        "ones are %s" % (unknown, list(_PLAN_KEYWORDS)))
    if itime is None and not _implicit_available(self):
        return None
    if itime == "implicit":
        _check_dmrg(self)
    order = opts.get("order", 2)
    beta_half = 1./(2.*T)
    steps = _plan(self, beta_half, opts.get("tau0"), opts.get("nsteps", 8),
                  order)
    if itime is None:
        ntdvp = max(1, int(np.ceil(beta_half/dbeta_half_step)))
        if len(steps)*(order//2) >= ntdvp: return None
    value, from_energy = _shift(self, opts.get("shift", "ground"))
    return dict(steps=[step_factors(dt, order) for dt in steps],
                shift=value, shift_energy=from_energy,
                nsweeps=int(opts.get("nsweeps", 1)),
                tol=float(opts.get("tol", 1e-10)),
                krylov=int(opts.get("krylov", 10)),
                maxiter=int(opts.get("maxiter", 100)))


def session_plan_kwargs(self, plan):
    """The keywords a session's metts_vev/metts_dynamical_correlator take
    for a plan: implicit= on "python", the itime_* ones on v3."""
    if plan is None: return {}
    if self.itensor_version == "python": return dict(implicit=plan)
    return dict(itime_steps=plan["steps"], itime_shift=plan["shift"],
                itime_shift_energy=plan["shift_energy"],
                itime_nsweeps=plan["nsweeps"], itime_tol=plan["tol"],
                itime_krylov=plan["krylov"], itime_maxiter=plan["maxiter"])
