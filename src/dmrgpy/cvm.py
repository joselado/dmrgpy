import contextlib
import numpy as np
from . import operatornames
from . import multioperator

def dynamical_correlator(self,es=np.linspace(0.,10.0,100),
        delta=1e-1,name="XX",i=0,j=0):
    """
    Compute the dynamical correlator using CVM method in DMRG.

    Everything that does not depend on the frequency is hoisted out of
    the loop over es (operator resolution, the ground state, and the
    right-hand side b = -eta*B|GS>, which cvm_correction_vector's linear
    system shares across all frequencies). Each frequency point is
    solved independently from the same cold start -- warm-starting a
    point's CG from the previous point's correction vector was tried and
    measured to be actively harmful (on a 12-site Heisenberg chain it
    left several points stagnating at a ~100x worse residual than the
    cold start reached, shifting the correlator by ~50%), so it is
    deliberately not done.
    """
    AB = operatornames.str2MO(self,name,i=i,j=j) # resolve operators once
    A,B = AB[0],AB[1]
    if _use_ddmrg(self,A,B): # variational solver (pyitensor only)
        return dynamical_correlator_ddmrg(self,A,B,es=es,delta=delta)
    wf0 = self.get_gs() # computes or returns the cached GS; also sets self.e0
    # sweep parameters used both by B*wf0 below and by every CG solve
    with _cvm_sweep_params(self):
        b = (-delta)*(B*wf0) # -eta*B|GS>, identical for every frequency
        out = [] # empty list
        for e in es: # loop over energies
            o,_,nit,res = cvm_correction_vector(self,A,B,e,delta,
                    tol=self.cvm_tol,max_it=int(self.cvm_nit),b=b)
            print("CVM in E = ",e," iterations = ",nit," residual = ",res)
            out.append(o) # store
    out = np.array(out)
#    from .inference import points2function
#    (es,out) = points2function(es,out)
    return (es,out) # return result




def _use_ddmrg(self,A,B):
    """Whether to route this correlator through the variational
    (Jeckelmann DDMRG) solver instead of the global conjugate gradient.

    Three conditions, all necessary:

    * `self.cvm_solver == "variational"` -- opt-in, so the default path
      stays byte-identical on every backend (see manybodychain.py).
    * `itensor_version == "python"` -- pyitensor is where the two-site
      environment machinery the solver reuses lives; the compiled
      backends would need it written again in C++.
    * A == B^dagger -- the variational principle needs a real quadratic
      form bounded below, i.e. the same operator on both sides. A general
      off-diagonal <GS|A (z-H)^-1 B|GS> has no such functional, so it
      falls back to CG rather than silently minimizing the wrong thing.
    """
    if getattr(self,"cvm_solver","cg")!="variational": return False
    if self.itensor_version!="python": return False
    # Unlike the CG path below, which only ever *applies* A and B (so the
    # already-built operators toMPO() returns work there unchanged), this
    # solver hands their symbolic term lists to ddmrg_correction_vector.
    # Checked here rather than in dynamical_correlator_ddmrg because the
    # A==B^dagger test just below already needs operator algebra a
    # StaticOperator does not have (operator_norm calls .simplify()).
    operatornames.require_symbolic(A,"cvm_solver='variational'")
    operatornames.require_symbolic(B,"cvm_solver='variational'")
    # Symbolically first, as the TD route decides the same question
    # (timedependent._pair_is_self_adjoint): the canonical form proves an
    # adjoint pair at any scale and refuses a name whose adjoint it does
    # not know. Only when that proof does not land (an alias pair such as
    # Sx against (Sp+Sm)/2, a same-site identity) does the numerical test
    # run, and is_zero_operator is relative to the difference's own
    # largest coefficient. It used to be the numerical test alone, with an
    # absolute 1e-4 on a squared norm, so any pair of operators below
    # about 1e-2 in size came through as adjoint and this solver, which
    # reads only B, returned C[B^dagger,B] for C[A,B]: 2.124 of the peak
    # off for (1e-2*Sz0, 1e-2*Sz3) (2026-09-25b audit, finding 16). "Not
    # adjoint" is always safe here, since it only sends the call to CG.
    from .multioperatortk import canonical
    if canonical.is_dagger_pair(A,B): return True
    return self.is_zero_operator(A.get_dagger()-B)


def dynamical_correlator_ddmrg(self,A,B,es=None,delta=1e-1):
    """Frequency sweep using the variational correction vector.

    Each point seeds its sweep with the previous point's converged
    correction vector. That is the opposite of what the CG path does, and
    deliberately so: warm-starting CG was measured to be actively harmful
    (see dynamical_correlator's docstring), because a global CG inherits
    the whole error history of whatever it starts from, while a
    variational sweep is non-increasing in W from *any* start -- a closer
    start can only reduce the number of sweeps needed.
    """
    self.get_gs() # ensure the ground state (and self.e0) exist
    out,x = [],None
    with _cvm_sweep_params(self):
        for e in es:
            o,x = self._session.ddmrg_correction_vector(
                    A.to_terms(),B.to_terms(),e,delta,self.e0,
                    int(self.cvm_maxm),int(self.cvm_nsweeps),x0=x)
            print("DDMRG in E = ",e," value = ",o)
            out.append(o)
    return (es,np.array(out))

@contextlib.contextmanager
def _cvm_sweep_params(self):
    """Point every MPO application/MPS truncation at cvm_maxm rather than
    the DMRG maxm, for the duration of a CVM computation. The C++
    backends do this via self._session.set_sweep_params(...) (no
    _session object exists to call on julia_live), so this temporarily
    overrides self.maxm instead -- mpsjulialive/mpo.py's MPO.__mul__ and
    mps.py's MPS.__add__ both read self.MBO.maxm directly, no other
    channel to override the bond dimension exists for that backend.
    A context manager (not a plain set-and-return-the-old-value function)
    so self.maxm is restored via `finally` even if the CG solve inside
    raises -- confirmed directly that the previous plain-function version
    left self.maxm permanently stuck at cvm_maxm on any exception, since
    its restore was a bare statement after the frequency loop. Safe to
    nest: cvm_correction_vector wraps itself in this too (so it stays
    correct when called standalone, not just from dynamical_correlator's
    loop), and each nested entry/exit only ever restores the self.maxm
    value that was live when *it* was entered."""
    maxm0 = self.maxm
    if self.itensor_version=="julia_live":
        self.maxm = self.cvm_maxm
    else:
        self._session.set_sweep_params(self.cvm_maxm,self.nsweeps,self.cutoff,self.noise)
        self._session.set_verbose(self.verbose)
        self._session.set_mpomaxm(max(self.cvm_maxm,self.mpomaxm))
    try:
        yield
    finally:
        if self.itensor_version=="julia_live":
            self.maxm = maxm0










def cvm_correction_vector(self,A,B,omega,eta,tol=1e-5,max_it=1000,
        b=None):
    """
    Correction Vector Method (Ramasesha; Kuhner & White 1999), returning
    dynamics.py's house convention -- the complex Lehmann density
    i*(G^R-G^A)/(2*pi) of the pair (A,B) -- computed by solving the
    Hermitian, positive-definite system

        [(H-omega-E0)^2 + eta^2] xc = -eta * B|GS> =: b

    by conjugate gradient, then recovering the correction vector as

        x = i*xc + (H-omega-E0)/eta * xc = (omega+E0+i*eta-H)^{-1} B|GS>

    This is the system the ED backend solves too
    (edtk/dynamics.py::solve_cv), carried out with MPS algebra
    (StaticOperator/MPS arithmetic) instead of dense matrices. Because
    (H-omega-E0)^2+eta^2 is Hermitian positive-definite whenever H is
    Hermitian, a Krylov method on it has no breakdown mode -- unlike the
    previous direct solve of the non-Hermitian (z-H) system via a
    hand-rolled BiCGSTAB (Chain::cvm_dynamical_correlator/bicstab in
    chain_session.h), which had no protection against BiCGSTAB's
    near-singular breakdown and could blow up.

    Stopping rule: ||r|| <= tol*||b||, relative to the right-hand side.
    It used to be ||r|| <= tol, absolute, while b carries the units of
    eta and of B, so wherever eta*||B|GS>|| was at or below tol (a small
    operator, a Hamiltonian in small units, or eta <= 2e-5 at unit scale,
    which includes get_kondo_spectrum's documented default delta=2e-6)
    the start already passed and every frequency returned the flat
    eta*<AB>/pi: 0.887 of the peak off at an operator scale of 1e-4, and
    a second-order Kondo dI/dV of 9.12e-09 against an exact 4.7124
    (2026-09-25b audit, finding 13). With b -> s*b and the system matrix
    -> s^2*(...) the relative residual is unchanged, and it bounds the
    answer: |dC|/(peak scale ||A|GS>|| ||B|GS>||/(pi*eta)) <= ||r||/||b||.
    The start is xc = 0 (r = b), not the xc = b the old BiCGSTAB used:
    b has the units of eta*B and the solution those of B/eta, so starting
    from b made the whole iteration depend on the units; from 0 it is
    exactly covariant.

    The best iterate and the two early exits read the CG functional
        phi(x) = <x|M|x>/2 - Re<b|x> = -Re(<x|b> + <x|r>)/2,
    M the system matrix and r = b - M x (the recurrence residual), and
    not the residual 2-norm they used to read. CG minimizes phi over a
    growing Krylov space, so in exact arithmetic phi falls at every
    iteration, by alpha*||r||^2/2, while ||r|| does not: on this system it
    rises by up to about the square root of the condition number, of order
    (level spacing/eta)^2, before it falls. Measured with float64 dense CG
    on a 6-site chain at eta=2e-3 on a line: 217 times its running
    minimum, and the MPS CG at full bond dimension tracked it digit for
    digit. So both exits fired on exact CG, and since the residual had
    never gone below its initial value, the "best iterate by ||r||" they
    returned was the initial guess: 1.5915e-04 against 13.2636 at
    eta=2e-3, and 7 of 121 points up to 0.989 of the peak off on a
    10-site chain at the test suite's own eta=0.05 (finding 14). Read off
    phi, "diverging" cannot happen on an untruncated solve and "no
    progress" only once phi has converged to its own rounding, the best
    iterate of one is the latest, and the first step always improves on
    the start. phi costs two MPS dot products per iteration.

    Why not the conjugate residual method, whose ||r|| is monotone too:
    it was tried and measured in the regime the exits exist for, and it
    loses them there. It carries M*p by recurrence as well as r, and
    under truncation its recurrence residual keeps falling while the
    iterate does not improve: on the 20-site Heisenberg chain at
    cvm_maxm=30 it ran 952 and 1000 iterations (228 s and 240 s per point
    against 11 s) and reported convergence at 1e-5 of ||b||, and on 14
    sites against exact ED it was no more accurate than CG (4.2e-2 against
    3.3e-2 at omega=0.3, eta=0.15, cvm_maxm=10) at up to three times the
    iterations. Truncated CG's recurrence does the opposite -- phi stops
    falling and ||r|| rises -- which is what the exits catch: on the same
    20-site chain (omega=0.3 and 1.0) this loop stops at 104 and 112
    iterations (41 s and 45 s single-core, against 50 iterations and 15 s
    and 21 s for the old exits, which returned the flat initial guess
    there), with a relative residual of 7.1 and 5.4, which the warning
    reports.

    b optionally supplies a precomputed right-hand side -eta*B|GS>,
    which is frequency-independent and can be shared across a sweep.
    Returns (value, best_xc, iterations_used, best_residual) so callers
    can report the solver effort and convergence quality per point;
    best_residual is the absolute ||r|| of best_xc, whose ratio to
    ||b|| is what tol is compared with. Implemented purely with
    already-exposed Python primitives (self.toMPO()/StaticOperator for a
    build-once/apply-many MPO, MPS +/-/scalar-* and .dot() for the rest),
    so no new C++/pybind11 bindings and no recompilation are needed to
    tune this further -- unlike the KPM/TDVP paths, which run their inner
    loop in C++.
    """
    wf0 = self.get_gs() # ground state (cheap/cached; also sets self.e0)
    with _cvm_sweep_params(self):
        # (H-omega-E0) as an MPO, rebuilt per omega: measured at ~1 ms, i.e.
        # negligible next to a single CG iteration (~2 MPO applications,
        # ~0.03-0.1 s at 12-20 sites). Applying the shift at MPS level
        # instead (HmE0*v - omega*v, one shared MPO) was benchmarked ~2x
        # SLOWER per application (extra truncated sums), so don't "optimize"
        # this line by hoisting it out of the frequency loop.
        Hshift = self.toMPO(self.hamiltonian-(self.e0+omega))
        def applyA(v): return Hshift*(Hshift*v) + (eta*eta)*v # (H-omega-E0)^2+eta^2
        if b is None: b = (-eta)*(B*wf0) # -eta*B|GS>
        bnorm = np.sqrt(abs(b.dot(b).real)) # ||b||, the scale of every residual
        target = tol*bnorm # relative stop; b = 0 stops at once on xc = 0
        xc = None # xc = 0, so r = b (covariant under a change of units)
        r = b
        p = r
        rs_old = bnorm**2 # ||r||^2, real since the system is Hermitian PD
        best_xc,best_res,best_phi = None,bnorm,0. # phi(0) = 0
        # Early-termination guards. The MPS truncation (maxdim=cvm_maxm)
        # puts a floor on what the recurrence can reach, and once it hits
        # it, it diverges (traced on a 20-site Heisenberg chain at
        # cvm_maxm=30: the running residual grew monotonically to ~3e4 by
        # iteration 1000), so continuing to iterate is waste. Stop (a) after
        # `patience` iterations without a new minimum of phi, or (b) once
        # phi is above its minimum and the running residual `blowup` times
        # above the best iterate's. Exact CG lowers phi at every iteration,
        # so neither can fire on an untruncated solve until phi has
        # converged to its own rounding (`slack`, far below any cvm_tol).
        # They used to read ||r|| alone -- (a) counted iterations without a
        # 0.1% gain in it -- and fired on exact CG, whose ||r|| is not
        # monotone (finding 14; see the docstring). The best iterate is the
        # latest of the lowest phi, ties within the rounding going to the
        # latest; a solve that reaches tol returns the iterate that did.
        # Tunable as self.cvm_patience / self.cvm_blowup (manybodychain.py).
        patience = int(getattr(self,'cvm_patience',50))
        blowup = float(getattr(self,'cvm_blowup',100.0))
        slack = 1e-12 # relative rounding of phi, which sums two dot products
        since_best = 0 # iterations since the last new minimum of phi
        niter = 0
        for k in range(max_it):
            if best_res<=target: break # b = 0 (xc = 0 is exact), or converged
            Ap = applyA(p)
            alpha = rs_old/p.dot(Ap).real # real: <p|A|p> is real for Hermitian A
            xc = alpha*p if xc is None else xc + alpha*p
            r = r - alpha*Ap
            rs_new = r.dot(r).real
            res = np.sqrt(abs(rs_new))
            niter = k+1
            phi = -0.5*(xc.dot(b) + xc.dot(r)).real # the CG functional
            if res<=target: # converged: certified by its own residual
                best_xc,best_res,best_phi = xc,res,phi
                break
            s = slack*abs(best_phi)
            if phi<best_phi-s: since_best = 0 # progress
            else: since_best += 1
            if phi<=best_phi+s: best_xc,best_res,best_phi = xc,res,phi
            if since_best>=patience: break # the recurrence stopped improving
            if phi>best_phi and res>blowup*best_res: break # diverging
            p = r + (rs_new/rs_old)*p
            rs_old = rs_new
        _warn_if_unconverged(self,omega,best_res,tol,bnorm)
        # The house convention (dynamics.py) needs both resolvents,
        #   i*(G^R-G^A)/(2*pi),   G^R/G^A = <GS|A (w+E0 -+ ... i*eta-H)^-1 B|GS>,
        # and both come out of this one solve for free. The system
        # matrix (H-w-E0)^2+eta^2 is even in eta, so flipping eta only
        # negates the right-hand side: xc(-eta) = -xc(+eta), and the
        # recovery formula turns into
        #   x(-eta) = -i*xc + (H-w-E0)*xc/eta  (vs  +i*xc + ... for +eta).
        # Hence G^R-G^A = <GS|A (2i*xc)> and the whole correlator is just
        #   i*(2i*<GS|A|xc>)/(2*pi) = -<GS|A|xc>/pi,
        # with the (H-w-E0)*xc/eta piece -- the dispersive part -- cancelling
        # identically. So this is also one MPO application cheaper than the
        # -Im(G^R)/pi it replaces, which built the full correction vector
        # x = i*xc + Hshift*xc/eta first.
        #
        # That -Im(G^R)/pi equals this only for a real Lehmann weight
        # M_n = <GS|A|n><n|B|GS>, i.e. only for A = B^dagger. For a general
        # pair the two differ by more than the correlator's own peak
        # (measured 0.5746 against a peak of 0.4257 on a 4-site
        # complex-hopping fermionic chain), and it was this backend, not
        # mode="ED", that sat off the convention: only the density
        # satisfies integral dw = <GS|A B|GS>. (2026-09 audit, finding #5.)
        if best_xc is None: # never left xc = 0 (b = 0, or max_it = 0)
            return 0., 0.*b, niter, best_res
        C = -wf0.dot(A*best_xc)/np.pi # <GS|A|xc>, xc = -eta*[(H-w-E0)^2+eta^2]^-1 B|GS>
        return C, best_xc, niter, best_res


_UNCONVERGED_WARNED = set()

def _warn_if_unconverged(self,omega,best_res,tol,bnorm,factor=100.):
    """Say so, loudly, when the solve stopped far short of `tol`.

    Relative to ||b||, as the stopping rule is: best_res/||b|| is also a
    bound on the error of the returned value relative to the spectrum's
    peak scale (see cvm_correction_vector's docstring). It used to compare
    the absolute residual with 100*tol, so it was silent exactly where the
    absolute stop had failed -- a residual that was never reduced but sat
    below 1e-3 because b itself did (2026-09-25b audit, findings 13 and
    14) -- and it named the MPS-truncation floor as the cause when the
    early exits had fired on an untruncated solve.

    This is not a cosmetic nicety. Whenever the MPS truncation actually
    binds -- `cvm_maxm` below the bond dimension the correction vector
    needs -- the recurrence stops improving, and the returned best iterate
    can be far from converged. Measured on a 14-site Heisenberg chain
    against an exact ED correction vector with the residual-keyed exits
    this module used before, the returned spectrum was wrong by 1.3e-1
    (against peak values of ~1.5e-1) at cvm_maxm = 10, 20 and 40 alike,
    at eta*<Sz Sz>/pi, the value of that loop's initial guess. Before
    this warning existed, that came back looking like an ordinary answer.

    Warned once per (chain, tolerance) rather than per frequency, so a
    300-point sweep does not print 300 identical paragraphs.
    """
    if not (best_res>factor*tol*bnorm): return
    key = (id(self),float(tol))
    if key in _UNCONVERGED_WARNED: return
    _UNCONVERGED_WARNED.add(key)
    import warnings
    warnings.warn(
        "CVM: the correction vector did not converge (relative residual "
        "%.3g at omega=%.4g, requested %.3g), and the returned spectrum "
        "can be off by up to about that fraction of its peak scale. If "
        "cvm_maxm is below the bond dimension the correction vector needs, "
        "the MPS truncation floors the solve: raise cvm_maxm, or -- on "
        "itensor_version='python' with a diagonal correlator -- set "
        "chain.cvm_solver='variational' to use the variational (Jeckelmann "
        "DDMRG) solver instead, which truncates inside the ansatz (see "
        "src/dmrgpy/pyitensor/ddmrg.py). If the residual was still falling, "
        "raise cvm_nit or cvm_patience."%(best_res/bnorm,omega,tol),
        RuntimeWarning,stacklevel=3)




def dynamical_correlator_analytic_continuation(self,name=None,
        delta=1e-1,es=np.linspace(0.,5.0,300)):
    """
    Compute the dynamical correlator using analytic continuation
    """
    A,B = name[0],name[1]
    wf = self.get_gs() # get the ground state
    wfa = A.get_dagger()*wf # apply A to the GS
    wfb = B*wf # apply B to the GS
    e0 = self.gs_energy() # ground state energy
    Hp = self.hamiltonian - e0
    def f(e): # function to compute
        wfi = self.applyinverse(-self.hamiltonian+(e0+e),wfa)
        return wfb.dot(wfi) # return result
#    return es,-np.array([f(e+1j*delta*10) for e in es]).imag*2/np.pi # brute force
    from .analyticcontinuation import imag2real
    xz = es*1j
    xz = np.linspace(delta*10,10.,100)*1j
    xz = np.concatenate([-xz,xz])
    xz = np.linspace(min(es),max(es),20) + delta*40*1j
#    xz = [np.random.random()-.5+1j*np.random.random()+0.5j for i in range(40)]
#    xz = 40.*np.array(xz)
    outz = np.array([f(z) for z in xz]) # complex axis
    esz,out = imag2real(xz,outz,x=es+1j*delta)
    out = -out.imag*2/np.pi
    return es,out



from .nonhermitian.dynamics import dynamical_correlator_cvm_explicit

#def dynamical_correlator_cvm_explicit(self,name=None,
#        delta=1e-1,es=np.linspace(0.,5.0,300)):
#    """
#    Compute the dynamical correlator using analytic continuation
#    """
#    ### So far this just works for onsite correlators
#    A,B = name[0],name[1]
#    if not self.is_zero_operator(A.get_dagger()-B): 
#        print("Only implemented for A^\dagger=B")
#        raise
#    wf = self.get_gs() # get the ground state
#    wfa = A.get_dagger()*wf # apply A to the GS
#    wfb = B*wf # apply B to the GS
#    e0 = self.gs_energy() # ground state energy
#    Hp = self.hamiltonian - e0
#    def f(e,delta): # function to compute
#        wfi = self.applyinverse(-self.hamiltonian+(e0+e+1j*delta),wfa)
#        return wfb.dot(wfi) # return result
#    from .analyticcontinuation import imag2real
#    outz = np.array([f(z,delta) - f(z,-delta) for z in es]) # complex axis
#    return es,1j*outz/np.pi





