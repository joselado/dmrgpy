# this library contains routines to perform thermal calculations
import numpy as np

from .spinchain import Spin_Chain
from . import multioperator

class Thermal_Spin_Chain():
    def __init__(self,sites,T=0.1,mode=None,**kwargs):
        if T<0.: # a negative temperature used to be silently treated as T=0
            raise ValueError("Thermal_Spin_Chain: T must be >= 0, got "
                             +repr(T))
        # mode is the wrapper's own, and get_gs() writes it onto MBChain;
        # that used to be a hardcoded "DMRG", which overwrote a forwarded
        # mode="ED" so that Thermal_Spin_Chain(...,mode="ED") ran DMRG. The
        # default is None, the chain's own default, rather than "DMRG": a
        # chain whose mode was "DMRG" used to answer an explicit mode="ED"
        # read of MBChain by DMRG (2026-09-25b hole hunt, finding 3)
        if mode is not None: # None leaves the chain on its default, as on any chain
            from .mode import _check_mode
            _check_mode(mode,"Thermal_Spin_Chain mode (mode=)")
        # every other keyword is a setting of MBChain, checked here under
        # this class's own name (the Spin_Chain below checks it again), so
        # a misspelled one is reported as Thermal_Spin_Chain()'s, not as
        # the Spin_Chain() the caller never wrote
        from .sites import check_settings
        check_settings(self,{k:v for k,v in kwargs.items()
                             if k!="itensor_version"})
        sitesT = []
        for s in sites: 
            sitesT += [s]
            sitesT += [s]
        n = len(sites) # number of sites
        # **kwargs was accepted and then dropped on the floor here, so
        # Thermal_Spin_Chain(...,itensor_version="python") silently built
        # the default (compiled) backend instead
        self.MBChain = Spin_Chain(sitesT,mode=mode,**kwargs) # get the chain
        self.computed_gs = False
        self.T = T # temperature
        Sx = [self.MBChain.Sx[2*i] for i in range(n)] 
        Sy = [self.MBChain.Sy[2*i] for i in range(n)] 
        Sz = [self.MBChain.Sz[2*i] for i in range(n)] 
        self.all_Sx = self.MBChain.Sx
        self.all_Sy = self.MBChain.Sy
        self.all_Sz = self.MBChain.Sz
        self.Sx = Sx
        self.Sy = Sy
        self.Sz = Sz
        self.wf0 = None
        self.mode = mode
    def get_gs(self):
        """Compute the ground state"""
        if self.computed_gs: 
            return self.wf0
        else:
            def terms():
                for i in range(len(self.Sx)):
                    yield self.all_Sx[2*i]*self.all_Sx[2*i+1]
                    yield self.all_Sy[2*i]*self.all_Sy[2*i+1]
                    yield self.all_Sz[2*i]*self.all_Sz[2*i+1]
            h = multioperator.msum(terms()) # initialize
            self.MBChain.mode = self.mode # overwrite the mode
#            wf = self.MBChain.random_mps()
            # float() so that a numpy scalar T (an np.linspace element)
            # takes exactly the path a Python float does
            T = float(self.T)
            if not T>=0.: # the setter bypasses __init__'s check
                raise ValueError("Thermal_Spin_Chain: T must be >= 0, got "
                                 +repr(self.T))
            # The switch to the plain ground state is in units of the
            # spectral width W, which anneal() needs anyway: it used to be
            # an absolute T>1e-5, so a Hamiltonian written below about
            # 1e-5 in absolute units got its ground state at every T
            method = self.anneal_method # checked before any solve
            if method not in (None,"taylor","implicit"):
                raise ValueError("Thermal_Spin_Chain: anneal_method must be "
                                 "None, 'taylor' or 'implicit', got "
                                 +repr(method))
            if method!="implicit" and self.anneal_options:
                raise ValueError("Thermal_Spin_Chain: anneal_options is "
                                 "used only with anneal_method='implicit'")
            erange = None
            if 0.<T<np.inf: erange = band_edges(self.MBChain,self.hamiltonian)
            if T>0. and (erange is None or T>1e-5*(erange[1]-erange[0])):
                self.MBChain.set_hamiltonian(h) # singlet Hamiltonian
                wf = self.MBChain.get_gs() # get the fully entangled WF
                if method is None: method = auto_anneal_method(self,wf,T,erange)
                if method=="implicit":
                    wf0 = anneal_implicit(self.MBChain,self.hamiltonian,wf,T,
                                          erange=erange,
                                          options=self.anneal_options)
                else:
                    wf0 = anneal(self.MBChain,self.hamiltonian,wf,T,
                                 step=self.anneal_step,
                                 order=self.anneal_order,erange=erange)
                wf0 = wf0.normalize()
                # Through the public setters, so MBChain holds the
                # physical Hamiltonian and the annealed state as a state
                # set by hand: the next read takes it unswept, and every
                # correlator measures it. Assigning MBChain.wf0 and
                # MBChain.hamiltonian directly left the singlet
                # Hamiltonian on the session and in the send-cache (and,
                # on mode="ED", the singlet ED object), so gs_energy()
                # returned the singlet energy (-2.25 against <wf|H|wf> =
                # -0.4164 on a 3-site chain at T=1), a correlator re-solved
                # over the annealed state and vev() on ED read the singlet
                # state (<Sz0 Sz1> -0.1667 and 0.0 against -0.0694)
                self.MBChain.set_hamiltonian(self.hamiltonian)
                self.MBChain.set_gs(wf0)
            else: # T<=1e-5*W: plain ground state. `else`, not another elif:
                # the two branches used to leave T exactly 1e-5 (and, before
                # the check in __init__, any negative T) falling through
                # both, so wf0 was never assigned and the next line raised
                # UnboundLocalError. MBChain's own solved state is this
                # one already, so nothing is written back to it
                self.MBChain.set_hamiltonian(self.hamiltonian)
                wf0 = self.MBChain.get_gs() # get the fully entangled WF
                wf0 = wf0.normalize()
            self.computed_gs = True
            self.wf0 = wf0 # store
            return wf0
    # The imaginary-time stepper get_gs() hands to anneal(): the largest
    # dimensionless step dtau*W (W the spectral width of the Hamiltonian)
    # and the Taylor order of each step. Class attributes, so that an
    # instance overrides them the way it sets T (tc.anneal_step = 1.0)
    # without a constructor keyword. See anneal() for what they cost and
    # buy: at these defaults the stepper's own thermal energy of the open
    # S=1/2 Heisenberg chain is within 7e-7 (relative) of the Boltzmann
    # value at every n from 3 to 12 and T from 0.05 to 10 (its closed form
    # over the exact spectrum), for about 2*beta*W factor applications.
    anneal_step = 2.0
    anneal_order = 8
    # anneal_method="implicit" takes the implicit steps on a logarithmic
    # grid instead (anneal_implicit(), arXiv:2606.02930), on
    # itensor_version=3 and "python", and on mode="ED" as the exact
    # e^{-tau H} on the same grid; anneal_options takes the keywords of
    # evolve_imaginary_time but the shift, which is the lower band edge.
    # anneal_step and anneal_order belong to the Taylor stepper only.
    # None, the default, picks one of the two (auto_anneal_method()).
    anneal_method = None
    anneal_options = None
    def set_hamiltonian(self,h):
        self.hamiltonian = (h + h.get_dagger())/2.
        self.computed_gs = False

                

                
def band_edges(sc,h):
    """Lowest and highest eigenvalue of the Hermitian operator h, from two
    ground-state solves on a clone of sc (its own backend and mode), as
    Many_Body_Chain.bandwidth() does"""
    mbc = sc.clone() # clone the object
    mbc.set_hamiltonian(h) ; emin = mbc.gs_energy()
    mbc.set_hamiltonian(-1*h) ; emax = -mbc.gs_energy()
    emin, emax = float(np.real(emin)), float(np.real(emax))
    return emin, max(emin,emax)


def taylor_exp_roots(order):
    """Roots r_j of the order-k Taylor polynomial of exp(-x),
    p_k(x) = sum_{j<=k} (-x)^j/j!, so that p_k(x) = prod_j (1-x/r_j)
    (p_k(0) = 1 fixes the normalization)"""
    from math import factorial
    order = int(order)
    if order<1: raise ValueError("anneal: order must be >= 1, got "+repr(order))
    # np.roots takes the coefficients from the highest power down
    return np.roots([(-1.)**j/factorial(j) for j in range(order,-1,-1)])


def anneal(sc,h,wf,T,dbeta=None,step=2.0,order=8,erange=None):
    """Return e^{-H/(2T)}|wf>, normalized, by imaginary-time steps.

    Purification convention: |Psi(beta)> = e^{-beta*H/2}|Psi(0)>, so that
    tracing out the ancilla out of |Psi(beta)><Psi(beta)| gives
    rho ~ e^{-beta*H}, i.e. the physical thermal state at the requested
    T = 1/beta. Applying the full beta here (as this used to) doses the
    wavefunction with e^{-beta*H} instead of e^{-beta*H/2}, which after
    tracing out the ancilla is the thermal state of T/2, not T -- confirmed
    numerically against exact ED thermal averages (see docs/user_guide).

    Each of the nst steps applies p_k(dtau*(H-E_ref)), the order-k Taylor
    polynomial of exp(-dtau*(H-E_ref)), with E_ref the middle of the
    spectrum [E_min, E_max] (erange, computed with band_edges() when not
    given) and dtau = beta_half/nst, nst = ceil(beta_half*W/step) and at
    least 1, W = E_max-E_min. So x = dtau*(E-E_ref) lies in
    [-step/2, step/2] for every eigenvalue E, in any units and at any
    constant offset, and the local inverse temperature of the purified
    weights p_k(x)^(2 nst) is off by a relative (step/2)^k/k! at worst,
    at the band edges, and by far less where the thermal weight sits:
    the error is set by the dimensionless step and the order alone, and it
    converges as step^k. The polynomial is applied as its k linear factors
    (1 - dtau*(H-E_ref)/r_j) over its complex roots, i.e. k applications
    of a MultiOperator per step and nothing else, so this runs on every
    backend. dbeta, if given, is an additional cap on dtau in the
    Hamiltonian's own units (what dbeta alone used to be).

    What this replaced, and why each part had to go (2026-09-25b audit,
    findings 17 to 19): nst = int(beta_half/0.1) first-order steps
    (1 - dtau*H) on the unshifted H. Its error was set by 0.1 times the
    absolute, extensive energy (12.3 per cent off at n=10 and T=0.5, and
    sign-inverting factors at a +20 offset); int() gave zero steps above
    T=5, a ZeroDivisionError for a float T and the T=infinity state for a
    numpy one; and it returned early once one step changed the state by
    less than 1e-7, a test in absolute units that said nothing about the
    imaginary time still to go (the T=infinity state for any Hamiltonian
    written below about 7e-3, and a manifold split by 1e-2 at unit scale
    left unresolved). There is no early return now: the step count is
    finite at every T>0, and on an exact eigenstate the remaining steps
    are no-ops."""
    T = float(T)
    if not T>0.: raise ValueError("anneal: T must be > 0, got "+repr(T))
    beta_half = 1./(2.*T) # half-beta part
    wf = wf.normalize() # normalize
    if beta_half==0.: return wf # T=infinity: the input is the purification
    if erange is None: erange = band_edges(sc,h)
    emin, emax = float(erange[0]), float(erange[1])
    width = emax-emin
    if not width>0.: return wf # h is a multiple of the identity
    dtau = step/width # largest step allowed by the dimensionless step
    if dbeta is not None: dtau = min(dtau,float(dbeta))
    x = beta_half/dtau
    if not np.isfinite(x):
        raise ValueError("anneal: T="+repr(T)+" is too small for the "
                         "spectral width "+repr(width))
    # ceil, not int: int() rounded a count like 2.9999999999999996 (0.3/0.1)
    # down; the slack keeps a count that is an integer up to rounding at
    # that integer, so that a rescaled Hamiltonian takes the same steps
    nst = max(1,int(np.ceil(x-1e-9)))
    dtau = beta_half/nst # the step actually taken
    eref = 0.5*(emin+emax) # middle of the spectrum, see above
    coefs = [dtau/r for r in taylor_exp_roots(order)]
    # (1 - a*(H-eref)) as one operator per root. On the ED backend each is
    # assembled into a matrix once, rather than on every application
    from .edtk.edchain import State
    if isinstance(wf,State):
        import scipy.sparse as sp
        hm = sp.csr_matrix(wf.MBO.MO2matrix(h))
        one = sp.identity(hm.shape[0],dtype=complex,format="csr")
        factors = [(1.+a*eref)*one - a*hm for a in coefs]
    else: factors = [(1.+a*eref) - a*h for a in coefs]
    for i in range(nst): # do as many steps
        wf1 = factors[0]*wf # first factor of this step
        # <wf|1-a(H-eref)|wf> = 1-a(<H>-eref), so the energy printed at
        # every step costs no extra application
        e = eref + (1.-wf.dot(wf1))/coefs[0]
        print("Annealing, energy",np.real(e))
        for f in factors[1:]: wf1 = f*wf1 # the rest of the polynomial
        wf = wf1.normalize() # normalize
    return wf









# The keywords anneal_options takes, and the defaults of the purification.
# The shift is not one of them: it is the lower band edge, which
# band_edges() has computed already.
_IMPLICIT_KEYWORDS = ("nsteps","tau0","order","nsweeps","tol","krylov",
                      "maxiter")
_IMPLICIT_DEFAULTS = dict(order=8,nsteps=2)
_AUTO_RATIO = 4. # Taylor applications per implicit solve to take the grid
_AUTO_TMAX = 0.05 # and the highest T, in units of the largest coupling


def anneal_counts(tc,T,erange):
    """(applications of anneal(), solves of anneal_implicit()) to reach
    T, at the stepper settings of tc and the defaults of the grid"""
    from . import imaginarytime
    beta_half = 1./(2.*float(T))
    width = float(erange[1])-float(erange[0])
    nst = max(1,int(np.ceil(beta_half*width/tc.anneal_step-1e-9)))
    u = imaginarytime._unit(tc.hamiltonian)
    tau0 = 0.01/u if u>0. else 0.01
    opts = _IMPLICIT_DEFAULTS
    steps = imaginarytime.log_grid(beta_half,tau0,opts["nsteps"])
    return nst*int(tc.anneal_order), len(steps)*(opts["order"]//2)


def auto_anneal_method(tc,wf,T,erange):
    """What anneal_method=None takes: "implicit" on itensor_version=3, when
    the implicit steps can run there, T is at most _AUTO_TMAX times the
    largest coupling of H and anneal() needs at least _AUTO_RATIO times as
    many MPO applications as the grid needs solves, and "taylor"
    otherwise, as on every other backend and on mode="ED".

    The temperature bound is the range the grid was measured in: much
    faster at T=0.05 and 0.02, and as accurate as anneal() from T=0.1 down
    (4.9e-6 against 1.5e-6 at T=0.5). The count alone would take the grid
    at high T on a long chain (about 280 applications against 64 solves
    at T=0.5 on 100 sites), where it is less accurate and was not timed.

    The ratio is where the grid was measured much faster (anneal_implicit's
    docstring): one complex solve costs about 3 applications on v3, so a
    count ratio of 4.9 (16 sites, T=0.05) gave 1.87 times the speed and
    10.3 (16 sites, T=0.02) 2.80 times, while 3.0 (10 sites, T=0.05) gave
    only 1.18 and the grid is also less accurate at high T. On "python" a
    solve costs more: a ratio of 4.9 gave 0.96. A chain whose anneal_step
    or anneal_order was set keeps the Taylor stepper it configured."""
    from . import imaginarytime
    from .mpsalgebra import wavefunction_mode
    sc = tc.MBChain
    if "anneal_step" in vars(tc) or "anneal_order" in vars(tc): return "taylor"
    if sc.itensor_version!=3 or wavefunction_mode(wf)!="DMRG": return "taylor"
    if erange is None or not 0.<float(T)<np.inf: return "taylor"
    # simplified, since set_hamiltonian stores (h+h^dagger)/2 term by term,
    # whose largest coefficient is half the coupling
    u = imaginarytime._unit(tc.hamiltonian.simplify())
    if not float(T)<=_AUTO_TMAX*u: return "taylor"
    if not float(erange[1])-float(erange[0])>0.: return "taylor"
    if not imaginarytime._implicit_available(sc): return "taylor"
    taylor, implicit = anneal_counts(tc,T,erange)
    return "implicit" if taylor>=_AUTO_RATIO*implicit else "taylor"


def anneal_implicit(sc,h,wf,T,erange=None,options=None):
    """Return e^{-H/(2T)}|wf>, normalized, by implicit steps on a
    logarithmic grid in imaginary time (imaginarytime.evolve,
    arXiv:2606.02930), from 0 to beta/2 with the shift at the lower band
    edge E_min (erange[0], computed with band_edges() when not given).

    The cost is a number of linear solves that grows as log(beta), against
    about 2*beta*W MPO applications for anneal(), W the spectral width,
    which is extensive. The default is the 4-stage Gauss-Legendre rule
    (order=8) with 2 steps per panel. The reason for this is that the
    purification starts from the T=infinity state, which holds every
    energy of the spectrum: the trapezoid rule with 8 steps per panel, the
    default of METTS, is off by 1e-3 (relative) in the thermal energy of a
    10-site Heisenberg chain at T=0.5, and order=8 with nsteps=2 by 7e-9,
    against 1.2e-7 for anneal() at its defaults (closed forms over the
    exact spectrum, tests/test_thermal_anneal_implicit.py).

    On MPS the floor is the truncation, not the time step: every factor is
    a sweep that truncates at every bond, four of them per step, so at
    cutoff=1e-12 the energy of a 4-site chain is 4.1e-6 off at T=1 and
    6e-8 at T=0.1, on both backends, against 3.1e-8 at T=1 for anneal().
    Measured against anneal() on itensor_version=3 at maxm=64 (one thread,
    one P-core), 10 sites: 0.45, 0.91, 1.18 and 1.88 times its speed at
    T=0.5, 0.1, 0.05 and 0.02, 4.9e-6 off the Boltzmann energy at T=0.5
    (anneal() 1.5e-6) and within 9e-7 below; 16 sites: 1.87 and 2.80 times
    at T=0.05 and 0.02, the two energies 2e-6 apart. One complex
    solve costs about 3.4 MPO applications there. On "python" a solve
    costs more against an application: 0.96 times the speed of anneal()
    on 8 sites at T=0.02. A longer first panel,
    tau0 = 2/W, halves the solves in closed form at the same accuracy but
    not the time on MPS: the GMRES of a long step takes more iterations,
    and one sweep cannot follow a long step out of the T=infinity state
    (8e-5 off at T=0.5), so tau0 stays at the default of
    evolve_imaginary_time.

    sc holds the singlet Hamiltonian on entry; h is set on it here, since
    the steps solve with the chain's own Hamiltonian."""
    from . import imaginarytime
    opts = dict(_IMPLICIT_DEFAULTS)
    opts.update(dict(options or {}))
    unknown = sorted(set(opts)-set(_IMPLICIT_KEYWORDS))
    if unknown:
        raise TypeError("Thermal_Spin_Chain.anneal_options: unknown "
                        "keyword(s) %s; the known ones are %s"
                        %(unknown,list(_IMPLICIT_KEYWORDS)))
    T = float(T)
    if not T>0.: raise ValueError("anneal_implicit: T must be > 0, got "+repr(T))
    beta_half = 1./(2.*T)
    wf = wf.normalize()
    if beta_half==0.: return wf # T=infinity: the input is the purification
    if erange is None: erange = band_edges(sc,h)
    emin, emax = float(erange[0]), float(erange[1])
    if not emax-emin>0.: return wf # h is a multiple of the identity
    sc.set_hamiltonian(h)
    wf = imaginarytime.evolve(sc,wf,beta_half,shift=emin,**opts)[2]
    return wf.normalize()
