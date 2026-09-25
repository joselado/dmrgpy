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
            erange = None
            if 0.<T<np.inf: erange = band_edges(self.MBChain,self.hamiltonian)
            if T>0. and (erange is None or T>1e-5*(erange[1]-erange[0])):
                self.MBChain.set_hamiltonian(h) # singlet Hamiltonian
                wf = self.MBChain.get_gs() # get the fully entangled WF
                wf0 = anneal(self.MBChain,self.hamiltonian,wf,T,
                             step=self.anneal_step,order=self.anneal_order,
                             erange=erange)
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







