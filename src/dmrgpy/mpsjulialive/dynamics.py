import numpy as np

from .. import multioperator
from ..algebra import kpm
from ..kpmdmrg import dynamical_correlator_from_moments


def _max_energy_bound(self,H):
    """Variational upper bound on H's spectrum, from a reduced-effort
    ground-state search on -H (mirrors pyitensor.chain.Chain's own
    _maximum_energy: a Chebyshev-window bound, never a physical result).

    Runs on a clone, the bandwidth()/lowest_eigenvalue() path of the other
    backends, so the chain's own state, Hamiltonian and solve marks are
    never touched: a Julia deepcopy keeps the ITensor Index ids, so the
    clone's sites are the chain's. It used to run in place on the chain
    and put everything back afterwards; the solve is the same one (a cold
    start at min(maxm,20) and min(nsweeps,5)), so the bound is unchanged."""
    mbc = self.clone()
    mbc.maxm = min(self.maxm,20)
    mbc.nsweeps = min(self.nsweeps,5)
    mbc.set_hamiltonian(-1*H) # restart=True: a cold start
    return -mbc.gs_energy()


def _min_energy(self,H):
    """H's lower band edge from a solve at the chain's own parameters, for
    the Chebyshev window when the chain's e0 is not marked as that edge
    (groundstate.mark_lower_edge()): a state set with set_gs() or
    set_initial_wf(), gs_energy_generalized()'s, or one swept from a warm
    start (wf0=, reconverge=True, restart=False). On a clone, like
    _max_energy_bound() above, so the chain's state stays as it is."""
    mbc = self.clone()
    mbc.set_hamiltonian(H) # restart=True: a cold start
    return mbc.gs_energy()


def _same_mps(vi,vj,maxm,cutoff):
    from .juliasession import Main as Mainjl
    return bool(Mainjl.same_mps(vi,vj,maxm,cutoff))


def _kpm_moments_full(jlmpo,vi,vj,n,kpmmaxm,kpmcutoff):
    """Chebyshev moments <vj|T_k(scaledH)|vi>, via the standard three-term
    recursion T_{k+1} = 2*scaledH*T_k - T_{k-1}. The whole loop runs
    natively in Julia (kpm.jl's kpm_moments_full, driven through
    mpsalgebra.jl's own apply_op()/summps) rather than one
    Python<->Julia round trip per Chebyshev step -- confirmed directly
    that the per-step-round-trip version this replaced was ~1.7x slower
    than the compiled ITensor v3 backend even with a warm Julia session
    (30-site Heisenberg chain), which this closes. apply_op() (rather
    than mpsalgebra.jl's applyoperator(), which does an extra, redundant
    MPO contraction per call -- see its own docstring) roughly halves
    the cost again on top of that, since this recursion is O(n moments)
    calls and dominates KPM's total runtime."""
    from .juliasession import Main as Mainjl
    out = Mainjl.kpm_moments_full(jlmpo,vi,vj,n,kpmmaxm,kpmcutoff)
    return list(out)


def _kpm_moments_accelerated(jlmpo,vi,n,kpmmaxm,kpmcutoff):
    """Same recursion as _kpm_moments_full, specialized for vi==vj (a
    diagonal correlator): the even/odd moment pair at each step is read
    off two consecutive Chebyshev vectors instead of two separate
    recursions, roughly halving the MPO applications. Native Julia loop,
    see _kpm_moments_full's note."""
    from .juliasession import Main as Mainjl
    out = Mainjl.kpm_moments_accelerated(jlmpo,vi,n,kpmmaxm,kpmcutoff)
    return list(out)


def kpm_moments(self,mi,mj,delta,n_scale):
    """The KPM moments <mj|T_k(H_scaled)|mi> of the chain's current state
    on julia_live, what the session backends' kpm_dynamical_correlator
    returns: (moments,emin,emax,scale,n), with the moments not yet cut to
    the calibrated count n (kpmdmrg.dynamical_correlator_moments does the
    cut and kpm_extrapolate for every backend). mi=B and mj=A^dagger, as
    kpmdmrg passes them; delta is the effective resolution. The ground
    state is read, not solved, so the caller has it current already.

    The moments are computed with kpm.jl's native recursions
    (kpm_moments_full, or kpm_moments_accelerated when B|GS> and
    A^dagger|GS> are the same vector and kpm_accelerate is on)."""
    H = self.hamiltonian
    e0 = self.e0
    wf0 = self.wf0
    # The state is measured from its own energy e0, and the window takes
    # its lower edge from e0 only when a plain solve marked e0 as H's lower
    # edge (groundstate.mark_lower_edge), and from a solve of its own
    # otherwise, as the session backends do (2026-09-24c audit, finding 1).
    from ..groundstate import e0_is_lower_edge
    emin = e0 if e0_is_lower_edge(self) else _min_energy(self,H)
    emax = _max_energy_bound(self,H)
    shift = -(emin+emax)/2.0
    scale = 1.0/((emax-emin)*self.kpm_scale)
    # the shared calibrated moment count, see
    # algebra/kpm.py::polynomials_for_broadening and kpmdmrg.py
    from ..algebra.kpm import polynomials_for_broadening
    n = polynomials_for_broadening(1.0/scale,delta,n_scale=n_scale)
    Hscaled_MO = (H+shift*multioperator.identity())*scale
    from .mpo import MPO
    from .juliasession import Main as Mainjl
    Hscaled = MPO(Hscaled_MO,MBO=self)
    Aop = MPO(mi,MBO=self)
    Bop = MPO(mj,MBO=self)
    psi1 = Mainjl.apply_op(Aop.jlmpo,wf0.jlmps,self.kpmmaxm,self.kpmcutoff)
    psi2 = Mainjl.apply_op(Bop.jlmpo,wf0.jlmps,self.kpmmaxm,self.kpmcutoff)
    if self.kpm_accelerate and _same_mps(psi1,psi2,self.kpmmaxm,self.kpmcutoff):
        moments = _kpm_moments_accelerated(Hscaled.jlmpo,psi1,
                n,self.kpmmaxm,self.kpmcutoff)
    else:
        moments = _kpm_moments_full(Hscaled.jlmpo,psi1,psi2,
                n,self.kpmmaxm,self.kpmcutoff)
    return np.array(moments),emin,emax,scale,n


def general_kpm_moments(self,X,wfa,wfb,num_p,accelerate):
    """Chebyshev moments <wfb|T_k(X)|wfa> of an already rescaled operator
    X, the julia_live counterpart of the sessions' general_kpm: the same
    kpm.jl recursions, at the chain's maxm and cutoff, the accelerated
    one when accelerate is on and the two states are the same vector."""
    from .mpo import MPO
    from .juliasession import Main as Mainjl
    Xop = MPO(X,MBO=self)
    vi,vj = wfa.jlmps,wfb.jlmps
    if accelerate and _same_mps(vi,vj,self.maxm,self.cutoff):
        out = _kpm_moments_accelerated(Xop.jlmpo,vi,int(num_p),
                self.maxm,self.cutoff)
    else:
        out = _kpm_moments_full(Xop.jlmpo,vi,vj,int(num_p),
                self.maxm,self.cutoff)
    return np.array(out)


def kpm_dynamical_correlator(self,name=None,delta=1e-1,kernel="jackson",
        es=np.linspace(-1.,10,500),hodc_order=6,hodc_eta=None,**kwargs):
    """submode="KPM" on julia_live: kpmdmrg.dynamical_correlator_moments,
    whose julia_live branch calls kpm_moments() above, then the shared
    reconstruction. **kwargs goes to dynamical_correlator_moments, which
    takes i=/j= and raises TypeError on anything else. Kept here, rather
    than in kpmdmrg.get_dynamical_correlator, so that the reconstruction
    stays this module's dynamical_correlator_from_moments."""
    from ..kpmdmrg import dynamical_correlator_moments
    mus,emin,emax,scale,n,delta = dynamical_correlator_moments(self,
            name=name,delta=delta,**kwargs)
    # the frequency origin is the measured state's own energy, as on the
    # session backends (2026-09-24c audit, finding 1)
    return dynamical_correlator_from_moments(mus,emin,emax,scale,n,es,
            kernel=kernel,delta=delta,
            hodc_order=hodc_order,hodc_eta=hodc_eta,origin=self.e0)
