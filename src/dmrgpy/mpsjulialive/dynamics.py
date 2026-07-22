import numpy as np

from .. import multioperator
from ..algebra import kpm
from ..algebra.kpm import generate_profile


def _max_energy_bound(self,H):
    """Variational upper bound on H's spectrum, from a reduced-effort
    ground-state search on -H (mirrors pyitensor.chain.Chain's own
    _maximum_energy: a Chebyshev-window bound, never a physical result).
    Runs in place on `self` and restores its Hamiltonian/ground-state
    cache afterward -- julia_live has no deepcopy-safe clone() for its
    live Julia session handles (self.jlsites/self.wf0.jlmps aren't part
    of Python's copy protocol, see manybodychain.py's __deepcopy__), so
    this can't go through bandwidth()/lowest_eigenvalue() the way the
    C++/pure-Python backends do."""
    maxm0,nsweeps0 = self.maxm,self.nsweeps
    self.maxm = min(self.maxm,20)
    self.nsweeps = min(self.nsweeps,5)
    self.restart()
    self.set_hamiltonian(-1*H,restart=False)
    emax = -self.gs_energy()
    self.maxm,self.nsweeps = maxm0,nsweeps0
    self.restart()
    self.set_hamiltonian(H,restart=False)
    return emax


def _same_mps(jlsites,vi,vj,maxm):
    from .juliasession import Main as Mainjl
    return bool(Mainjl.same_mps(vi,vj,maxm))


def _kpm_moments_full(jlsites,jlmpo,vi,vj,n,kpmmaxm,kpmcutoff):
    """Chebyshev moments <vj|T_k(scaledH)|vi>, via the standard three-term
    recursion T_{k+1} = 2*scaledH*T_k - T_{k-1}. The whole loop runs
    natively in Julia (kpm.jl's kpm_moments_full, driven through
    mpsalgebra.jl's applyoperator/summps) rather than one Python<->Julia
    round trip per Chebyshev step -- confirmed directly that the
    per-step-round-trip version this replaced was ~1.7x slower than the
    compiled ITensor v3 backend even with a warm Julia session (30-site
    Heisenberg chain), which this closes."""
    from .juliasession import Main as Mainjl
    out = Mainjl.kpm_moments_full(jlsites,jlmpo,vi,vj,n,kpmmaxm,kpmcutoff)
    return list(out)


def _kpm_moments_accelerated(jlsites,jlmpo,vi,n,kpmmaxm,kpmcutoff):
    """Same recursion as _kpm_moments_full, specialized for vi==vj (a
    diagonal correlator): the even/odd moment pair at each step is read
    off two consecutive Chebyshev vectors instead of two separate
    recursions, roughly halving the MPO applications. Native Julia loop,
    see _kpm_moments_full's note."""
    from .juliasession import Main as Mainjl
    out = Mainjl.kpm_moments_accelerated(jlsites,jlmpo,vi,n,kpmmaxm,kpmcutoff)
    return list(out)


def get_dynamical_correlator(self,submode="KPM",**kwargs):
    """Dispatch a dynamical correlator computation on the Julia backend.
    submode="KPM" (see docs/documentation.md §4.6) and submode="CVM" are
    implemented; other submodes (TD -- use
    timedependent.dynamical_correlator/evolution_DC instead, which
    already dispatch to julia_live -- TDZ, EX, maxent, ...) are only
    implemented for the C++/pure-Python backends (dynamics.py).

    Deliberately takes submode + **kwargs only (mirroring the top-level
    dynamics.py::get_dynamical_correlator's own signature), not named
    parameters like name=/es=/delta=: an earlier version of this function
    declared those as its own named parameters with KPM-flavored
    defaults, which silently captured the caller's es=/name=/delta=
    kwargs out of **kwargs before they could reach cvm.dynamical_correlator
    (whose own, different defaults would then be used unless the caller's
    values happened to be explicitly re-forwarded) -- confirmed directly,
    this caused CVM to silently run on a totally different frequency grid
    than the one requested. Each submode's own function applies its own
    defaults now, exactly as the C++/pure-Python dispatch already does."""
    if submode=="CVM":
        # cvm.py's CG solve is already backend-agnostic MPS/MPO algebra
        # (self.toMPO()/MPS +,-,scalar-*,.dot()), which already works for
        # julia_live -- the only thing that needed fixing was its own
        # self._session.set_sweep_params(...) calls, see
        # cvm.py::_set_cvm_sweep_params.
        from .. import cvm
        return cvm.dynamical_correlator(self,**kwargs)
    if submode!="KPM":
        raise NotImplementedError(
            "itensor_version='julia_live' only implements submode='KPM'/"
            "'CVM' for get_dynamical_correlator, got submode=%r"%submode)
    return _kpm_dynamical_correlator(self,**kwargs)


def _kpm_dynamical_correlator(self,n=1000,
             name=None,delta=1e-1,kernel="jackson",
             es=np.linspace(-1.,10,500),deconvolve=None,
             **kwargs):
    """The submode="KPM" implementation, factored out of
    get_dynamical_correlator so its own named-parameter defaults can't
    shadow another submode's kwargs (see that function's docstring).
    Mirrors kpmdmrg.py::get_dynamical_correlator (the C++/pure-Python KPM
    entry point) and its Chebyshev-moment post-processing exactly, but
    computes the moments with the mpsjulialive Julia MPS/MPO primitives
    instead of self._session.kpm_dynamical_correlator (julia_live has no
    such session object)."""
    if delta<0.0: raise
    if self.kpm_extrapolate: delta = delta*self.kpm_extrapolate_factor
    if type(name[0])!=multioperator.MultiOperator: raise
    mi = name[1] # first operator
    mj = name[0].get_dagger() # second operator
    H = self.hamiltonian
    e0 = self.gs_energy() # compute ground state (also sets self.e0/self.wf0)
    wf0 = self.wf0
    emin = self.e0
    emax = _max_energy_bound(self,H)
    self.wf0,self.e0,self.computed_gs = wf0,e0,True # restore GS cache
    shift = -(emin+emax)/2.0
    scale = 1.0/((emax-emin)*self.kpm_scale)
    n = int(round((emax-emin)/delta))*self.kpm_n_scale
    Hscaled_MO = (H+shift*multioperator.identity())*scale
    from .mpo import MPO
    from .juliasession import Main as Mainjl
    Hscaled = MPO(Hscaled_MO,MBO=self)
    Aop = MPO(mi,MBO=self)
    Bop = MPO(mj,MBO=self)
    psi1 = Mainjl.applyoperator(self.jlsites,Aop.jlmpo,wf0.jlmps,
            self.kpmmaxm,self.kpmcutoff)
    psi2 = Mainjl.applyoperator(self.jlsites,Bop.jlmpo,wf0.jlmps,
            self.kpmmaxm,self.kpmcutoff)
    if self.kpm_accelerate and _same_mps(self.jlsites,psi1,psi2,self.kpmmaxm):
        moments = _kpm_moments_accelerated(self.jlsites,Hscaled.jlmpo,psi1,
                n,self.kpmmaxm,self.kpmcutoff)
    else:
        moments = _kpm_moments_full(self.jlsites,Hscaled.jlmpo,psi1,psi2,
                n,self.kpmmaxm,self.kpmcutoff)
    mus = np.array(moments)
    if self.kpm_extrapolate:
        mus = kpm.extrapolate_moments(mus,fac=self.kpm_extrapolate_factor,
                extrapolation_mode=self.kpm_extrapolate_mode)
    xs = 0.99*np.linspace(-1.0,1.0,int(n*10),endpoint=False) # energies
    ys = generate_profile(mus,xs,use_fortran=False,kernel=kernel) # generate the DOS
    xs /= scale # scale back the energies
    xs += (emin+emax)/2. -emin # shift the energies
    ys *= scale # renormalize the y values
    from scipy.interpolate import interp1d
    fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
    fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
    return (es,fr(es)+1j*fi(es)) # interpolate
