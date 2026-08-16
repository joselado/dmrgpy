from __future__ import print_function
from . import operatornames
import numpy as np
from scipy.interpolate import interp1d
from . import multioperator
from .edtk import timedependent as tded



# The exact substring both bond_hamiltonians() implementations raise on a
# term spanning 3+ sites (pyitensor/tebd.py's NotImplementedError,
# mpscpp3/tebd.h's ITError -- surfaced to Python as a RuntimeError since
# ITError derives from std::runtime_error and pybind11 translates that
# automatically). tevol_method="AUTO" matches on this specific text so it
# only swallows the "not nearest-neighbor" condition it exists to handle,
# not an unrelated bug that happens to also raise NotImplementedError/
# RuntimeError from somewhere inside the TEBD call.
_TEBD_NN_ERROR_MARKER = "nearest-neighbor"


def _is_tebd_nn_error(exc):
    """True if `exc` is TEBD's own rejection of a non-nearest-neighbor
    Hamiltonian, as opposed to some other failure raised while attempting
    it."""
    return _TEBD_NN_ERROR_MARKER in str(exc)


def _tebd_or_tdvp(tebd_call,tdvp_call):
    """Backs tevol_method="AUTO": try `tebd_call` (cheaper -- no per-step
    Krylov/Lanczos, see tevol_method's docstring in manybodychain.py) and
    transparently fall back to `tdvp_call` if the Hamiltonian turns out
    not to be strictly nearest-neighbor. Both backends already do this
    check once, up front, before touching the wavefunction (see
    bond_hamiltonians() in pyitensor/tebd.py / mpscpp3/tebd.h), so retrying
    with TDVP here costs at most one discarded MPO build, not a discarded
    partial time evolution."""
    try:
        return tebd_call()
    except (NotImplementedError,RuntimeError) as exc:
        if not _is_tebd_nn_error(exc): raise
        return tdvp_call()


def evolution_DC(self,mode="DMRG",**kwargs):
    if mode=="DMRG":  return evolution_dmrg_DC(self,**kwargs)
    if mode=="ED": 
        edobj = self.get_ED_obj() # get the ED object
        return tded.evolution_DC(edobj,h=self.hamiltonian,**kwargs)



def evolution_dmrg_DC(self,name="XX",nt=10000,dt=0.1,restart=True,**kwargs):
    """
    Real-time quench dynamical correlator via the in-process pybind11
    extension.

    **Measure before evolving.** Every backend's evolution loop
    (Chain::quench*/evolve_and_measure* in mpscpp2/mpscpp3, the same
    methods in pyitensor/chain.py, mpsjulialive/{tdvp,tebd}.jl, and
    edtk/timedependent.py) records its observable at the *top* of the
    step loop, so `correlator[k]` is the value at `t = k*dt` and lines up
    with the `ts = [0, dt, ..., (nt-1)*dt]` grid returned just below.

    This used to be the other way round -- evolve first, measure after --
    which made `correlator[k]` the value at `(k+1)*dt` while still
    labelling it `k*dt`. That was not cosmetic. It put a spurious
    `exp(i*omega*dt)` phase on the Fourier transform, and, worse, it
    dropped `C(0)` -- the single largest sample -- from the Riemann sum in
    `_fourier_transform_correlator` entirely, leaving an O(dt*C(0)) error
    in every submode="TD" spectrum that did *not* vanish as `nt` grew at
    fixed total time. Confirmed directly on an L=10 Heisenberg chain:
    `correlator[0]` came back as 0.2409/0.2477/0.2494 for dt=0.2/0.1/0.05
    against an exact C(0)=<A B>=0.25 (the imaginary part halving exactly
    with dt), and the resulting spectral weight moved 75% across
    dt=0.1/0.05/0.025 (-0.105/-0.159/-0.184) for a C(t) that is itself
    dt-independent to 5e-5. With the measurement moved to the top of the
    loop the same sweep gives -0.2047/-0.2061/-0.2065, i.e. converging
    rather than drifting, and `correlator[0]` reproduces C(0) to 1e-13 on
    all of itensor_version 2, 3, "python" and mode="ED".

    Because the ED path uses the identical convention, every
    DMRG-vs-ED cross-check in tests/ and examples/ compares like with
    like; the fix changes both sides together.

    Defaults to TDVP (mpscpp3/chain_session.h's Chain::quench_tdvp(), see
    TDVP/ and self.tevol_method) for itensor_version=3 or "python" (the
    pure-Python backend has its own TDVP, see pyitensor/tdvp.py); falls
    back to the legacy MPO-Taylor Chain::quench() otherwise
    (itensor_version=2, or self.tevol_method="MPO" explicitly).
    self.tevol_method="TDVP_GSE" instead runs one-site TDVP with Krylov
    global subspace expansion (Chain::quench_tdvp_gse(), arXiv:2005.06104)
    for the first self.tdvp_gse_sweeps steps -- same itensor_version
    support as "TDVP" (3 or "python" only). A v2-API port
    (mpscpp2/TDVP/) was attempted and briefly landed here but was removed:
    it was numerically correct (verified against ED and against v3/
    "python") but had a severe, unresolved performance regression at
    n>~10 sites (the dynamical-correlator step didn't finish in 25
    minutes at n=12, versus under a second for the same computation on
    v3/"python") that couldn't be root-caused in the time available.

    self.tevol_method="TEBD" instead runs 2nd-order-Trotter TEBD (gates
    built once from the bare nearest-neighbor bond Hamiltonians, reused
    unchanged every step, no per-step Krylov/Lanczos at all) --
    itensor_version=3 or "python" only (itensor_version="python" via
    pyitensor/tebd.py's TEBDEvolver, itensor_version=3 via
    mpscpp3/tebd.h's bond_hamiltonians()/build_tebd_gates()/tebd_step(),
    a C++ port onto ITensor's own BondGate primitive), and only for a
    strictly nearest-neighbor self.hamiltonian (both backends raise --
    NotImplementedError in Python, ITError in C++ -- for any term
    spanning 3+ distinct sites; fall back to "TDVP" for longer-range
    models). self.tevol_method="AUTO" makes that fallback automatic
    (see _tebd_or_tdvp() above): try "TEBD" first, and transparently
    retry as "TDVP" if it turns out self.hamiltonian isn't nearest-
    neighbor -- "TEBD" itself stays a hard opt-in (it still raises) so a
    caller who explicitly asked for it is told when its assumption
    doesn't hold, rather than silently getting a different integrator.

    fit_td is hardcoded False in the MPO fallback, not read from
    self.fit_td: the removed file-based backend wrote it to tasks.in under
    the key "tevol_fit", but time_evolution.h actually read
    "tevol_fit_td" (a pre-existing, unrelated key-name mismatch) -- so the
    fitApplyMPO branch there was unreachable regardless of self.fit_td,
    and False reproduces that actual behavior rather than the
    intended-but-never-taken one.

    "restart" has no effect: quench()'s C++ implementation always starts
    from get_gs() regardless of its value.
    """
    if self.itensor_version=="julia_live":
        from .mpsjulialive import timedependent as tdjl
        return tdjl.evolution_dmrg_DC(self,name=name,nt=nt,dt=dt,**kwargs)
    name = operatornames.str2MO(self,name,**kwargs)
    name[0] = name[0].get_dagger()
    A,B = name[0],name[1]
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    if self.itensor_version in (3,"python") and self.tevol_method=="TDVP":
        correlator,_wf = self._session.quench_tdvp(
                self.hamiltonian.to_terms(),A.to_terms(),B.to_terms(),
                int(nt),dt)
    elif self.itensor_version in (3,"python") and self.tevol_method=="TEBD":
        correlator,_wf = self._session.quench_tebd(
                self.hamiltonian.to_terms(),A.to_terms(),B.to_terms(),
                int(nt),dt)
    elif self.itensor_version in (3,"python") and self.tevol_method=="AUTO":
        correlator,_wf = _tebd_or_tdvp(
                lambda: self._session.quench_tebd(
                    self.hamiltonian.to_terms(),A.to_terms(),B.to_terms(),
                    int(nt),dt),
                lambda: self._session.quench_tdvp(
                    self.hamiltonian.to_terms(),A.to_terms(),B.to_terms(),
                    int(nt),dt))
    elif self.itensor_version in (3,"python") and self.tevol_method=="TDVP_GSE":
        correlator,_wf = self._session.quench_tdvp_gse(
                self.hamiltonian.to_terms(),A.to_terms(),B.to_terms(),
                int(nt),dt,self.tdvp_gse_sweeps,self.tdvp_gse_krylov_order,
                self.tdvp_gse_cutoff)
    else:
        correlator,_wf = self._session.quench(
                self.hamiltonian.to_terms(),A.to_terms(),B.to_terms(),
                int(nt),dt,False)
    cs = np.array(correlator)
    ts = np.array([dt*ii for ii in range(nt)])
    return ts,cs.real-1j*cs.imag



def evolve_and_measure(self,mode="DMRG",**kwargs):
    """Evolve and measure"""
    if mode=="DMRG": return evolve_and_measure_dmrg(self,**kwargs)
    elif mode=="ED": 
        edobj = self.get_ED_obj() # get the ED object
        h = self.hamiltonian # get the ED object
        return tded.evolve_and_measure(edobj,h,**kwargs)



def evolve_and_measure_dmrg(self,operator=None,nt=1000,h=None,
        dt=1e-2,wf=None,return_wf=False,**kwargs):
    """
    Real-time evolution + measurement via the in-process pybind11
    extension.

    Defaults to TDVP (mpscpp3/chain_session.h's
    Chain::evolve_and_measure_tdvp(), see TDVP/ and self.tevol_method) for
    itensor_version=3 or "python"; falls back to the legacy MPO-Taylor
    Chain::evolve_and_measure() otherwise (itensor_version=2, or
    self.tevol_method="MPO" explicitly). self.tevol_method="TDVP_GSE"
    instead runs one-site TDVP with Krylov global subspace expansion
    (Chain::evolve_and_measure_tdvp_gse(), arXiv:2005.06104) for the first
    self.tdvp_gse_sweeps steps -- see evolution_dmrg_DC's docstring.
    self.tevol_method="TEBD"/"AUTO" behave exactly as documented there too.

    fit_td is hardcoded False in the MPO fallback, for the same reason as
    evolution_dmrg_DC (see its docstring): the "tevol_fit"/"tevol_fit_td"
    key-name mismatch meant the old file-based backend's fitApplyMPO
    branch was unreachable regardless of self.fit_td.

    return_wf=True additionally returns the final wavefunction (wrapped as
    an mps.MPS, see mpsalgebra.py's exponential_dmrg() for the same
    cpp_handle-wrapping pattern) as a third element -- e.g. to chain a
    forward evolution into a subsequent backward one for a round-trip
    fidelity check where ED isn't feasible (see
    examples/tdvp_VS_ED_time_evolution/benchmark_scaling.py).
    """
    if self.itensor_version=="julia_live":
        from .mpsjulialive import timedependent as tdjl
        return tdjl.evolve_and_measure_dmrg(self,operator=operator,nt=nt,
                h=h,dt=dt,wf=wf,return_wf=return_wf,**kwargs)
    if h is None: h = self.hamiltonian # Hamiltonian
    if wf is None: wf = self.wf0 # get ground state
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    if self.itensor_version in (3,"python") and self.tevol_method=="TDVP":
        correlator,_wf = self._session.evolve_and_measure_tdvp(
                h.to_terms(),operator.to_terms(),wf.cpp_handle,
                int(nt),dt)
    elif self.itensor_version in (3,"python") and self.tevol_method=="TEBD":
        correlator,_wf = self._session.evolve_and_measure_tebd(
                h.to_terms(),operator.to_terms(),wf.cpp_handle,
                int(nt),dt)
    elif self.itensor_version in (3,"python") and self.tevol_method=="AUTO":
        correlator,_wf = _tebd_or_tdvp(
                lambda: self._session.evolve_and_measure_tebd(
                    h.to_terms(),operator.to_terms(),wf.cpp_handle,
                    int(nt),dt),
                lambda: self._session.evolve_and_measure_tdvp(
                    h.to_terms(),operator.to_terms(),wf.cpp_handle,
                    int(nt),dt))
    elif self.itensor_version in (3,"python") and self.tevol_method=="TDVP_GSE":
        correlator,_wf = self._session.evolve_and_measure_tdvp_gse(
                h.to_terms(),operator.to_terms(),wf.cpp_handle,
                int(nt),dt,self.tdvp_gse_sweeps,self.tdvp_gse_krylov_order,
                self.tdvp_gse_cutoff)
    else:
        correlator,_wf = self._session.evolve_and_measure(
                h.to_terms(),operator.to_terms(),wf.cpp_handle,
                int(nt),dt,False)
    cs = np.array(correlator)
    ts = np.array([dt*ii for ii in range(int(nt))])
    if return_wf:
        from . import mps as mpsmod
        wf_final = mpsmod.MPS(self,cpp_handle=_wf).copy()
        return ts,cs.real-1j*cs.imag,wf_final
    return ts,cs.real-1j*cs.imag


def evolution_ABA(self,A=None,B=None,mode="DMRG",wf=None,**kwargs):
    """Apply an operator, evolve and measure"""
    if A is None: A = multioperator.identity()
    if B is None: B = multioperator.identity()
    if mode=="DMRG":
        if wf is None: wf = self.get_gs() # get ground state
        wfA = A*wf # apply the operator
        return evolve_and_measure_dmrg(self,wf=wfA,operator=B,**kwargs)
    elif mode=="ED":
        edobj = self.get_ED_obj() # get the ED object
        return tded.evolution_ABA(edobj,h=self.hamiltonian,A=A,B=B,wf=wf,
                **kwargs)






def dynamical_correlator(self,window=[-1,10],es=None,dt=0.1,
        nt=None,factor=1,delta=5e-2,damping_periods=6,damping="exp",
        predict=True,lp_order=None,lp_extend_factor=10,
        lp_fit_start_fraction=0.5,lp_max_pole_radius=1.0,
        **kwargs):
    """
    Compute a dynamical correlator from real-time evolution + Fourier
    transform (submode="TD", TDVP-backed for itensor_version=3).

    The raw finite-time correlator C(t) is windowed before the FFT
    (`damping`, see `_fourier_transform_correlator`'s docstring for the
    available choices and the tradeoff between them). The default,
    `damping="exp"`, applies exp(-delta*t); this is what actually turns
    `delta` into a Lorentzian broadening of width `delta` in the resulting
    spectral function -- matching what `delta` means in the KPM/CVM
    submodes -- and lets the required total evolution time follow directly
    from the decay itself: `damping_periods`/delta e-foldings of exp(-delta*t)
    make the truncation error exp(-damping_periods) negligible (default 6
    -> ~0.25%), instead of the previous undamped 100/delta default, which
    had no explicit broadening mechanism and relied on brute-force long
    evolution (and the resulting rectangular-window ringing) to get a
    comparably resolved spectrum. This cuts the number of time steps -- and
    thus wall-clock cost, the dominant cost of this submode -- by more than
    an order of magnitude for the same `delta`, making it competitive with
    KPM. The Fourier sum is normalized as a Riemann sum (factor dtnew) to
    match the analytic Fourier-transform convention of the other submodes,
    replacing the previous ad hoc 1/sqrt(nt) scaling that was tied to the
    old undamped/long-time convention.

    predict=True (the default) additionally extrapolates C(t) via linear
    prediction before windowing (see `_fourier_transform_correlator`'s
    `predict=` kwarg and dynamicstk/linearprediction.py) -- confirmed
    empirically (see docs/td_dynamical_correlator_sharpening_plan.md) to
    give a measurably narrower, better-centered peak than plain `"exp"`
    damping alone, at no extra real-TDVP cost, which is why it is the
    default rather than opt-in; pass predict=False to recover the old
    behavior exactly. `damping="exp"` (unchanged) stays the default taper
    -- pairing prediction with "gaussian" was checked too and came out
    *worse* (its wider intrinsic FWHM at fixed delta partly cancels
    prediction's own narrowing), so it is offered but not defaulted to.
    lp_order=None (default) auto-picks a safe AR order
    (`min(20, max(4, nt//10))`) so this stays robust even if `nt` ends up
    small (e.g. from an unusually large `delta` or an explicit small
    `nt`), rather than risking linear_predict_extend's own order-vs-length
    ValueError at the previous fixed default of 20.
    """
    self.get_gs() # get the ground state
    if nt is None: nt=int(damping_periods/delta/dt)
    if lp_order is None: lp_order=min(20,max(4,nt//10))
    (ts,cs) = evolution_DC(self,dt=dt,nt=nt,**kwargs) # get correlator
    return _fourier_transform_correlator(ts,cs,dt,es=es,window=window,
            delta=delta,factor=factor,damping=damping,predict=predict,
            lp_order=lp_order,lp_extend_factor=lp_extend_factor,
            lp_fit_start_fraction=lp_fit_start_fraction,
            lp_max_pole_radius=lp_max_pole_radius)


def _damping_window(ts,delta,damping="exp"):
    """
    Time-domain taper applied to C(t) before the FFT, selecting the
    lineshape/tail behavior of the resulting spectral function -- the
    time-domain analogue of `algebra/kpm.py`'s
    `kernel="jackson"/"lorentz"/"plain"` choice for the KPM submode.

    `delta` keeps the same "characteristic broadening width" meaning
    across all choices (matching KPM/CVM's own `delta`), but the shape of
    the resulting line differs:

    - "exp" (default, unchanged behavior): exp(-delta*t), the previous
      hardcoded choice. Exact Lorentzian broadening in frequency space,
      i.e. a `1/(omega-omega0)**2` algebraic tail -- this is the "long
      tail" reported for submode="TD" vs KPM's default Jackson-kernel
      reconstruction, which decays much faster away from a peak (see
      Weisse, Wellein, Alvermann & Fehske, Rev. Mod. Phys. 78, 275 (2006),
      and `algebra/kpm.py::jackson_kernel`).
    - "gaussian": exp(-(delta*t)**2/2), a Gaussian taper. Its Fourier
      transform is itself a Gaussian, decaying as exp(-omega**2) --
      dramatically faster far from the peak than the Lorentzian's
      algebraic 1/omega**2 tail -- at the cost of the usual Lorentzian-
      vs-Gaussian lineshape tradeoff: at the same `delta`, the Gaussian's
      FWHM (2*sqrt(2*ln(2))*delta ~ 2.35*delta) is actually slightly
      *wider* than the Lorentzian's (2*delta), i.e. the peak itself looks
      marginally broader/shorter even as the far tail is suppressed by
      orders of magnitude. Standard NMR/spectroscopy apodization choice;
      also offered by TeNPy's SpectralSimulation
      ("gaussian windowing", GPL-3.0, same license as this project;
      https://tenpy.readthedocs.io/en/v1.0.2/reference/tenpy.simulations.time_evolution.SpectralSimulation.html).
    - "parzen": the Parzen window, a smooth taper that goes to zero (with
      vanishing derivative) exactly at the truncation time `Tmax=max(ts)`,
      independent of `delta`'s decay rate. This targets a different
      artifact than the peak-broadening tradeoff above: the Gibbs ringing
      from abruptly truncating C(t) at a finite Tmax (the implicit
      rectangular window every choice here still has, since the FFT only
      ever sees `ts` up to Tmax) -- a taper that is exactly zero at both
      ends removes that discontinuity. Reported in the windowed-FT
      literature for real-time tensor-network correlators as the
      standard fix for that specific artifact, independent of the choice
      above (see docs/td_dynamical_correlator_sharpening_plan.md for the
      literature pointers). Still combined multiplicatively with "exp"'s
      exp(-delta*t) so `delta` keeps controlling the overall resolution.
    """
    if damping=="exp":
        return np.exp(-delta*ts)
    elif damping=="gaussian":
        return np.exp(-0.5*(delta*ts)**2)
    elif damping=="parzen":
        tmax = np.max(ts)
        if tmax<=0: return np.ones_like(ts)
        x = ts/tmax # in [0,1]
        w = np.where(x<=0.5,
                1.-6.*x**2*(1.-x),
                2.*(1.-x)**3)
        return w*np.exp(-delta*ts)
    else:
        raise ValueError("Unknown damping: "+str(damping))


def _fourier_transform_correlator(ts,cs,dt,es=None,window=[-1,10],
        delta=5e-2,factor=1,damping="exp",predict=False,lp_order=20,
        lp_extend_factor=10,lp_fit_start_fraction=0.5,
        lp_max_pole_radius=1.0):
    """
    Shared time-domain -> frequency-domain tail: optional linear-
    prediction extrapolation (`predict`), a damping/window taper (see
    `_damping_window`'s docstring for the available choices and the
    tradeoff between them), interpolation onto a uniform (optionally
    oversampled by `factor`) grid, a Riemann-sum-normalized FFT, and
    interpolation onto the requested frequencies `es`. Factored out of
    dynamical_correlator (submode "TD") so other time-domain submodes
    (e.g. "TDZ", see tdz.py) and `sxt_to_skomega` (per k-point) can reuse
    it unchanged instead of duplicating the FFT/windowing/extrapolation
    convention.

    predict=True runs `dynamicstk.linearprediction.linear_predict_extend`
    on the raw `(ts,cs)` first, extending the effective simulated time
    well beyond what was actually evolved -- done here, before damping,
    so the (now much longer) extrapolated series is what the damping
    window and FFT actually see (see
    docs/td_dynamical_correlator_sharpening_plan.md). `lp_order`/
    `lp_extend_factor`/`lp_fit_start_fraction`/`lp_max_pole_radius` are
    passed straight through to `linear_predict_extend` -- see its own
    docstring.
    """
    if predict:
        from .dynamicstk.linearprediction import linear_predict_extend
        ts,cs = linear_predict_extend(ts,cs,order=lp_order,
                extend_factor=lp_extend_factor,
                fit_start_fraction=lp_fit_start_fraction,
                max_pole_radius=lp_max_pole_radius)
    cs = cs*_damping_window(ts,delta,damping=damping) # damping/window taper
    # interpolate the time evolution
    ftr = interp1d(ts,cs.real,fill_value=0.0,bounds_error=False)
    fti = interp1d(ts,cs.imag,fill_value=0.0,bounds_error=False)
    # interpolate the time evolution
    tnew = np.linspace(np.min(ts),np.max(ts),len(ts)*factor) # ten times
    cnew = ftr(tnew) + 1j*fti(tnew)
    ts = tnew.copy() # overwrite
    cs = cnew.copy() # overwrite
    dtnew = dt/factor
    # do the fourier transform
    ss = np.fft.fft(cs)*dtnew # fourier transform (Riemann-sum normalization)
    ws = np.fft.fftfreq(len(cs),d=dtnew)*2.*np.pi # fourier frequencies
    fr = interp1d(ws,ss.real,fill_value=0.0,bounds_error=False)
    fi = interp1d(ws,ss.imag,fill_value=0.0,bounds_error=False)
    if es is None:
        es = np.linspace(window[0],window[1],800)
    gr = fr(es)+ 1j*fi(es) # advanced
    ga = np.conjugate(gr) # retarded
#    gp = fr(es) - fr(-es) + 1j*fi(es) + 1j*fi(-es)
    return (es,gr)


def sxt_to_skomega(ts,xs,S,dt,ks=None,es=None,window=[-1,10],
        delta=5e-2,factor=1,damping="exp",predict=False,lp_order=20,
        lp_extend_factor=10,lp_fit_start_fraction=0.5,
        lp_max_pole_radius=1.0):
    """S(k,omega) from a real-space/real-time correlator S(x,t) (`S`
    shaped `(len(ts),len(xs))`): a spatial DFT
    (`S(k,t)=sum_x e^{-ikx}S(x,t)`) followed by `_fourier_transform_correlator`
    -- factored out so every real-time-evolution-based dynamical-correlator
    submode that produces an `S(x,t)` array reduces it to `S(k,omega)` the
    same way, instead of duplicating this loop per submode/backend (this
    was previously inlined separately in both
    `pyitensor.idmrg_window.dynamical_correlator_komega` and
    `infinitechain.py`'s own `itensor_version=3` dispatch for
    `td_dynamical_correlator` -- confirmed via code review to be exact
    duplicates, now unified here so a future fix to this reduction can't
    silently apply to only one of them).

    `ks` defaults to 200 points in `[-pi,pi]` (the first Brillouin zone,
    since `x` is measured in physical sites); `es`/`window`/`delta`/
    `factor`/`damping`/`predict`/`lp_*` are passed straight through to
    `_fourier_transform_correlator` (applied independently per k-point)
    -- see its own docstring. Returns `(ks, es, Skw)`, `Skw` shaped
    `(len(ks), len(es))`."""
    if ks is None:
        ks = np.linspace(-np.pi,np.pi,200)
    ks = np.asarray(ks)
    xs = np.asarray(xs)

    Skw = None
    es_out = es
    for ik,k in enumerate(ks):
        phase = np.exp(-1j*k*xs)
        Skt = S@phase
        es_k,gk = _fourier_transform_correlator(ts,Skt,dt,es=es_out,
                                                  window=window,delta=delta,
                                                  factor=factor,damping=damping,
                                                  predict=predict,lp_order=lp_order,
                                                  lp_extend_factor=lp_extend_factor,
                                                  lp_fit_start_fraction=lp_fit_start_fraction,
                                                  lp_max_pole_radius=lp_max_pole_radius)
        if Skw is None:
            es_out = es_k
            Skw = np.zeros((len(ks),len(es_k)),dtype=complex)
        Skw[ik] = gk
    return ks,es_out,Skw



def generic_evolution(H,wf,normalize=True,dt=1e-2,nt=100,A=None):
    """Perform a time evolution and project onto itself,
    assuming U = e^tH """
    wf0 = wf.copy() # copy wavefunction
    wf1 = wf.copy() # copy wavefunction
    out = []
    for i in range(int(nt)): # loop
        wf1 = wf1 + dt*H*wf1
        if normalize:  wf1 = wf1*(1./np.sqrt(wf1.dot(wf1)))
  #      out.append(wf0.dot(wf1)) # compute
        out.append(wf1.dot(A*wf1)) # compute
        print(i)
    return np.array(range(int(nt)))*dt,np.array(out) # retunr result



