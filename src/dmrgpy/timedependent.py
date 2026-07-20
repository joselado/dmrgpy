from __future__ import print_function
from . import operatornames
import numpy as np
from scipy.interpolate import interp1d
from . import multioperator
from .edtk import timedependent as tded



def evolution_DC(self,mode="DMRG",**kwargs):
    if mode=="DMRG":  return evolution_dmrg_DC(self,**kwargs)
    if mode=="ED": 
        edobj = self.get_ED_obj() # get the ED object
        return tded.evolution_DC(edobj,h=self.hamiltonian,**kwargs)



def evolution_dmrg_DC(self,name="XX",nt=10000,dt=0.1,restart=True,**kwargs):
    """
    Real-time quench dynamical correlator via the in-process pybind11
    extension.

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
    if h is None: h = self.hamiltonian # Hamiltonian
    if wf is None: wf = self.wf0 # get ground state
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    if self.itensor_version in (3,"python") and self.tevol_method=="TDVP":
        correlator,_wf = self._session.evolve_and_measure_tdvp(
                h.to_terms(),operator.to_terms(),wf.cpp_handle,
                int(nt),dt)
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
        nt=None,factor=1,delta=5e-2,damping_periods=6,**kwargs):
    """
    Compute a dynamical correlator from real-time evolution + Fourier
    transform (submode="TD", TDVP-backed for itensor_version=3).

    The raw finite-time correlator C(t) is windowed with an exponential
    decay exp(-delta*t) before the FFT. This is what actually turns
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
    """
    self.get_gs() # get the ground state
    if nt is None: nt=int(damping_periods/delta/dt)
    (ts,cs) = evolution_DC(self,dt=dt,nt=nt,**kwargs) # get correlator
    return _fourier_transform_correlator(ts,cs,dt,es=es,window=window,
            delta=delta,factor=factor)


def _fourier_transform_correlator(ts,cs,dt,es=None,window=[-1,10],
        delta=5e-2,factor=1):
    """
    Shared time-domain -> frequency-domain tail: exponential-decay
    windowing (turns `delta` into a Lorentzian broadening, see
    dynamical_correlator's docstring), interpolation onto a uniform
    (optionally oversampled by `factor`) grid, a Riemann-sum-normalized
    FFT, and interpolation onto the requested frequencies `es`. Factored
    out of dynamical_correlator (submode "TD") so other time-domain
    submodes (e.g. "TDZ", see tdz.py) can reuse it unchanged instead of
    duplicating the FFT/windowing convention.
    """
    cs = cs*np.exp(-delta*ts) # damping window -> Lorentzian broadening "delta"
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



