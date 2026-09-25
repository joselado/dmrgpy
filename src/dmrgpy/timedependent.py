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


TEVOL_METHODS = ("TDVP","TDVP_GSE","TEBD","AUTO","MPO")


def check_tevol_method(self):
    """Reject an unrecognized self.tevol_method instead of silently
    running a different integrator.

    Every dispatch below tests `itensor_version in (3,"python") and
    tevol_method=="<literal>"` and ends in a bare `else` that runs the
    legacy MPO-Taylor path. That else-branch is doing two jobs at once:
    it is the documented, intended fallback for a backend that cannot
    honor the request (itensor_version=2 has no TDVP/ module at all, see
    docs/user_guide.md), and it was also catching every misspelling. So
    tevol_method="tdvp" (wrong case), "TVDP", "TEBD " or any typo ran
    MPO-Taylor without a word: measured against ED on a 6-site chain,
    5.9e-5 error for the TDVP that was asked for versus 3.5e-2 for what
    was actually run, bit-identical to an explicit "MPO".

    Only the name is checked here. A recognized method that this backend
    cannot run still falls back exactly as documented -- that case is a
    property of the backend, the typo is a property of the caller."""
    if self.tevol_method not in TEVOL_METHODS:
        raise ValueError(
            "Unrecognized tevol_method "+repr(self.tevol_method)+"; expected "
            +", ".join(repr(m) for m in TEVOL_METHODS)
            +". (An unrecognized method used to run the legacy MPO-Taylor "
            "integrator silently.)")


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
    check_tevol_method(self) # reject a typo instead of running MPO-Taylor
    if self.itensor_version=="julia_live":
        from .mpsjulialive import timedependent as tdjl
        return tdjl.evolution_dmrg_DC(self,name=name,nt=nt,dt=dt,**kwargs)
    name = operatornames.str2MO(self,name,
            require_symbolic_for="submode='TD'",**kwargs)
    name[0] = name[0].get_dagger()
    A,B = name[0],name[1]
    from .groundstate import send_hamiltonian
    send_hamiltonian(self) # quench*/evolve_and_measure* start from the
    # session's own get_gs(), which aborts the process if set_hamiltonian
    # was never called -- reachable from a freshly built chain, see
    # groundstate.send_hamiltonian's docstring.
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



def evolve_and_measure(self,mode="DMRG",nt=1000,dt=1e-2,h=None,**kwargs):
    """Evolve and measure <psi(t)|operator|psi(t)>, psi(t) = e^{-iht}|wf>,
    with h the chain's Hamiltonian unless given.

    nt, dt and h are named here and forwarded to both modes, so one call
    gives the same time grid and the same Hamiltonian on DMRG and on ED:
    the ED route used to default to nt=100 against DMRG's 1000, and to
    raise "got multiple values for argument 'h'" on an h= that DMRG
    honours (2026-09-24c audit, finding 16)."""
    if mode=="DMRG": return evolve_and_measure_dmrg(self,nt=nt,dt=dt,h=h,**kwargs)
    elif mode=="ED": 
        edobj = self.get_ED_obj() # get the ED object
        if h is None: h = self.hamiltonian
        return tded.evolve_and_measure(edobj,h,nt=nt,dt=dt,**kwargs)



def evolve_and_measure_dmrg(self,operator=None,nt=1000,h=None,
        dt=1e-2,wf=None,return_wf=False):
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

    What comes back is <psi(t)|O|psi(t)>, psi(t) = e^{-iHt}|wf>, exactly
    the list the session measures (every session method takes
    <psi|O|psi>), so at t=0 it is vev(O) on the same state, and for a
    non-Hermitian O it keeps the sign of its imaginary part. It used to
    be returned conjugated, a line copied from evolution_dmrg_DC above,
    which made this <psi(t)|O^dagger|psi(t)>: invisible on a Hermitian O,
    whose imaginary part is roundoff, and exactly minus the imaginary part
    otherwise, -0.5i at t=0 for O = Sz_0 + i*Sx_0 on the +x state where
    vev(O) is +0.5i, on every backend and integrator (2026-09-24b audit,
    finding 10). evolution_dmrg_DC keeps its conjugation, and for a
    different reason: it returns a time correlator, not an expectation
    value, and conjugating the session's <GS|B^dagger e^{-i(H-E_0)t}
    A^dagger|GS> is what gives sum_n M_n e^{+i D_n t}, the series whose
    one-sided transform puts the lines of dynamics.py's house convention
    at omega = +D_n.

    There is no **kwargs: it used to take one that nothing read, so a
    misspelled keyword (DT=0.2) ran silently at the defaults, dt=1e-2
    (2026-09-24c audit, finding 15), where mode="ED" raises.
    """
    check_tevol_method(self) # reject a typo instead of running MPO-Taylor
    if self.itensor_version=="julia_live":
        from .mpsjulialive import timedependent as tdjl
        return tdjl.evolve_and_measure_dmrg(self,operator=operator,nt=nt,
                h=h,dt=dt,wf=wf,return_wf=return_wf)
    if h is None: h = self.hamiltonian # Hamiltonian
    if wf is None: wf = self.wf0 # get ground state
    from .groundstate import send_hamiltonian
    send_hamiltonian(self) # quench*/evolve_and_measure* start from the
    # session's own get_gs(), which aborts the process if set_hamiltonian
    # was never called -- reachable from a freshly built chain, see
    # groundstate.send_hamiltonian's docstring.
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
    # <psi(t)|O|psi(t)> as measured, not conjugated: see the docstring for
    # why evolution_dmrg_DC conjugates and this does not
    cs = np.array(correlator)
    ts = np.array([dt*ii for ii in range(int(nt))])
    if return_wf:
        from . import mps as mpsmod
        wf_final = mpsmod.MPS(self,cpp_handle=_wf).copy()
        return ts,cs,wf_final
    return ts,cs


def evolution_ABA(self,A=None,B=None,mode="DMRG",wf=None,nt=1000,dt=1e-2,
        h=None,**kwargs):
    """Apply A, evolve and measure: <wf|A^dagger e^{iht} B e^{-iht} A|wf>.
    nt, dt and h are forwarded to both modes, as in evolve_and_measure."""
    if A is None: A = multioperator.identity()
    if B is None: B = multioperator.identity()
    if mode=="DMRG":
        if wf is None: wf = self.get_gs() # get ground state
        wfA = A*wf # apply the operator
        return evolve_and_measure_dmrg(self,wf=wfA,operator=B,nt=nt,dt=dt,
                h=h,**kwargs)
    elif mode=="ED":
        edobj = self.get_ED_obj() # get the ED object
        if h is None: h = self.hamiltonian
        return tded.evolution_ABA(edobj,h=h,A=A,B=B,wf=wf,nt=nt,dt=dt,
                **kwargs)






def lehmann_density_from_one_sided(self,name,transform,i=0,j=0):
    """Assemble the house dynamical correlator from one-sided real-time
    transforms, for the two submodes that are built on one.

    `transform(pair)` must return `(es,y)` for the operator pair it is
    given, where `y` is the one-sided Fourier transform
    `(1/pi) int_0^inf dt e^{-i w t} C(t)` that
    `_fourier_transform_correlator` produces; `name` is the pair the
    caller asked for, in any form `operatornames.str2MO` understands,
    with `i`/`j` the sites of a string name. Both are resolved here, so
    the two have to reach this function: they used not to, and
    get_dynamical_correlator_MB(name="ZZ",i=1,j=1) returned C[Sz_0,Sz_0]
    bit for bit on both submodes (2026-09-24 audit, finding 10).

    What this is for. The house quantity (see dynamics.py's module
    docstring) is the complex Lehmann density
    `C_AB(w) = sum_n M_n L_delta(w-D_n)`, which is the *two*-sided
    transform of the time correlator `C(t) = sum_n M_n e^{-i D_n t}`,
    damped by `e^{-delta|t|}`. A real-time run only ever produces the
    `t>=0` half, and that half alone is a resolvent: its real part is
    `C_AB` when every `M_n` is real, while its imaginary part is the
    dispersive term `C_AB` does not have, measured at 70% of the
    correlator's own peak even on a Hermitian pair. That is what made
    submode="TD"/"TDZ" the two routes off the convention, recorded as
    open item O1 of the 2026-09 audit.

    The missing half is not a second simulation backwards in time. For
    `t>0`,

        C(-t) = <GS|B e^{+i(H-E_0)t} A|GS>
              = conj( <GS|A^dagger e^{-i(H-E_0)t} B^dagger|GS> ),

    so the backward half of the pair `(A,B)` is the conjugate of the
    *forward* half of the pair `(B^dagger, A^dagger)`, which the same
    machinery computes with no change at all. With the transform kernel
    satisfying `K(-t) = conj(K(t))` for real w, the two halves combine as

        C_AB = ( F[(A,B)] + conj(F[(B^dagger,A^dagger)]) ) / 2,

    and the 1/pi the one-sided transform already carries turns into the
    1/(2*pi) the two-sided one needs. When `A` is provably `B^dagger`
    the adjoint pair *is* the original pair and the formula collapses to
    `Re F`, so that case costs one evolution rather than two; this is
    every example in the documentation. `canonical.is_dagger_pair` is
    the test, and it refuses rather than guesses when an operator name
    has no known adjoint (a parafermionic `Sig`, a caller's own name),
    which costs a second evolution. That second evolution is built from
    `get_dagger()`, so it is right exactly when `get_dagger()` knows the
    name's adjoint: a name it leaves untouched is treated as Hermitian,
    which for an anti-Hermitian one (the raw backend name `ISy`, i*Sy,
    while get_dagger() passed it through unchanged) returns exactly minus
    the correlator (2026-09-24 audit, finding 2). Refusing the proof
    protects the one-run shortcut, not the adjoint the second run is
    built on.

    Measured on the 2026-09 audit's own seeded 4-site complex-hopping
    chain (`A = Cdag_0`, `B = C_2`, `max|Im M_n| = 0.27`, exact peak
    0.2313, delta=0.4, dt=0.1), against an exact Lehmann sum built by
    dense diagonalization outside dmrgpy:

        raw one-sided transform (what this replaced)  1.16e-01
        its real part alone                           2.15e-01
        this combination                              2.57e-04

    i.e. the real part alone is not a partial fix on a complex-weight
    pair, it is worse than leaving the transform whole, while the
    combination lands at the TD method's own discretization error. These
    three are submode="TD" at its default predict=True, with the
    frequency stage of the time; since that stage evaluates the sum at
    each requested frequency instead of interpolating an FFT grid (see
    _fourier_transform_correlator), the combination is at 3.2e-05, with
    TD at predict=False and TDZ both at 5.4e-04."""
    # No require_symbolic_for here, deliberately: mode="ED" consumes the
    # already-built operators toMPO(mode="ED") returns, and always has.
    # The DMRG side rebuilds its operators from to_terms() inside the
    # backend and rejects a compiled one in evolution_dmrg_DC, which is
    # where that restriction belongs and where it is pinned.
    pair = operatornames.str2MO(self,name,i=i,j=j)
    (xs,ys) = transform([pair[0],pair[1]])
    if _pair_is_self_adjoint(pair):
        # the adjoint pair is this pair, so the second run would return
        # the same array and the combination below is exactly its real part
        return xs,np.asarray(ys).real+0.0j
    (_xs,ys2) = transform(_adjoint_pair(pair))
    return xs,(np.asarray(ys)+np.conjugate(ys2))/2.0


def _dagger_operator(o):
    """The adjoint of whatever `name=` was given, symbolic or compiled."""
    if hasattr(o,"get_dagger"): return o.get_dagger()
    so = getattr(o,"SO",None) # a compiled ED operator carries its matrix
    if so is not None:
        out = o.copy(); out.SO = so.conj().transpose(); return out
    raise TypeError(
        "submode='TD'/'TDZ': cannot take the adjoint of an operator of "
        "type "+type(o).__name__+", which the correlator needs for a pair "
        "that is not its own adjoint (see "
        "lehmann_density_from_one_sided). Pass the MultiOperators "
        "themselves rather than a compiled operator.")


def _adjoint_pair(pair):
    """(B^dagger, A^dagger): the pair whose forward-time correlator is the
    conjugate of this pair's backward-time one."""
    return [_dagger_operator(pair[1]),_dagger_operator(pair[0])]


def _pair_is_self_adjoint(pair):
    """True when A is provably B^dagger, so that the two-term combination
    collapses to a real part and one evolution suffices. False is "not
    proven" and costs a second evolution, which is built from
    get_dagger() and so is right exactly when get_dagger() knows the
    adjoint of every name in the pair (see
    lehmann_density_from_one_sided).

    Symbolically this is canonical.is_dagger_pair, which refuses rather
    than guesses for a name with no known adjoint. A compiled operator
    has no canonical form but does carry its own matrix, so there the
    same question is a numerical one and is answered exactly."""
    from .multioperatortk import canonical
    from .multioperator import MultiOperator
    A,B = pair[0],pair[1]
    if isinstance(A,MultiOperator) and isinstance(B,MultiOperator):
        return canonical.is_dagger_pair(A,B)
    sa,sb = getattr(A,"SO",None),getattr(B,"SO",None)
    if sa is not None and sb is not None:
        try: return abs(sa-sb.conj().transpose()).max()<1e-12
        except Exception: return False
    return False


def dynamical_correlator(self,window=[-1,10],es=None,dt=0.1,
        nt=None,factor=1,delta=5e-2,damping_periods=6,damping="exp",
        predict=True,lp_order=None,lp_extend_factor=10,
        lp_fit_start_fraction=0.5,lp_max_pole_radius=1.0,
        **kwargs):
    """
    Compute a dynamical correlator from real-time evolution + Fourier
    transform (submode="TD", TDVP-backed for itensor_version=3).

    The raw finite-time correlator C(t) is windowed before the Fourier sum
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
    behavior exactly. That comparison was made on the old frequency
    stage, and at the default time window (delta*T = damping_periods =
    6) most of the narrowing it saw was that stage's FFT grid, which
    prediction's tenfold longer series made ten times finer (see
    `_fourier_transform_correlator`). Re-measured on the same 4-site
    Heisenberg chain at delta=0.05: at the default window the line is at
    the exact width with or without prediction (FWHM 0.0975 on a 0.0025
    grid, as the exact density), and what prediction still buys there is
    accuracy, max|y - exact| 2.55e-03 -> 3.3e-06, by carrying the series
    past the e^-6 cut; where the window is short enough that truncation
    sets the width it also narrows the line, 0.2125 -> 0.0975 at nt=200
    (delta*T=1), 0.3875 -> 0.0975 at nt=100, and moves a displaced peak
    back onto the gap. `damping="exp"` (unchanged) stays the default taper
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
    name = kwargs.pop("name","XX")
    # i/j are the sites of a string name, and lehmann_density_from_one_sided
    # is where the name is resolved, so they go there and nowhere else.
    # Only these two are popped: every other keyword rides on into
    # evolution_DC, and on the DMRG side into str2MO, which is what makes a
    # misspelled one raise TypeError instead of being ignored.
    i,j = kwargs.pop("i",0),kwargs.pop("j",0)
    def transform(pair):
        (ts,cs) = evolution_DC(self,dt=dt,nt=nt,name=pair,**kwargs)
        return _fourier_transform_correlator(ts,cs,dt,es=es,window=window,
                delta=delta,factor=factor,damping=damping,predict=predict,
                lp_order=lp_order,lp_extend_factor=lp_extend_factor,
                lp_fit_start_fraction=lp_fit_start_fraction,
                lp_max_pole_radius=lp_max_pole_radius)
    # one evolution for a pair whose adjoint is itself, two otherwise --
    # see lehmann_density_from_one_sided for the identity and what this
    # used to return instead
    return lehmann_density_from_one_sided(self,name,transform,i=i,j=j)


def _damping_window(ts,delta,damping="exp"):
    """
    Time-domain taper applied to C(t) before the Fourier sum, selecting the
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
      rectangular window every choice here still has, since the Fourier sum only
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


# Upper bound on the number of entries of the exp(-i w t) block the direct
# frequency evaluation below builds at once (2**21 complex entries is 32
# MB); the requested frequencies are taken in chunks of that size, so the
# memory stays bounded whatever len(es)*len(ts) is.
_DIRECT_FT_MAX_ELEMENTS = 2**21


def _damped_sum_at(cs,dtnew,es):
    """(dtnew/pi) sum_k cs[k] exp(-i w k dtnew), evaluated at every w in
    `es`, chunked over `es`. The time origin is the first sample, k=0,
    which is the phase the FFT puts on sample k, so at a frequency on the
    FFT grid the two agree to rounding. Outside the band the FFT covers,
    [min(fftfreq), max(fftfreq)], the result is 0, as the interpolation
    of the FFT returned there: out of band the sum only repeats an
    in-band frequency (aliasing), it carries no new information."""
    es_arr = np.asarray(es,dtype=float)
    ws = es_arr.reshape(-1)
    tt = dtnew*np.arange(len(cs))
    gr = np.zeros(len(ws),dtype=np.complex128)
    rows = max(1,_DIRECT_FT_MAX_ELEMENTS//max(1,len(tt)))
    for s in range(0,len(ws),rows):
        gr[s:s+rows] = np.exp(-1j*np.outer(ws[s:s+rows],tt))@cs
    gr = gr*dtnew/np.pi
    n = len(cs)
    wmin = -(n//2)*2.*np.pi/(n*dtnew)
    wmax = ((n-1)//2)*2.*np.pi/(n*dtnew)
    gr[(ws<wmin) | (ws>wmax)] = 0.0
    return gr.reshape(es_arr.shape)


def _fourier_transform_correlator(ts,cs,dt,es=None,window=[-1,10],
        delta=5e-2,factor=1,damping="exp",predict=False,lp_order=20,
        lp_extend_factor=10,lp_fit_start_fraction=0.5,
        lp_max_pole_radius=1.0,_evaluation="direct"):
    """
    Shared time-domain -> frequency-domain tail: optional linear-
    prediction extrapolation (`predict`), a damping/window taper (see
    `_damping_window`'s docstring for the available choices and the
    tradeoff between them), interpolation onto a uniform (optionally
    oversampled by `factor`) grid, and the trapezoid sum of the damped
    series evaluated at each requested frequency in `es`. Factored out of
    dynamical_correlator (submode "TD") so other time-domain submodes
    (e.g. "TDZ", see tdz.py) and `sxt_to_skomega` (per k-point) can reuse
    it unchanged instead of duplicating the Fourier/windowing/extrapolation
    convention.

    predict=True runs `dynamicstk.linearprediction.linear_predict_extend`
    on the raw `(ts,cs)` first, extending the effective simulated time
    well beyond what was actually evolved -- done here, before damping,
    so the (now much longer) extrapolated series is what the damping
    window and the Fourier sum actually see (see
    docs/td_dynamical_correlator_sharpening_plan.md). `lp_order`/
    `lp_extend_factor`/`lp_fit_start_fraction`/`lp_max_pole_radius` are
    passed straight through to `linear_predict_extend` -- see its own
    docstring.

    How the frequencies are evaluated. The sum used to be taken with an
    FFT, i.e. only on the grid of spacing 2*pi/(nt*dt), and then
    interpolated linearly onto `es`. At the default time window that
    spacing is 2*pi*delta/damping_periods, 1.05*delta at
    damping_periods=6, so a Lorentzian of half-width delta was sampled
    about once per delta and the interpolation between samples was the
    dominant error of every predict=False spectrum: on the 2026-09
    audit's seeded 4-site complex-hopping chain (A=Cdag_0, B=C_2,
    delta=0.4, dt=0.1, exact peak 0.2313) it was 3.31e-02 for TDZ at its
    defaults and for TD at predict=False alike, on DMRG and on ED, an
    error the records had attributed to TDZ's complex-time contour, whose
    own share is 1e-06 (2026-09-24 audit, finding 7). The trapezoid sum
    is now evaluated at each requested frequency directly
    (`_damped_sum_at`), which at a frequency on the old grid gives the
    FFT's value to rounding and in between gives the sum itself, so the
    error left is the one set by the finite time window and the method:
    5.4e-04 on that chain for TDZ and TD at predict=False, and 3.2e-05
    for TD at its default predict=True, which went through the same
    interpolation on a grid ten times finer (2.57e-04 before).

    `_evaluation="fft"` keeps the old FFT-plus-interpolation stage as a
    reference, and no production route calls it any more: its one test
    (tests/test_audit_2026_09_24_realtime.py::
    test_direct_evaluation_is_the_fft_on_its_own_grid) checks the direct
    sum against it on the FFT's own grid. The infinite-chain
    `sxt_to_skomega` was the last caller, deferred on the grounds that the
    direct sum was unmeasured at the short windows its defaults give; the
    deferral did not depend on the window, so a caller who converged it
    (delta*T = 6) still got the interpolation, 2.1e-01 of the peak off
    the exact damped transform of a single-magnon series with peak
    heights down to 0.79 of exact, where the direct sum is 2.5e-04 off
    it and equals the closed-form damped trapezoid sum to 1e-13
    (2026-09-24b audit, finding 7). At the short windows (delta*T of 1
    and below) neither stage is converged, and what that regime needs is
    a longer `nt`, not a different frequency stage.
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
    # Trapezoidal, not rectangular: giving the t=0 sample its full weight
    # dtnew (as a plain Riemann sum does) leaves a real, frequency-
    # independent offset Re C(0)*dtnew/2 across the whole spectrum --
    # measured at exactly that value for dt=0.1/0.05/0.025 on a 4-site
    # Heisenberg chain. The tail sample is halved for the same reason;
    # after damping it is ~0 anyway, so that half costs nothing.
    cs = cs.copy()
    cs[0] = cs[0]*0.5
    cs[-1] = cs[-1]*0.5
    # do the fourier transform. The 1/pi is what makes the two-sided
    # combination in lehmann_density_from_one_sided come out with the
    # 1/(2*pi) of the house convention, the complex Lehmann density C_AB
    # (dynamics.py's module docstring), whose sum rule int C_AB dw = <A B>
    # every other submode satisfies. It used to be missing, and "TD" (and
    # "TDZ", which shares this tail) came out a factor of pi too large --
    # 1.284 against an exact 0.250 on a 4-site chain. Uniform in w, so no
    # peak position or width ever moved and nothing in tests/ could see
    # it; it showed up only when the absolute weight was integrated.
    if es is None:
        es = np.linspace(window[0],window[1],800)
    if _evaluation=="direct": # the sum itself at each es, see the docstring
        return (es,_damped_sum_at(cs,dtnew,es))
    if _evaluation!="fft":
        raise ValueError("Unknown _evaluation: "+repr(_evaluation))
    ss = np.fft.fft(cs)*dtnew/np.pi # fourier transform (trapezoid + 1/pi)
    ws = np.fft.fftfreq(len(cs),d=dtnew)*2.*np.pi # fourier frequencies
    fr = interp1d(ws,ss.real,fill_value=0.0,bounds_error=False)
    fi = interp1d(ws,ss.imag,fill_value=0.0,bounds_error=False)
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
    `(len(ks), len(es))`.

    The sign of the frequency. `S(x,t)` comes in as the infinite-window
    routes produce it, `<psi|A_x e^{-i(H-E_0)t} B_0|psi> = sum_n M_n(x)
    e^{-i D_n t}`, and the momentum series is conjugated after the
    spatial sum, `conj(S(k,t))`, before the time transform. That is the
    step the finite TD route takes at the end of evolution_dmrg_DC, which
    returns `sum_n M_n e^{+i D_n t}`, and it is what puts the lines at
    omega = +D_n, D_n = E_n - E_0 > 0, as on every finite route and in
    dynamics.py's house convention. Without it every infinite-chain
    S(k,omega), on both backends, came out mirrored, with the lines at
    -D_n and mostly below the default window: on the transverse-field
    paramagnet 1.4*Sz + Sx*Sx with A = B = Sx, the magnon at k=pi, where
    eps = 0.9, peaked at omega = -1.045 with 0.928 of its weight below
    zero (2026-09-24b audit, finding 8). The conjugation goes after the
    sum and not on each x series before it: that would be the transform
    of the conjugate at -k, which only a parity-symmetric model cannot
    tell apart, and on a chiral free-fermion ring it puts the occupied
    band at the wrong k.

    What comes back is still the one-sided transform, not the two-sided
    density `lehmann_density_from_one_sided` assembles for submode="TD",
    since this reduction never sees the operator pair. Its real part is
    the house density `sum_n W_n(k) delta/(pi((omega-D_n)^2+delta^2))`,
    `W_n(k) = sum_x e^{-ikx} M_n(x)`, whenever the weights `W_n(k)` are
    real, which a pair with `A_x = B_x^dagger` guarantees on a
    translation-invariant state when `xs` is symmetric about 0 (the
    default `x_values` of both routes) or runs over a full period; on a
    one-sided or ragged `xs` they are complex even for such a pair. Its
    imaginary part is the dispersive term the density does not have.

    The frequency stage is the direct per-frequency sum, the same as on
    every finite route. It used to be hard-wired to the FFT-plus-
    interpolation stage, whatever `nt` the caller gave, which at a
    converged window (delta*T = 6) is 2.1e-01 of the peak off exact on a
    single-magnon series, see `_fourier_transform_correlator`'s docstring
    (2026-09-24b audit, finding 7)."""
    if ks is None:
        ks = np.linspace(-np.pi,np.pi,200)
    ks = np.asarray(ks)
    xs = np.asarray(xs)

    Skw = None
    es_out = es
    for ik,k in enumerate(ks):
        phase = np.exp(-1j*k*xs)
        # conjugate after the spatial sum, see the docstring
        Skt = np.conj(S@phase)
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



