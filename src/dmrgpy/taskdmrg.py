"""submode="TASK": the tangent-space Krylov dynamical correlator
(Kovalska, von Delft and Gleis, arXiv:2510.07279).

The Lanczos recursion runs in the one-site tangent space of the ground-state
MPS rather than on compressed MPS, so the Krylov basis stays orthogonal and
the result is a set of discrete poles w_alpha with weights M_alpha, broadened
here with the house Lorentzian (dynamics.py's module docstring). The method
itself is in pyitensor/task.py, together with its one approximation, the
projection of H onto the tangent space, and the two-site variance that
measures it pole by pole. On itensor_version="python", 3 and
"julia_live"; v3 and julia_live hand their ground-state MPS and MPO site
tensors over as numpy arrays (julia_live through mpsjulialive/arrays.jl)
and the same Python code runs on them.
"""
import numpy as np
from . import operatornames


BACKENDS = ("python",3,"julia_live")

def _check_backend(self,what="get_dynamical_correlator"):
    if self.itensor_version not in BACKENDS:
        raise NotImplementedError(
            "%s: submode=\"TASK\" is implemented on itensor_version="
            "\"python\", 3 and \"julia_live\" only (got itensor_version=%r)"
            %(what,self.itensor_version))


def ground_state(self):
    """Make sure the chain's current ground state is the one the arrays
    are exported from: groundstate.ground_state_on_session on the session
    backends, get_gs() on julia_live, which has no session."""
    if self.itensor_version=="julia_live":
        self.get_gs()
        return
    from .groundstate import ground_state_on_session
    ground_state_on_session(self)


def _julia_mpo(self,op):
    """The Julia MPO of a MultiOperator (mpsjulialive/mpo.py's get_MPO, the
    route of every julia_live operator), or of a term list from
    to_terms(), whose names are already Jordan-Wigner transformed and are
    passed to the same AutoMPO reader as they are."""
    from .multioperator import MultiOperator
    from .mpsjulialive.juliasession import Main, to_julia_strvec
    if isinstance(op,MultiOperator):
        from .mpsjulialive.mpo import get_MPO
        return get_MPO(op,self)
    if len(op)==0: raise ValueError("submode=\"TASK\": an operator with "
            "no terms")
    ls = [str(len(op))]
    for c,factors in op:
        c = complex(c)
        ls.append(str(len(factors)))
        ls.append("  ".join([repr(c.real),repr(c.imag)]
                  +["%s  %d"%(name,int(i)) for name,i in factors])+"  ")
    return Main.toMPO(to_julia_strvec(ls),self.jlsites)


def _julia_arrays(jlarrays):
    return [np.ascontiguousarray(np.asarray(a),dtype=complex)
            for a in jlarrays]


def _julia_state(self):
    """The julia_live ground state, and its site arrays
    (mpsjulialive/arrays.jl), exported as they are, as on v3."""
    from .mpsjulialive.juliasession import Main
    wf = self.wf0
    if wf is None or getattr(wf,"jlmps",None) is None:
        raise RuntimeError("submode=\"TASK\" called before gs_energy")
    return _julia_arrays(Main.mps_site_arrays(wf.jlmps)),wf


def _julia_vev(self,op,wf):
    """<wf|O|wf>/<wf|wf> on julia_live."""
    from .mpsjulialive.juliasession import Main
    num = complex(Main.mps_mpo_mps(wf.jlmps,_julia_mpo(self,op),wf.jlmps))
    return num/complex(Main.inner(wf.jlmps,wf.jlmps)).real


def _julia_mpo_arrays(self,op):
    from .mpsjulialive.juliasession import Main
    return _julia_arrays(Main.mpo_site_arrays(_julia_mpo(self,op)))


def _v3_arrays(session,terms):
    return session.mpo_site_arrays(session.build_operator(terms))


def _python_state(s):
    """The ground state of a "python" session as site arrays, and the
    state itself (the arrays of the chain's state, chain._mps_arrays_lpr)."""
    from .pyitensor.chain import _mps_arrays_lpr
    if s.wf0 is None:
        raise RuntimeError("submode=\"TASK\" called before gs_energy")
    return _mps_arrays_lpr(s.wf0),s.wf0


def session_poles(self,terms_i,terms_j,terms_idag,nkry,variance=False,
        terms_ij=None,pad_bonds=None,expand_bonds=None,expand_rounds=1):
    """(w, weights, elastic, variances, dim) on the chain's session, which
    must hold the ground state (groundstate.ground_state_on_session).
    With terms_ij, the terms of the product A*B, a sixth entry is
    <GS|A B|GS> on the same state, the sum rule of the weights.
    pad_bonds (None, an int or "full") pads the bonds of the state before
    the recursion, expand_bonds (None or an int) expands them with
    directions from the two-site discarded space of B|GS> and of the first
    Krylov step (pyitensor/task.py's module docstring); not both.

    Both backends call task.correlator_poles here with their exported
    arrays: on "python" the arrays of the chain's state
    (chain._mps_arrays_lpr), on 3 those of the session's own export, on
    julia_live those of mpsjulialive/arrays.jl. On julia_live the operator
    arguments may be MultiOperators instead of terms (_julia_mpo)."""
    s = self._session
    from .pyitensor import task
    task.check_expansion_keywords(pad_bonds,expand_bonds,expand_rounds)
    if self.itensor_version=="julia_live":
        Ms,wf = _julia_state(self)
        Ws = _julia_mpo_arrays(self,self.hamiltonian)
        Wb = _julia_mpo_arrays(self,terms_j)
        Wadag = (None if terms_idag is None
                 else _julia_mpo_arrays(self,terms_idag))
        ws,wt,var,dim = task.correlator_poles(Ms,Ws,Wb,Wadag,nkry,variance,
                pad_bonds=pad_bonds,expand_bonds=expand_bonds,
                expand_rounds=expand_rounds)
        elastic = _julia_vev(self,terms_i,wf)*_julia_vev(self,terms_j,wf)
        if terms_ij is None: return ws,wt,complex(elastic),var,dim
        return (ws,wt,complex(elastic),var,dim,
                complex(_julia_vev(self,terms_ij,wf)))
    if self.itensor_version=="python":
        Ms,wf = _python_state(s)
        Ws = task.mpo_arrays(s.H)
        Wb = task.mpo_arrays(s._mpo(terms_j))
        Wadag = (None if terms_idag is None
                 else task.mpo_arrays(s._mpo(terms_idag)))
    else:
        wf = s.gs_wavefunction()
        Ms = s.mps_site_arrays(wf)
        Ws = s.mpo_site_arrays(s.hamiltonian_mpo())
        Wb = _v3_arrays(s,terms_j)
        Wadag = None if terms_idag is None else _v3_arrays(s,terms_idag)
    ws,wt,var,dim = task.correlator_poles(Ms,Ws,Wb,Wadag,nkry,variance,
            pad_bonds=pad_bonds,expand_bonds=expand_bonds,
            expand_rounds=expand_rounds)
    elastic = s.vev(terms_i,wf)*s.vev(terms_j,wf)
    if terms_ij is None: return ws,wt,complex(elastic),var,dim
    return ws,wt,complex(elastic),var,dim,complex(s.vev(terms_ij,wf))


# A "python" ground state of a real H comes in a complex gauge (its random
# start and its tensors are complex), so its TASK spectrum ran in complex
# arithmetic all the way: task.real_starts gave up at its first line.
# session_spectrum rotates the export copy to a real gauge before the
# tangent space is built (real_gauge_state), and the basis is then real.
# REAL_GAUGE = False keeps the complex gauge, to compare the two on one
# state.
REAL_GAUGE = True
# The largest part of the state, relative to its norm, that the rotation
# may drop; above it the complex gauge is kept.
REAL_GAUGE_DROP = 1e-6


def _mps_overlap(As,Bs):
    """<A|B> of two MPS given by their (left, d, right) site arrays."""
    E = np.ones((1,1),dtype=complex)
    for A,B in zip(As,Bs):
        E = np.tensordot(np.tensordot(E,np.conj(A),axes=(0,0)),B,
                         axes=([0,1],[0,1]))
    return complex(E[0,0])


def real_gauge_state(Ms,Ws):
    """(arrays, H arrays, dropped): the site arrays Ms of a state rotated
    to a real gauge (nrgmps.real_gauge, weighted by the bond matrices) and
    the MPO arrays Ws of H in a real gauge (task.real_mpo, the same
    operator to rounding), and the norm of the part of the state the
    rotation dropped, relative to its norm, sqrt(1 - |<real|psi>|^2) for
    the normalized states. That part is the imaginary part the solve
    left, since a solved state is real only to the accuracy of the solve;
    it is resolved down to about 1e-8, the square root of the rounding of
    the overlap. Arrays that are real already give (Ms, Ws, 0.0).
    (Ms, Ws, None), the arrays as given, keeps the complex gauge: with
    REAL_GAUGE off, for a complex H, for a state that is not real up to a
    phase (a member of a degenerate manifold is in general a complex
    mixture of two real states, and real_gauge refuses it), and for a
    dropped part above REAL_GAUGE_DROP."""
    from .pyitensor import backend
    host = [np.asarray(backend.to_host(M)) for M in Ms]
    if all(not np.iscomplexobj(M) or not np.any(M.imag) for M in host):
        return Ms,Ws,0.0
    if not REAL_GAUGE: return Ms,Ws,None
    from .pyitensor import task
    from .nrgmps import real_gauge
    Wr = task.real_mpo(Ws)
    if Wr is None: return Ms,Ws,None # a complex H
    As,_,Cs = task.mixed_gauge(host)
    As = [np.asarray(backend.to_host(A)) for A in As]
    Cs = [np.asarray(backend.to_host(C)) for C in Cs]
    try:
        Mr,_ = real_gauge(As,Cs)
    except RuntimeError: # not real up to a phase
        return Ms,Ws,None
    kept = abs(_mps_overlap(Mr,host))**2/(_mps_overlap(Mr,Mr).real
                                         *_mps_overlap(host,host).real)
    dropped = float(np.sqrt(max(0.,1.-kept)))
    if dropped>REAL_GAUGE_DROP: return Ms,Ws,None
    return Mr,Wr,dropped


def _exported_real(Ms):
    """0.0 when the exported arrays are real, None otherwise: v3 and
    julia_live hand their state over as it is."""
    return 0.0 if all(not np.iscomplexobj(M) or not np.any(np.asarray(M).imag)
                      for M in Ms) else None


# How many spectra a chain keeps (SPECTRUM_CACHE_SIZE of them, the most
# recent ones): enough for the members of a small n_gs manifold on a few
# sites. Each is a few matrices of the size of the Krylov space.
SPECTRUM_CACHE_SIZE = 8


def _spectrum_key(Ms,Ws,ops,settings):
    """A digest of everything task.spectrum_from_arrays reads: the site
    arrays of the state, of H and of every operator, and the settings. The
    key is the content and not a flag of the chain, so a state set with
    set_gs, a new Hamiltonian or a re-solve that moved the state each give
    a new key, and nothing has to be invalidated by hand."""
    import hashlib
    from .pyitensor import backend
    h = hashlib.blake2b(digest_size=20)
    h.update(repr((settings,backend.backend_name())).encode())
    def feed(arrays):
        h.update(b"|%d" % len(arrays))
        for a in arrays:
            a = np.ascontiguousarray(backend.to_host(a))
            h.update(repr((a.shape,a.dtype.str)).encode())
            h.update(a.data)
    feed(Ms)
    feed(Ws)
    for W,Wd in ops:
        feed(W)
        feed(Wd)
    return h.hexdigest()


def session_spectrum(self,ops_terms,nkry,pad_bonds=None,expand_bonds=None,
        expand_rounds=1,residual=False,column=False,gauge=False):
    """(e, mats, dim) of task.projected_spectrum on the chain's session,
    ops_terms a list of (terms of O, terms of O^dagger). pad_bonds and
    expand_bonds as in session_poles. With residual=True a fourth entry
    is the residual norm of every state (task.projected_spectrum), zero
    for the ones the Krylov space has converged. column=True builds row
    and column 0 of every matrix only, the tangent-tangent block being
    zero (task.projected_spectrum); it is kept under its own key, so a
    full build is never served one of these.

    On "python" the state is first rotated to a real gauge when H is real
    (real_gauge_state), so that the spectrum runs in real arithmetic; the
    part of the state this drops is at most REAL_GAUGE_DROP of its norm.
    With gauge=True a last entry gives it: a float (0.0 for a state that
    is real as exported, on every backend), or None when the complex gauge
    was kept.

    The result is kept on the chain under a digest of the arrays it was
    computed from (_spectrum_key), so a second call on the same state, the
    same Hamiltonian and the same operators costs the export and the digest
    only. This is what makes a scan of get_kondo_spectrum over its coupling
    constants or its bias grid cost one Krylov space: the spectrum does not
    depend on them. Copies are returned, the stored arrays are never handed
    out."""
    s = self._session
    from .pyitensor import task
    task.check_expansion_keywords(pad_bonds,expand_bonds,expand_rounds)
    if self.itensor_version=="julia_live":
        Ms,_ = _julia_state(self)
        Ws = _julia_mpo_arrays(self,self.hamiltonian)
        ops = [(_julia_mpo_arrays(self,t),_julia_mpo_arrays(self,td))
               for t,td in ops_terms]
        dropped = _exported_real(Ms)
    elif self.itensor_version=="python":
        Ms,_ = _python_state(s)
        Ws = task.mpo_arrays(s.H)
        ops = [(task.mpo_arrays(s._mpo(t)),task.mpo_arrays(s._mpo(td)))
               for t,td in ops_terms]
        Ms,Ws,dropped = real_gauge_state(Ms,Ws) # the digest is of these
    else:
        Ms = s.mps_site_arrays(s.gs_wavefunction())
        Ws = s.mpo_site_arrays(s.hamiltonian_mpo())
        ops = [(_v3_arrays(s,t),_v3_arrays(s,td)) for t,td in ops_terms]
        dropped = _exported_real(Ms)
    key = _spectrum_key(Ms,Ws,ops,(int(nkry),pad_bonds,expand_bonds,
            int(expand_rounds),bool(column)))
    cache = self.__dict__.setdefault("_task_spectrum_cache",{})
    if key in cache:
        cache[key] = cache.pop(key) # most recent last
    else:
        cache[key] = task.spectrum_from_arrays(Ms,Ws,ops,nkry,
                pad_bonds=pad_bonds,expand_bonds=expand_bonds,
                expand_rounds=expand_rounds,residual=True,
                column=bool(column))
        while len(cache)>SPECTRUM_CACHE_SIZE:
            cache.pop(next(iter(cache)))
    e,mats,dim,res = cache[key]
    out = (e.copy(),[m.copy() for m in mats],dim)
    if residual: out = out+(res.copy(),)
    if gauge: out = out+(dropped,)
    return out


def get_poles(self,name=None,i=0,j=0,nkry=100,variance=False,pad_bonds=None,
        expand_bonds=None,expand_rounds=1):
    """The discrete poles of <GS|A delta(w-H+E0) B|GS> in the tangent space.

    Returns a dict with "omega" (poles relative to E0), "weight" (complex,
    <GS|A|alpha><alpha|B|GS>), "elastic" (the n=GS weight <A><B>, at w=0),
    "variance" (the two-site variance of each pole, or None),
    "dimension" (of the tangent space) and "missing_weight", the sum rule
    <GS|A B|GS> minus the elastic weight minus the pole weights.

    The variance measures how well each pole is an eigenstate; it does
    not see the part of B|GS> outside the tangent space, which is dropped
    from the start. That part is zero for a sum of one-site operators,
    and for a fermion operator on a state of definite parity, but not for
    a bond operator or for C/Cdag on a parity mixture (a degenerate
    Majorana pair), where a truncated maxm lost up to all of the weight.
    missing_weight measures it, and a RuntimeWarning is raised when it is
    above 1e-6 of the weight.

    pad_bonds enlarges the tangent space by padding every bond of the
    state with zero-weight directions before the recursion, up to
    min(pad_bonds, structural bound), or to the structural bound for
    "full" (the exact Lanczos of ED on a short chain); None (the default)
    leaves the state as it is. It is the cure for a state of small
    Schmidt rank, a product state or a state solved in one sector of a
    degenerate pair, whose tangent space misses part of B|GS>; the
    padded directions below the bound are arbitrary, see
    pyitensor/task.py's module docstring. With "full" the variances are
    zero; with an int below the bound they are inflated by the padded
    directions and are no longer a sharp error bar of the poles. A padded
    tangent space above task.MAX_PADDED_DIMENSION raises a ValueError.
    The keyword is unrelated to pyitensor's backend.set_pad_bonds, which
    freezes bond shapes for JAX.

    expand_bonds (None, the default, or a positive int k) is the
    controlled version: every bond grows by at most k directions taken
    from the two-site discarded space of B|GS> and A^dagger|GS> and of
    the first Lanczos step H P B|GS>, the largest singular directions of
    what the recursion would otherwise drop. The state does not change.
    Fewer than k are added where that discarded weight is below 1e-8 of
    its largest value over the chain. expand_rounds (default 1) repeats
    the expansion from the expanded state, so that b and H b see the
    directions of the round before; it stops early when a round adds
    nothing. It is what a product state needs: at J = 0 the first round
    brings back the whole weight of the composite vertex, and with k of 4
    or more the second one the poles of ED. pad_bonds and expand_bonds cannot be combined, and an expanded tangent space above
    task.MAX_PADDED_DIMENSION raises a ValueError."""
    _check_backend(self)
    from .pyitensor.task import check_expansion_keywords
    check_expansion_keywords(pad_bonds,expand_bonds,expand_rounds)
    from .multioperatortk import canonical
    A,B = operatornames.str2MO(self,name,i=i,j=j,
            require_symbolic_for="submode=\"TASK\"")
    if (isinstance(nkry,bool) or not isinstance(nkry,(int,np.integer))
            or nkry<1):
        raise ValueError("nkry must be a positive integer, got %r"%(nkry,))
    ground_state(self)
    if self.itensor_version=="julia_live": # MultiOperators, see _julia_mpo
        terms = lambda m: m
    else: terms = lambda m: m.to_terms()
    Adag = None if canonical.is_dagger_pair(A,B) else terms(A.get_dagger())
    ws,wt,el,var,dim,ab = session_poles(self,terms(A),terms(B),Adag,
            int(nkry),variance=variance,terms_ij=terms(A*B),
            pad_bonds=pad_bonds,expand_bonds=expand_bonds,
            expand_rounds=expand_rounds)
    missing = ab-el-complex(np.sum(wt))
    scale = max(abs(ab),abs(el)+float(np.sum(np.abs(wt))))
    if abs(missing)>1e-6*scale:
        import warnings
        warnings.warn(
            "submode=\"TASK\": the poles and the elastic term miss %.3g of "
            "the sum rule <A B> = %.6g. The part of B|GS> outside the "
            "one-site tangent space is dropped, which the pole variances "
            "do not show: a bond operator, a fermion operator on a state "
            "of mixed parity (a degenerate Majorana pair), or, for A not "
            "B^dagger, a Krylov space too small for A^dagger|GS>, or a "
            "state of small Schmidt rank (a product state). Raise maxm (or "
            "nkry), project the state onto one parity sector, pass "
            "expand_bonds (with expand_rounds=2 for a product state) to "
            "get_poles, dynamical_correlator, ysr.get_stm_spectrum or "
            "zaverage.stm_poles, or "
            "pad_bonds=\"full\" on a short chain, which makes the tangent "
            "space the whole Hilbert space."
            %(abs(missing),ab.real if ab.imag==0 else abs(ab)),
            RuntimeWarning,stacklevel=3)
    return dict(omega=ws,weight=wt,elastic=el,variance=var,dimension=dim,
            missing_weight=missing)


def dynamical_correlator(self,es=np.linspace(-1.,10,100),delta=1e-1,
        name=None,i=0,j=0,nkry=100,return_poles=False,variance=None,
        pad_bonds=None,expand_bonds=None,expand_rounds=1):
    """Dynamical correlator from the tangent-space Krylov poles, each one
    broadened into a Lorentzian of half width delta.

    With return_poles=True returns (es, correlator, info), info being
    get_poles()'s dict; the two-site variance of every pole is computed
    then unless variance=False, since it is the error bar of the poles
    (not with an int pad_bonds below the bound, where the padded
    directions inflate it). pad_bonds, expand_bonds and expand_rounds are
    get_poles()'s."""
    if delta<=0.0: raise ValueError("delta must be positive, got %s"%delta)
    if variance is None: variance = return_poles
    if es is None: raise ValueError("es must be an array of frequencies")
    es = np.atleast_1d(np.asarray(es,dtype=float))
    if es.ndim!=1: raise ValueError("es must be a one-dimensional array "
            "of frequencies, got shape %s"%(es.shape,))
    p = get_poles(self,name=name,i=i,j=j,nkry=nkry,variance=variance,
            pad_bonds=pad_bonds,expand_bonds=expand_bonds,
            expand_rounds=expand_rounds)
    ws = np.concatenate([p["omega"],[0.]])
    wt = np.concatenate([p["weight"],[p["elastic"]]])
    lor = delta/(np.pi*((es[:,None]-ws[None,:])**2+delta**2))
    if return_poles: return es,lor@wt,p
    return es,lor@wt
