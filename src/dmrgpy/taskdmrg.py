"""submode="TASK": the tangent-space Krylov dynamical correlator
(Kovalska, von Delft and Gleis, arXiv:2510.07279).

The Lanczos recursion runs in the one-site tangent space of the ground-state
MPS rather than on compressed MPS, so the Krylov basis stays orthogonal and
the result is a set of discrete poles w_alpha with weights M_alpha, broadened
here with the house Lorentzian (dynamics.py's module docstring). The method
itself is in pyitensor/task.py, together with its one approximation, the
projection of H onto the tangent space, and the two-site variance that
measures it pole by pole. On itensor_version="python" and 3; v3 hands its
ground-state MPS and MPO site tensors over as numpy arrays and the same
Python code runs on them.
"""
import numpy as np
from . import operatornames


BACKENDS = ("python",3)

def _check_backend(self,what="get_dynamical_correlator"):
    if self.itensor_version not in BACKENDS:
        raise NotImplementedError(
            "%s: submode=\"TASK\" is implemented on itensor_version="
            "\"python\" and 3 only (got itensor_version=%r)"
            %(what,self.itensor_version))


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
    (chain._mps_arrays_lpr), on 3 those of the session's own export."""
    s = self._session
    from .pyitensor import task
    task.check_expansion_keywords(pad_bonds,expand_bonds,expand_rounds)
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


def session_spectrum(self,ops_terms,nkry,pad_bonds=None,expand_bonds=None,
        expand_rounds=1):
    """(e, mats, dim) of task.projected_spectrum on the chain's session,
    ops_terms a list of (terms of O, terms of O^dagger). pad_bonds and
    expand_bonds as in session_poles."""
    s = self._session
    from .pyitensor import task
    task.check_expansion_keywords(pad_bonds,expand_bonds,expand_rounds)
    if self.itensor_version=="python":
        Ms,_ = _python_state(s)
        ops = [(task.mpo_arrays(s._mpo(t)),task.mpo_arrays(s._mpo(td)))
               for t,td in ops_terms]
        return task.spectrum_from_arrays(Ms,task.mpo_arrays(s.H),ops,nkry,
                pad_bonds=pad_bonds,expand_bonds=expand_bonds,
            expand_rounds=expand_rounds)
    ops = [(_v3_arrays(s,t),_v3_arrays(s,td)) for t,td in ops_terms]
    return task.spectrum_from_arrays(s.mps_site_arrays(s.gs_wavefunction()),
            s.mpo_site_arrays(s.hamiltonian_mpo()),ops,nkry,
            pad_bonds=pad_bonds,expand_bonds=expand_bonds,
            expand_rounds=expand_rounds)


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
    from .groundstate import ground_state_on_session
    A,B = operatornames.str2MO(self,name,i=i,j=j,
            require_symbolic_for="submode=\"TASK\"")
    if (isinstance(nkry,bool) or not isinstance(nkry,(int,np.integer))
            or nkry<1):
        raise ValueError("nkry must be a positive integer, got %r"%(nkry,))
    ground_state_on_session(self)
    Adag = None if canonical.is_dagger_pair(A,B) else A.get_dagger().to_terms()
    ws,wt,el,var,dim,ab = session_poles(self,A.to_terms(),B.to_terms(),Adag,
            int(nkry),variance=variance,terms_ij=(A*B).to_terms(),
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
