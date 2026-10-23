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


def session_poles(self,terms_i,terms_j,terms_idag,nkry,variance=False,
        terms_ij=None):
    """(w, weights, elastic, variances, dim) on the chain's session, which
    must hold the ground state (groundstate.ground_state_on_session).
    With terms_ij, the terms of the product A*B, a sixth entry is
    <GS|A B|GS> on the same state, the sum rule of the weights."""
    s = self._session
    if self.itensor_version=="python":
        out = s.task_dynamical_correlator(terms_i,terms_j,terms_idag,nkry,
                variance=variance)
        if terms_ij is None: return out
        return out+(complex(s.vev(terms_ij,s.wf0)),)
    from .pyitensor import task
    wf = s.gs_wavefunction()
    Wadag = None if terms_idag is None else _v3_arrays(s,terms_idag)
    ws,wt,var,dim = task.correlator_poles(s.mps_site_arrays(wf),
            s.mpo_site_arrays(s.hamiltonian_mpo()),_v3_arrays(s,terms_j),
            Wadag,nkry,variance)
    elastic = s.vev(terms_i,wf)*s.vev(terms_j,wf)
    if terms_ij is None: return ws,wt,complex(elastic),var,dim
    return ws,wt,complex(elastic),var,dim,complex(s.vev(terms_ij,wf))


def session_spectrum(self,ops_terms,nkry):
    """(e, mats, dim) of task.projected_spectrum on the chain's session,
    ops_terms a list of (terms of O, terms of O^dagger)."""
    s = self._session
    if self.itensor_version=="python":
        return s.task_projected_spectrum(ops_terms,nkry)
    from .pyitensor import task
    ops = [(_v3_arrays(s,t),_v3_arrays(s,td)) for t,td in ops_terms]
    return task.spectrum_from_arrays(s.mps_site_arrays(s.gs_wavefunction()),
            s.mpo_site_arrays(s.hamiltonian_mpo()),ops,nkry)


def get_poles(self,name=None,i=0,j=0,nkry=100,variance=False):
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
    above 1e-6 of the weight."""
    _check_backend(self)
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
            int(nkry),variance=variance,terms_ij=(A*B).to_terms())
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
            "B^dagger, a Krylov space too small for A^dagger|GS>. Raise "
            "maxm (or nkry), or project the state onto one parity sector."
            %(abs(missing),ab.real if ab.imag==0 else abs(ab)),
            RuntimeWarning,stacklevel=3)
    return dict(omega=ws,weight=wt,elastic=el,variance=var,dimension=dim,
            missing_weight=missing)


def dynamical_correlator(self,es=np.linspace(-1.,10,100),delta=1e-1,
        name=None,i=0,j=0,nkry=100,return_poles=False,variance=None):
    """Dynamical correlator from the tangent-space Krylov poles, each one
    broadened into a Lorentzian of half width delta.

    With return_poles=True returns (es, correlator, info), info being
    get_poles()'s dict; the two-site variance of every pole is computed
    then unless variance=False, since it is the error bar of the poles."""
    if delta<=0.0: raise ValueError("delta must be positive, got %s"%delta)
    if variance is None: variance = return_poles
    if es is None: raise ValueError("es must be an array of frequencies")
    es = np.atleast_1d(np.asarray(es,dtype=float))
    if es.ndim!=1: raise ValueError("es must be a one-dimensional array "
            "of frequencies, got shape %s"%(es.shape,))
    p = get_poles(self,name=name,i=i,j=j,nkry=nkry,variance=variance)
    ws = np.concatenate([p["omega"],[0.]])
    wt = np.concatenate([p["weight"],[p["elastic"]]])
    lor = delta/(np.pi*((es[:,None]-ws[None,:])**2+delta**2))
    if return_poles: return es,lor@wt,p
    return es,lor@wt
