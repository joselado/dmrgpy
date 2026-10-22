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


def session_poles(self,terms_i,terms_j,terms_idag,nkry,variance=False):
    """(w, weights, elastic, variances, dim) on the chain's session, which
    must hold the ground state (groundstate.ground_state_on_session)."""
    s = self._session
    if self.itensor_version=="python":
        return s.task_dynamical_correlator(terms_i,terms_j,terms_idag,nkry,
                variance=variance)
    from .pyitensor import task
    wf = s.gs_wavefunction()
    Wadag = None if terms_idag is None else _v3_arrays(s,terms_idag)
    ws,wt,var,dim = task.correlator_poles(s.mps_site_arrays(wf),
            s.mpo_site_arrays(s.hamiltonian_mpo()),_v3_arrays(s,terms_j),
            Wadag,nkry,variance)
    elastic = s.vev(terms_i,wf)*s.vev(terms_j,wf)
    return ws,wt,complex(elastic),var,dim


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
    "variance" (the two-site variance of each pole, or None) and
    "dimension" (of the tangent space)."""
    _check_backend(self)
    from .multioperatortk import canonical
    from .groundstate import ground_state_on_session
    A,B = operatornames.str2MO(self,name,i=i,j=j,
            require_symbolic_for="submode=\"TASK\"")
    if nkry<1: raise ValueError("nkry must be a positive integer, got %r"%nkry)
    ground_state_on_session(self)
    Adag = None if canonical.is_dagger_pair(A,B) else A.get_dagger().to_terms()
    ws,wt,el,var,dim = session_poles(self,A.to_terms(),B.to_terms(),Adag,
            int(nkry),variance=variance)
    return dict(omega=ws,weight=wt,elastic=el,variance=var,dimension=dim)


def dynamical_correlator(self,es=np.linspace(-1.,10,100),delta=1e-1,
        name=None,i=0,j=0,nkry=100,return_poles=False,variance=None):
    """Dynamical correlator from the tangent-space Krylov poles, each one
    broadened into a Lorentzian of half width delta.

    With return_poles=True returns (es, correlator, info), info being
    get_poles()'s dict; the two-site variance of every pole is computed
    then unless variance=False, since it is the error bar of the poles."""
    if delta<=0.0: raise ValueError("delta must be positive, got %s"%delta)
    if variance is None: variance = return_poles
    p = get_poles(self,name=name,i=i,j=j,nkry=nkry,variance=variance)
    es = np.asarray(es)
    ws = np.concatenate([p["omega"],[0.]])
    wt = np.concatenate([p["weight"],[p["elastic"]]])
    lor = delta/(np.pi*((es[:,None]-ws[None,:])**2+delta**2))
    if return_poles: return es,lor@wt,p
    return es,lor@wt
