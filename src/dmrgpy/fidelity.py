import numpy as np
from .algebra import algebra

def get_fidelity(MBO,H0,H1,l=1.0,delta=1e-4,
      dl = 1e-1,n=2,
      fmode="derivative",**kwargs):
    """Given two Hamiltonians H0 and H1, compute the fidelity as a function
    of lambda.
    -fmode: "derivative" (default) takes the second difference in lambda
     of the ground-state overlap, with step dl. When the ground state is
     degenerate within delta (searched among n states) it takes the second
     difference of |det U| of the overlap matrix of the manifold, which is
     the sum of chi over the degenerate states. "PT" sums the
     perturbative formula over the n lowest states, with delta as the
     broadening of the denominator."""
    MBO = MBO.copy() # copy the object
    if fmode=="PT": # perturbation theory mode
        H = H0 + l*H1 # full Hamiltonian
        MBO.set_hamiltonian(H) # set the Hamiltonian
        # n states enter the sum (the ground state and n-1 excited ones);
        # n used to be dropped here, so the sum always ran over the
        # default two states whatever n the caller passed
        (es,wfs) = MBO.get_excited_states(n=n,purify=True,**kwargs) # get excited states
        chi = 0. # initialize
        e0,wf0 = es[0],wfs[0] # ground state
        for k in range(1,len(es)): # loop over states
            den = np.abs(es[k]-e0)**2 + delta**2
            wij = np.abs(wf0.dot(H1*wfs[k]))
            chi = chi + wij**2/den
        return chi
    elif fmode=="derivative": # use explicit derivative
        Hp = H0 + (l+dl)*H1 # full Hamiltonian
        Hm = H0 + (l-dl)*H1 # full Hamiltonian
        H = H0 + l*H1 # full Hamiltonian
        MBO.set_hamiltonian(H) # set the Hamiltonian
        wfs = MBO.get_gs_manifold(n=n,tol=delta,**kwargs) # get GS manifold
        ngs = len(wfs) # ground state degeneracy
        MBO.set_hamiltonian(Hp) # set the Hamiltonian
        (e,wfsp) = MBO.get_excited_states(n=ngs,purify=False,
                                      **kwargs) # get GS manifold
        MBO.set_hamiltonian(Hm) # set the Hamiltonian
        (e,wfsm) = MBO.get_excited_states(n=ngs,purify=False,
                                      **kwargs) # get GS manifold
        if ngs==1: return _chi_single(wfs[0],wfsp[0],wfsm[0],dl)
        return _chi_manifold(wfs,wfsp,wfsm,dl)
    else: raise ValueError("Unknown fmode "+str(fmode))


def _chi_single(w,wp,wm,dl):
    """Fidelity susceptibility of a single ground state w, from the ground
    states wp and wm at lambda+dl and lambda-dl: the second difference of
    the overlap, |<w|wp> + <w|wm> - 2|/dl^2, with the phases of wp and wm
    fixed so that both overlaps are real and positive."""
    w,wp,wm = [x.copy().normalize() for x in (w,wp,wm)]
    wp = wp*np.exp(1j*np.angle(wp.dot(w)))
    wm = wm*np.exp(1j*np.angle(wm.dot(w)))
    dw2 = (wp - w) - (w - wm) # difference
    return np.abs(w.dot(dw2))/(dl**2)


def _chi_manifold(wfs,wfsp,wfsm,dl):
    """Fidelity susceptibility of a ground-state manifold of g states, the
    matrix generalization of _chi_single: the second difference of |det U|,
    with U_ij = <w_i|w'_j> the overlap matrix of the two manifolds,
    |(|det U+| + |det U-| - 2|det U0|)|/dl^2. For g=1 it is exactly
    _chi_single. For g>1 it is the trace of the quantum metric of the
    manifold, sum_{a in GS, m not in GS} |<m|H1|a>|^2/(E_m-E_0)^2, so it
    is a sum over the g states and not an average. The determinant of the
    second difference of the matrix, used before 2026-10-15, scaled as
    chi^g instead."""
    from .algebra.arnolditk import gram_smith
    wfs = gram_smith(wfs) # orthogonalize
    wfsp = gram_smith(wfsp) # orthogonalize
    wfsm = gram_smith(wfsm) # orthogonalize
    ngs = len(wfs)
    def Uij(ws1,ws2):
        m = np.zeros((ngs,ngs),dtype=np.complex128) # empty matrix
        for i in range(ngs):
          for j in range(ngs):
              m[i,j] = ws1[i].dot(ws2[j]) # scalar product
        return m # return matrix
    dp = np.abs(algebra.det(Uij(wfs,wfsp)))
    dm = np.abs(algebra.det(Uij(wfs,wfsm)))
    d0 = np.abs(algebra.det(Uij(wfs,wfs)))
    return np.abs(dp + dm - 2*d0)/dl**2
