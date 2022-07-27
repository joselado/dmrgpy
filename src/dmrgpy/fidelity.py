import numpy as np
from .algebra import algebra

def get_fidelity(MBO,H0,H1,l=1.0,delta=1e-4,
      dl = 1e-2,n=2,
      fmode="derivative",**kwargs):
    """Given two Hamiltonians H0 and H1, compute the fidelity as a function
    of lambda.
    -fmode: mode of the calculation"""
    MBO = MBO.copy() # copy the object
    if fmode=="PT": # perturbation theory mode
        H = H0 + l*H1 # full Hamiltonian
        MBO.set_hamiltonian(H) # set the Hamiltonian
        (es,wfs) = MBO.get_excited_states(purify=True,**kwargs) # get excited states
        print("Excited states",np.round(es,2))
        chi = 0. # initialize
        e0,wf0 = es[0],wfs[0] # ground state
        for n in range(1,len(es)): # loop over states
            den = np.abs(es[n]-e0)**2 + delta**2
            wij = np.abs(wf0.dot(H1*wfs[n]))
            print(np.round(wij,2))
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
        (e,wfsp) = MBO.get_excited_states(n=ngs,**kwargs) # get GS manifold
        MBO.set_hamiltonian(Hm) # set the Hamiltonian
        (e,wfsm) = MBO.get_excited_states(n=ngs,**kwargs) # get GS manifold
        if ngs>1: # not properly implemented, this is just a workaround
            from .algebra.algebra import smooth_gauge
            wfsp = smooth_gauge(wfs,wfsp) # smooth gauge
            wfsm = smooth_gauge(wfs,wfsm) # smooth gauge
            def Uij(ws1,ws2):
                m = np.zeros((ngs,ngs),dtype=np.complex) # empy matrix
                for i in range(ngs):
                  for j in range(ngs):
                      m[i,j] = ws1[i].dot(ws2[j]) # scalar product
                return m # return matrix
            Um = Uij(wfs,wfsm) # below
            Up = Uij(wfs,wfsp) # above
            U = Uij(wfs,wfs) # above
            Ud = ((Up-U)-(U-Um))/dl**2 # difference
            chi = np.abs(algebra.det(Ud))
            return chi
        else: # just one wave
            w = wfs[0]
            wp = wfsp[0]
            wm = wfsm[0]
            # factor the phase just in case
            wp = wp*np.exp(1j*np.angle(wp.dot(w)))
            wm = wm*np.exp(1j*np.angle(wm.dot(w)))
            dw2 = (wp - w) - (w - wm) # difference
            chi = np.abs(w.dot(dw2))/(dl**2)
            return chi





