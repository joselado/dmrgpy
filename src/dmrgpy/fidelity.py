import numpy as np

def get_fidelity(MBO,H0,H1,l=1.0,delta=1e-3,
      dl = 1e-3,n=2,
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
        wfs = MBO.get_gs_manifold(n=n,**kwargs) # get GS manifold
        ngs = len(wfs) # ground state degeneracy
        MBO.set_hamiltonian(Hp) # set the Hamiltonian
        (e,wfsp) = MBO.get_excited_states(n=ngs,**kwargs) # get GS manifold
        MBO.set_hamiltonian(H) # set the Hamiltonian
        (e,wfsm) = MBO.get_excited_states(n=ngs,**kwargs) # get GS manifold
        if len(wfs)>1: # not implemented
            from .algebra.algebra import smooth_gauge
            wfsp = smooth_gauge(wfs,wfsp) # smooth gauge
            wfsm = smooth_gauge(wfs,wfsm) # smooth gauge
            w = wfs[0]
            wp = wfsp[0]
            wm = wfsm[0]
        else: # just one wave
            w = wfs[0]
            wp = wfsp[0]
            wm = wfsm[0]
            # factor the phase just in case
            wp = wp*np.exp(1j*np.angle(wp.dot(w)))
            wm = wm*np.exp(1j*np.angle(wm.dot(w)))
        chi = np.abs(w.dot(wp)-w.dot(wm))/(dl**2)
        return chi





