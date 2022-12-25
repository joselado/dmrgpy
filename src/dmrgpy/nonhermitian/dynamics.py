

import numpy as np





def dynamical_correlator_cvm_explicit(self,name=None,
        delta=1e-1,es=np.linspace(0.,5.0,300)):
    """
    Compute the dynamical correlator using analytic continuation
    """
    print("Non Hermitian mode")
    ### So far this just works for onsite correlators
    A,B = name[0],name[1]
    if not self.is_zero_operator(A.get_dagger()-B):
        print("Only implemented for A^\dagger=B")
        raise
    wf = self.get_gs() # get the ground state
    wfa = A.get_dagger()*wf # apply A to the GS
    wfb = B*wf # apply B to the GS
    e0 = self.gs_energy() # ground state energy
    Hp = self.hamiltonian - e0
    def f(e,delta): # function to compute
        wfi = self.applyinverse(-self.hamiltonian+(e0+e+1j*delta),wfa)
        return wfb.dot(wfi) # return result
#    from .analyticcontinuation import imag2real
    outz = np.array([f(z,delta) - f(z,-delta) for z in es]) # complex axis
    return es,np.abs(1j*outz/np.pi)




dynamical_correlator_non_hermitian = dynamical_correlator_cvm_explicit






