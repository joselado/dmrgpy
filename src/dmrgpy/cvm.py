import numpy as np
from . import operatornames
from . import multioperator

def dynamical_correlator(self,es=np.linspace(0.,10.0,100),
        delta=1e-1,name="XX",i=0,j=0):
    """
    Compute the dynamical correlator using CVM method in DMRG
    """
    if not self.computed_gs: self.get_gs() # compute ground state
    out = [] # empty list
    for e in es: # loop over energies
        print("CVM in E = ",e)
        o = cvm_dmrg(self,name=name,i=i,j=j,delta=delta,e=e)
        out.append(o) # store
    out = np.array(out)
#    from .inference import points2function
#    (es,out) = points2function(es,out)
    return (es,out) # return result










def cvm_dmrg(self,name="XX",delta=1e-1,e=0.0,**kwargs):
    """
    Return the dynamical correlator for a single energy, via the
    in-process pybind11 extension (mpscpp2/chain_session.h's
    Chain::cvm_dynamical_correlator). Parameter names differ from the C++
    side's (omega, eta, energy) to match this file's existing naming:
    e -> omega, delta -> eta, self.e0 -> energy.
    """
    name = operatornames.str2MO(self,name,**kwargs)
    A = name[0]
    B = name[1]
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    return self._session.cvm_dynamical_correlator(
            A.to_terms(),B.to_terms(),e,delta,self.e0,
            self.cvm_tol,int(self.cvm_nit))




def dynamical_correlator_analytic_continuation(self,name=None,
        delta=1e-1,es=np.linspace(0.,5.0,300)):
    """
    Compute the dynamical correlator using analytic continuation
    """
    A,B = name[0],name[1]
    wf = self.get_gs() # get the ground state
    wfa = A.get_dagger()*wf # apply A to the GS
    wfb = B*wf # apply B to the GS
    e0 = self.gs_energy() # ground state energy
    Hp = self.hamiltonian - e0
    def f(e): # function to compute
        wfi = self.applyinverse(-self.hamiltonian+(e0+e),wfa)
        return wfb.dot(wfi) # return result
#    return es,-np.array([f(e+1j*delta*10) for e in es]).imag*2/np.pi # brute force
    from .analyticcontinuation import imag2real
    xz = es*1j
    xz = np.linspace(delta*10,10.,100)*1j
    xz = np.concatenate([-xz,xz])
    xz = np.linspace(min(es),max(es),20) + delta*40*1j
#    xz = [np.random.random()-.5+1j*np.random.random()+0.5j for i in range(40)]
#    xz = 40.*np.array(xz)
    outz = np.array([f(z) for z in xz]) # complex axis
    esz,out = imag2real(xz,outz,x=es+1j*delta)
    out = -out.imag*2/np.pi
    return es,out



from .nonhermitian.dynamics import dynamical_correlator_cvm_explicit

#def dynamical_correlator_cvm_explicit(self,name=None,
#        delta=1e-1,es=np.linspace(0.,5.0,300)):
#    """
#    Compute the dynamical correlator using analytic continuation
#    """
#    ### So far this just works for onsite correlators
#    A,B = name[0],name[1]
#    if not self.is_zero_operator(A.get_dagger()-B): 
#        print("Only implemented for A^\dagger=B")
#        raise
#    wf = self.get_gs() # get the ground state
#    wfa = A.get_dagger()*wf # apply A to the GS
#    wfb = B*wf # apply B to the GS
#    e0 = self.gs_energy() # ground state energy
#    Hp = self.hamiltonian - e0
#    def f(e,delta): # function to compute
#        wfi = self.applyinverse(-self.hamiltonian+(e0+e+1j*delta),wfa)
#        return wfb.dot(wfi) # return result
#    from .analyticcontinuation import imag2real
#    outz = np.array([f(z,delta) - f(z,-delta) for z in es]) # complex axis
#    return es,1j*outz/np.pi





