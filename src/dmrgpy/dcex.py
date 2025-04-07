import numpy as np
from . import operatornames
from . import taskdmrg

def dynamical_correlator(self,name="XX",i=0,j=0,delta=2e-2,
        nex=20,es=np.linspace(-1.0,10.0,1000),scale=20.0,**kwargs):
    """
    Compute dynamical correlator with excited states with DMRG
    """
#    self.gs_energy()
#    name = operatornames.str2MO(self,name)
#    task = {"dynamical_correlator_excited":"true",
#            "scale_lagrange_excited":str(scale),
#            "nexcited":str(nex),
#            }
#    self.task = task # override tasks
#    name[0] = name[0].get_dagger()
#    self.execute(lambda: name[0].write(name="dc_multioperator_i.in"))
#    self.execute(lambda: name[1].write(name="dc_multioperator_j.in"))
#    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
#    self.execute( lambda : self.run()) # run calculation
#    # now read the data
#    eex = self.get_file("EXCITED.OUT").T[0] # excitation energies
#    eex = eex[1:len(eex)] - eex[0] # substract GS energy
#    cs = self.get_file("EXCITED_OVERLAPS.OUT") # overlaps with GS
#    c1 = cs[:,0] + 1j*cs[:,1] # first overlap
#    c2 = cs[:,2] - 1j*cs[:,3] # second overlap
    # compute the two correlators
    esex,wsex = self.get_excited_states(n=nex,**kwargs)
    A,B = name[0].get_dagger(),name[1] # operators
    wf0 = wsex[0] # ground state
    from .algebra.arnolditk import gram_smith
    wsex = gram_smith(wsex) # orthogonalize
    Aop = [[wfi.dot(A*wfj) for wfj in wsex] for wfi in wsex] # representation
    Bop = [[wfi.dot(B*wfj) for wfj in wsex] for wfi in wsex] # representation
    h = self.hamiltonian # Hamiltonian
    Hop = [[wfi.dot(h*wfj) for wfj in wsex] for wfi in wsex] # representation
    # transform to arrays
    Aop = np.array(Aop)
    Bop = np.array(Bop)
    Hop = np.array(Hop)
    # now get eigenvalues and eigenvectors
    from scipy.linalg import eigh
    (esex,wsex) = eigh(Hop) ; wsex = wsex.T # transpose 
    # from now on we operate with numpy arrays
    wf0 = wsex[0]
    wfa = Aop@wf0 # A times ground state
    wfb = Bop@wf0 # B times ground state
    c1 = [np.conjugate(wfa).dot(wfi) for wfi in wsex] # matrix element
    c2 = [np.conjugate(wfi).dot(wfb) for wfi in wsex] # matrix element
    eex = esex - esex[0] # difference
    es,adv = dcex(eex,c1,c2,es=es,delta=delta) # return correlator
    es,ret = dcex(eex,c2,c1,es=es,delta=-delta) # return correlator
    return es,1j*(adv-ret)/np.pi # return correlator





def dcex(eex,c1,c2,es=np.linspace(-1.0,10.0,300),delta=1e-1):
    """
    Compute a dynamical correlator with excitation energies and
    matrix elements
    """
    if len(eex)!=len(c1): raise # something wrong
    if len(eex)!=len(c2): raise # something wrong
    if es is None: es = np.linspace(-1.0,np.max(eex),2000)
    out = np.zeros(len(es),dtype=np.complex128) # output
    out = [] # empty list
    out = 0.0j # initialize
    for i in range(len(eex)): # loop over excited states
        out = out + c1[i]*c2[i]/(es-eex[i]+1j*delta) # dynamical correlator
    return es,out # return 








