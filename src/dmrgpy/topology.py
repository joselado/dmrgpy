import numpy as np


def berry_phase_old(f,nk=100):
    """Compute the Berry phase, requires a function
    that returns a manybody Hamiltonian"""
    ks = np.linspace(0.,1.,nk,endpoint=False) # kpoints
    f(0).get_gs() # first iteration
    wfs = [] # empty list
    for k in ks:
        sc = f(k) # get object
        sc.sites_from_file = True # use the old sites
        wfs.append(sc.get_gs()) # get wavefunction
    berry = 1.0 + 0.0j # initialize
    for i in range(nk-1): # loop
      fac = wfs[i].dot(wfs[i+1])
      print(np.abs(fac))
#      if np.abs(fac)>0.9: 
      berry *= fac/np.abs(fac)
    fac = wfs[nk-1].dot(wfs[0]) # last one
    berry *= fac # last factor
    berry = np.arctan2(berry.imag,berry.real)/np.pi # berry phase
    return berry



import numpy as np



def berry_phase_matrix(hkgen,nk=20,**kwargs):
    """ Calculates the Berry phase of a matrix"""
    ks = np.linspace(0.,1.,nk,endpoint=False) # kpoints
    wf0 = occupied_states(hkgen,ks[0],**kwargs) # get occupied states, first k-point
    wfold = wf0.copy() # copy
    m = np.matrix(np.identity(len(wf0))) # initialize as the identity matrix
    for ik in range(1,len(ks)): # loop over k-points, except first one
      wf = occupied_states(hkgen,ks[ik],**kwargs)  # get waves
      m = m@uij(wfold,wf)   # get the uij   and multiply
      wfold = wf.copy() # this is the new old
    m = m@uij(wfold,wf0)   # last one
    d = lg.det(m) # calculate determinant
    phi = np.arctan2(d.imag,d.real)
    return phi # return Berry phase

import scipy.sparse.linalg as slg
import scipy.linalg as lg

def occupied_states(hkgen,k,max_waves=None,states=None):
    """ Returns the WF of the occupied states in a 2d hamiltonian"""
    hk = hkgen(k) # get hamiltonian
    if max_waves is None: es,wfs = lg.eigh(hk) # diagonalize all waves
    else:  es,wfs = slg.eigsh(csc_matrix(hk),k=max_waves,which="SA",
                        sigma=0.0,tol=arpack_tol,maxiter=arpack_maxiter)
    wfs = np.conjugate(wfs.transpose()) # wavefunctions
    if states is None: states = range(0,len(wfs)//2) # states
    occwf = []
    ii = 0
    for (ie,iw) in zip(es,wfs):  # loop over states
        if ii in states: # store
          occwf.append(iw)  # add to the list
        ii += 1
    return np.array(occwf)




def uij(wf1,wf2):
  """ Calcultes the matrix product of two sets of input wavefunctions"""
  out =  np.matrix(np.conjugate(wf1))@(np.matrix(wf2).T)  # faster way
  return out

