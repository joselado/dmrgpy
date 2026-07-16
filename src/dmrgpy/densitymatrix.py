# routines to compute density matrices
import numpy as np

def reduced_dm(self,i=0):
    """
    Compute the reduced density matrix
    """
    wf = self.get_gs() # compute ground state
    return self._session.reduced_dm(wf.cpp_handle,i+1)



def reduced_dm_projective(self,wf,i=0,j=None):
    """Compute the reduced density matrix using a brute force approach"""
    from .fermionchain import Fermionic_Chain
    from .fermionchain import Spinful_Fermionic_Chain
    from .spinchain import Spin_Chain
    def projectors(k):
        """Build the projectors onto the different single-site components"""
        site = self.sites[k] # take this site
        if type(self)==Spin_Chain: # spin chain object
          if site==2: # S=1/2 site
              Szk = self.Sz[k] # Sz operator
              P01 = self.Sx[k]+1j*self.Sy[k] # project and rotate
              # we need to project on up/dn and rotate to up!
              return [Szk+0.5,P01] 
          elif site==3: # S=1 site
              raise # not finished
              Szk = self.Sz[k]
              return [Szk*(Szk+1)/2.,-(Szk-1)*(Szk+1),Szk*(1-Szk)/2.] 
          else: raise # not implemented
        elif type(self)==Fermionic_Chain: # spin chain object
          N = self.N[k] # density operator
          P01 = self.Cdag[k] # create an electron
          # we need to project on up/dn and rotate to up!
          return [N,P01] # return the projectors
        elif type(self)==Fermionic_Chain: # spin chain object
          raise # not implemented yet
        else: raise # not implemented
    Pi = projectors(i) # projectors for site i
    if j is not None: Pj = projectors(j) # projectors for site j
    else: Pj = [1] # workaround for a single site
    Pk = [] # projectors in the ij subspace
    for pi in Pi:
        for pj in Pj: Pk.append(pi*pj) # store the projector in this subspace
    # now compute the density matrix in this subspace
    n = len(Pk) # number of components
    dm = np.zeros((n,n),dtype=np.complex128) # initialize
    for i in range(n):
        wfi = Pk[i]*wf # projector
        for j in range(n):
            dm[i,j] = wfi.aMb(Pk[j],wf) # project
#    print(dm.real)
    return dm # return density matrix



#
#def explicit_dm(sc,wf,inds=[0]):
#    """Compute the density matrix explicitly by summing
#    over all the vectors. This is a very heavy procedure, but good
#    for benchmarking and debugging"""
#    sc = sc.copy() # make a copy
#    for site in sc.sites: # check that you only have S=1/2
#        if site !=2: 
#            print("Only implemented for S=1/2")
#            raise # stop
#    # now loop over all the sites
#    for ii in inds: # loop over sites where you want the entropy
#        for Bz in [-1,1]: # the two combinations
#            Hi = sc.Sz[ii] # the two magnetic fields
#
#
#

