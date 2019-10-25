from .manybodychain import Many_Body_Hamiltonian
import numpy as np
from .dmrgpy2pychain import correlator as correlatorpychain
from .algebra import algebra
from . import effectivehamiltonian
from . import pychainwrapper

class Coupling():
  def __init__(self,i,j,g):
    self.i = i
    self.j = j
    self.g = g

Spin_Hamiltonian = Many_Body_Hamiltonian

class Spin_Hamiltonian(Many_Body_Hamiltonian):
    """Class for spin Hamiltonians"""
    def __init__(self,sites):
        Many_Body_Hamiltonian.__init__(self,sites)
        # default exchange constants
        self.set_exchange(lambda i,j: abs(i-j)==1*1.0)
    def set_exchange(self,fun):
      """Set the exchange coupling between sites"""
      one = np.matrix(np.identity(3))
      self.computed_gs = False # say that GS has not been computed
      self.exchange = [] # empty list
      for i in range(self.ns): # loop
        for j in range(self.ns):  # loop
          g = fun(i,j).real # call the function
          if np.sum(np.abs(fun(i,j)-fun(j,i)))>1e-5: raise # something wrong
          g = g*one # multiply by the identity
          if np.sum(np.abs(g))!=0.0:
            g = g/2. # normalize
            c = Coupling(i,j,g) # create class
            self.exchange.append(c) # store
    def gs_energy(self,mode="DMRG",**kwargs):
        """
        Return the ground state energy
        """
        if mode=="DMRG":
          return Many_Body_Hamiltonian.gs_energy(self,**kwargs)
        elif mode=="ED":
          return pychainwrapper.gs_energy(self,**kwargs)
    def get_pychain(self):
        return pychainwrapper.get_pychain(self)
    def get_full_hamiltonian(self):
        """Return the full Hamiltonian"""
        from . import pychainwrapper
        return pychainwrapper.get_full_hamiltonian(self)
    def sisj_edge(self):
        if not self.computed_gs: self.get_gs() # compute ground state
        pairs = [(0,i) for i in range(self.ns)] # create pairs
        cs = self.correlator(pairs=pairs,mode="DMRG") # compute
        ns = range(len(cs)) # number of correlators
        return np.array([ns,cs])
    def get_magnetization(self,mode="DMRG",**kwargs):
        if mode=="DMRG": # using DMRG
            pairs = [(i,i) for i in range(self.ns)]
            mx = self.get_correlator(pairs=pairs,name="X").real
            my = self.get_correlator(pairs=pairs,name="Y").real
            mz = self.get_correlator(pairs=pairs,name="Z").real
            np.savetxt("MAGNETIZATION.OUT",np.matrix([mx,my,mz]).T)
            return (mx,my,mz)
        elif mode=="ED": # using ED
            return pychainwrapper.get_magnetization(self,**kwargs)
        else: raise
    def get_dynamical_correlator(self,mode="DMRG",**kwargs):
        """
        Compute the dynamical correlator
        """
        if mode=="DMRG": # Use MPS
            return Many_Body_Hamiltonian.get_dynamical_correlator_MB(self,
                    **kwargs)
        else: # use exact diagonalization
            from . import pychainwrapper
            return pychainwrapper.get_dynamical_correlator(self,
                    **kwargs)
    def get_excited(self,mode="DMRG",n=10,**kwargs):
        """
        Compute excited state energies
        """
        if mode=="DMRG":
            return Many_Body_Hamiltonian.get_excited(self,n=n,
                    **kwargs)
        elif mode=="ED":
            h = self.get_full_hamiltonian() # get the Hamiltonian
            from . import pychain
            return pychain.spectrum.eigenstates(h,k=n) # return energies
        else: raise
    def get_dos(self,mode="DMRG",**kwargs):
        """
        Compute the overall density of states
        """
        if mode=="DMRG":
            return Many_Body_Hamiltonian.get_dos(self,**kwargs)
        elif mode=="ED":
            h = self.get_full_hamiltonian() # get the Hamiltonian
            from .pychain import dos
            return dos.dos_kpm(h,**kwargs)

    def get_correlator(self,pairs=[[]],mode="DMRG",**kwargs):
        """Return the correlator"""
        if mode=="DMRG": # using DMRG
            return Many_Body_Hamiltonian.get_correlator(self,pairs=pairs,
                    **kwargs)
        elif mode=="ED": # using exact diagonalization
            self.to_folder()
            m = correlatorpychain.correlator(self,pairs=pairs,**kwargs)
            self.to_origin() # go to main folder
            return m
        else: raise
    def get_effective_hamiltonian(self,**kwargs):
        """Return the effective Hamiltonian"""
        return effectivehamiltonian.get_effective_hamiltonian(self,
                    name="XX",**kwargs)
