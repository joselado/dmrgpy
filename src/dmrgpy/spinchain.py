from .manybodychain import Many_Body_Chain
import numpy as np
from .dmrgpy2pychain import correlator as correlatorpychain
from .algebra import algebra
from . import effectivehamiltonian
from . import pychainwrapper
from . import multioperator

class Coupling():
  def __init__(self,i,j,g):
    self.i = i
    self.j = j
    self.g = g

Spin_Chain = Many_Body_Chain

class Spin_Chain(Many_Body_Chain):
    """Class for spin Hamiltonians"""
    def __init__(self,sites):
        Many_Body_Chain.__init__(self,sites)
        # default exchange constants
        self.use_ampo_hamiltonian = True # use ampo
        self.pychain_object = None # pychain object
        self.Sx = [self.get_operator("Sx",i) for i in range(self.ns)]
        self.Sy = [self.get_operator("Sy",i) for i in range(self.ns)]
        self.Sz = [self.get_operator("Sz",i) for i in range(self.ns)]
        self.Si = [self.Sx,self.Sy,self.Sz]
    def SS(self,i,j):
        return self.Sx[i]*self.Sx[j] + self.Sy[i]*self.Sy[j] + self.Sz[i]*self.Sz[j]
    def set_fields(self,fun):
        h = 0
        for i in range(self.ns):
            b = fun(i)
            for j in range(3):  h = h + b[j]*self.Si[j][i]
        self.fields = h
        self.hamiltonian = self.exchange + self.fields # update Hamiltonian
    def set_exchange(self,fun):
      """Set the exchange coupling between sites"""
      h = 0
      for i in range(self.ns): # loop
        for j in range(self.ns):  # loop
          g = fun(i,j).real # call the function
          if np.sum(np.abs(fun(i,j)-fun(j,i)))>1e-5: raise # something wrong
          one = np.identity(3) # identity matrix
          g = g*one # multiply by the identity
          for ii in range(3):
            for jj in range(3):
                h = h + g[ii,jj]*self.Si[ii][i]*self.Si[jj][j]
      self.exchange = h # exchange matrix
      self.hamiltonian = self.exchange + self.fields # update Hamiltonian
#    def vev(self,MO,mode="DMRG",**kwargs):
#        """ Return a vaccum expectation value"""
#        if mode=="DMRG":
#            return self.vev_MB(MO,**kwargs)
#        elif mode=="ED":
#            SC = self.get_pychain() # get the object
#            SC.hamiltonian = self.get_hamiltonian() # store the MO
#            return SC.vev(MO,**kwargs) # return overlap
#    def gs_energy(self,mode="DMRG",**kwargs):
#        """
#        Return the ground state energy
#        """
#        if mode=="DMRG":
#          return Many_Body_Chain.gs_energy(self,**kwargs)
#        elif mode=="ED":
#          return pychainwrapper.gs_energy(self,**kwargs)
    def get_ED_obj(self):
        return pychainwrapper.get_pychain(self)
    def get_pychain(self):
        return pychainwrapper.get_pychain(self)
    def get_full_hamiltonian(self):
        """Return the full Hamiltonian"""
        from . import pychainwrapper
        return pychainwrapper.get_full_hamiltonian(self)
    def get_magnetization(self,**kwargs):
        mx = [self.vev(self.Sx[i],**kwargs) for i in range(self.ns)]
        my = [self.vev(self.Sy[i],**kwargs) for i in range(self.ns)]
        mz = [self.vev(self.Sz[i],**kwargs) for i in range(self.ns)]
        np.savetxt("MAGNETIZATION.OUT",np.array([mx,my,mz]).T)
        return np.array([mx,my,mz]).real
#    def get_excited(self,mode="DMRG",n=10,**kwargs):
#        """
#        Compute excited state energies
#        """
#        if mode=="DMRG":
#            return Many_Body_Chain.get_excited(self,n=n,
#                    **kwargs)
#        elif mode=="ED":
#            h = self.get_full_hamiltonian() # get the Hamiltonian
#            from . import pychain
#            return pychain.spectrum.eigenstates(h,k=n) # return energies
#        else: raise
#    def get_dos(self,mode="DMRG",**kwargs):
#        """
#        Compute the overall density of states
#        """
#        if mode=="DMRG":
#            return Many_Body_Chain.get_dos(self,**kwargs)
#        elif mode=="ED":
#            h = self.get_full_hamiltonian() # get the Hamiltonian
#            from .pychain import dos
#            return dos.dos_kpm(h,**kwargs)

#    def get_correlator(self,pairs=[[]],mode="DMRG",**kwargs):
#        """Return the correlator"""
#        if mode=="DMRG": # using DMRG
#            return Many_Body_Chain.get_correlator(self,pairs=pairs,
#                    **kwargs)
#        elif mode=="ED": # using exact diagonalization
#            self.to_folder()
#            m = correlatorpychain.correlator(self,pairs=pairs,**kwargs)
#            self.to_origin() # go to main folder
#            return m
#        else: raise
    def get_effective_hamiltonian(self,**kwargs):
        """Return the effective Hamiltonian"""
        return effectivehamiltonian.get_effective_hamiltonian(self,
                    name="XX",**kwargs)
    def get_hamiltonian(self):
        """Return Hamiltonian as a multioperator"""
        if self.hamiltonian is not None: return self.hamiltonian
        else: # conventional way
            sxs = [self.get_operator("Sx",i) for i in range(self.ns)]
            sys = [self.get_operator("Sy",i) for i in range(self.ns)]
            szs = [self.get_operator("Sz",i) for i in range(self.ns)]
            ss = [sxs,sys,szs]
        out = multioperator.zero() # initialize
        for c in self.exchange: # exchange coupling
            for i in range(3):
              for j in range(3):
                out = out + c.g[i,j]*ss[i][c.i]*ss[j][c.j]
        out.clean()
        if len(self.fields)>0:
            for i in range(len(self.fields)):
                b = self.fields[i]
                for j in range(3):
                  out = out + b[j]*ss[j][i]
        # still have to add the fields!!
        return out # return multioperator

Spin_Hamiltonian = Spin_Chain # backwards compatibility
