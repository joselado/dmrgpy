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

# dictionary for the sites, with a more readable nomenclature
label2site = dict() # dictionary
label2site["1/2"] = 2
label2site["S=1/2"] = 2
label2site[2] = 2
label2site["1"] = 3
label2site["S=1"] = 3
label2site[3] = 3
label2site["3/2"] = 4
label2site["S=3/2"] = 4
label2site[4] = 4
label2site["2"] = 5
label2site["S=2"] = 5
label2site[5] = 5
label2site["5/2"] = 6
label2site["S=5/2"] = 6
label2site["S=3"] = 7
label2site[6] = 6


def get_logdimension(self):
    """Return the logarithm of the dimension"""
    return np.sum(np.log(np.array(self.sites))) # return dimension



def get_site(label):
    if label in label2site: return label2site[label]
    else: return None


class Spin_Chain(Many_Body_Chain):
    """Class for spin Hamiltonians"""
    def __init__(self,sites,**kwargs):
        sites = [label2site[s] for s in sites]
        Many_Body_Chain.__init__(self,sites,**kwargs)
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
    def test(self,ntries=3,**kwargs):
        """Check the anticommunation relations"""
        Sx = self.Sx
        Sy = self.Sy
        Sz = self.Sz
        for ii in range(ntries):
            i = np.random.randint(self.ns)
            j = np.random.randint(self.ns)
            op = Sx[i]*Sy[j] - Sy[j]*Sx[i]
            if i==j: op = op - 1j*Sz[i]
            if not self.is_zero_operator(op,**kwargs): raise
    def get_logdimension(self):
        return get_logdimension(self)
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
    def get_ED_obj(self):
        if self.has_ED_obj: 
            return self.ED_obj
        else:
            self.ED_obj = pychainwrapper.get_pychain(self)
            return self.ED_obj
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
    def get_full_SS_correlator(self,**kwargs):
        """Return the full spin correlator"""
        from .dynamicstk import spincorrelators
        return spincorrelators.get_full_SS_correlator(self,**kwargs)
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
