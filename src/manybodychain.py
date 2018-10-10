from __future__ import print_function
import numpy as np
import os
import kpmdmrg
import pychainwrapper
import pychain.spectrum
import mps

dmrgpath = os.environ["DMRGROOT"] # path for the program
one = np.matrix(np.identity(3))



class Coupling():
  def __init__(self,i,j,g):
    self.i = i
    self.j = j
    self.g = g




class Many_Body_Hamiltonian():
  def __init__(self,sites):
    self.sites = sites # list of the sites
    self.path = os.getcwd()+"/.mpsfolder/" # folder of the calculations
    self.inipath = os.getcwd() # original folder
    self.ns = len(sites) # number of sites
    self.couplings = [Coupling(i,i+1,one) for i in range(self.ns-1)] # empty list
    self.fields = [] # empty list
    self.hoppings = dict() # empty dictionary
    self.spinful_hoppings = dict() # empty dictionary
    self.pairing = dict() # empty dictionary
    self.hubbard = dict() # empty dictionary
#    self.couplings.append(Coupling(0,self.ns-1,one)) # closed boundary
    # additional arguments
    self.kpmmaxm = 50 # bond dimension in KPM
    self.maxm = 30 # bond dimension in wavefunctions
    self.nsweeps = 15 # number of sweeps
    self.kpmcutoff = 1e-8 # cutoff in KPM
    self.cutoff = 1e-8 # cutoff in ground state
    self.kpmscale = 10.0 # this is meaningless
    self.restart = False # restart the calculation
    self.gs_from_file = False # start from a random wavefunction
    self.wf0 = None # no initial WF
    self.starting_file_gs = "starting_psi_GS.mps" # initial file for GS
    self.sites_from_file = False # read sites from the file
    self.computed_gs = False # computed the GS already
    os.system("mkdir -p "+self.path) # create folder for the calculations
  def to_folder(self): os.chdir(self.path) # go to calculation folder
  def to_origin(self): 
    if os.path.isfile(self.path+"/ERROR"): raise # something wrong
    os.chdir(self.inipath) # go to original folder
  def clean(self): os.system("rm -rf "+self.path) # clean temporal folder
  def set_exchange(self,fun):
    """Set the exchange coupling between sites"""
    self.computed_gs = False # say that GS has not been computed
    self.couplings = [] # empty list
    for i in range(self.ns): # loop
      for j in range(i+1,self.ns):  # loop
        g = fun(i,j).real # call the function
        g = g*one # multiply by the identity
        if np.sum(np.abs(g))!=0.0: 
          c = Coupling(i,j,g) # create class
          self.couplings.append(c) # store
    # now the onsite ones
    for i in range(self.ns): # loop
        g = fun(i,i).real # call the function
        g = g*one # multiply by the identity
        g = (g + g.H)/2. # the onsite one must be hermitian
        if np.sum(np.abs(g))!=0.0: # if nonzero 
          c = Coupling(i,i,g) # create class
          self.couplings.append(c) # store
  def set_hoppings(self,fun):
      """Add the spin independent hoppings"""
      self.computed_gs = False # say that GS has not been computed
      self.hoppings = dict()
      for i in range(self.ns): # loop
          for j in range(self.ns): # loop
              if self.sites[i]==1 and self.sites[j]==1:
                  c = fun(i,j)
                  if np.abs(c)>0.0:
                      self.hoppings[(i,j)] = Coupling(i,j,c) # store
  def set_spinful_hoppings(self,fun):
      """Add the spin independent hoppings"""
      self.computed_gs = False # say that GS has not been computed
      if callable(fun):
        self.spinful_hoppings = dict()
        for i in range(self.ns): # loop
            for j in range(self.ns): # loop
                if self.sites[i]==1 and self.sites[j]==1:
                    c = fun(i,j)
                    if np.abs(c)>0.0:
                        self.spinful_hoppings[(i,j)] = Coupling(i,j,c) # store
      else: # assume it is a matrix
          self.spinful_hoppings = np.matrix(fun)
  def set_pairing(self,fun):
      """Add the up/down pairing"""
      self.computed_gs = False # say that GS has not been computed
      self.pairing = dict()
      for i in range(self.ns): # loop
          for j in range(self.ns): # loop
              if self.sites[i]==1 and self.sites[j]==1:
                  c = fun(i,j)
                  if np.abs(c)>0.0:
                      self.pairing[(i,j)] = Coupling(i,j,c) # store
  def set_hubbard(self,fun):
      self.computed_gs = False # say that GS has not been computed
      self.hubbard = dict()
      for i in range(self.ns): # loop
          for j in range(self.ns): # loop
              if self.sites[i]==1 and self.sites[j]==1:
                  c = fun(i,j)
                  if np.abs(c)>0.0:
                      self.hubbard[(i,j)] = Coupling(i,j,c) # store
  def set_fields(self,fun):
    self.computed_gs = False # say that GS has not been computed
    self.fields = [fun(i) for i in range(self.ns)] # fields
  def setup_sweep(self,mode="default"):
    setup_sweep(self,mode=mode)
  def setup_task(self,mode="GS",task=dict()):
    from taskdmrg import setup_task
    setup_task(self,mode=mode,task=task)
  def write_hamiltonian(self):
      from writemps import write_hamiltonian
      write_hamiltonian(self)
  def run(self,automatic=False): 
    os.system(dmrgpath+"/mpscpp/mpscpp.x > status.txt") # run the DMRG calculation
  def entropy(self,n=1):
    """Return the entanglement entropy"""
    self.setup_sweep()
    self.setup_task("entropy")
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.run() # perform the calculation
    return np.genfromtxt("ENTROPY.OUT")
  def get_full_hamiltonian(self):
    return pychainwrapper.get_full_hamiltonian(self)
  def get_pychain(self):
    return pychainwrapper.get_pychain(self)
  def get_dos(self,i=0,delta=0.1,window=5.0):
    from dos import get_dos
    return get_dos(self,i=i,delta=delta,window=window)
  def get_spismj(self,n=1000,mode="DMRG",i=0,j=0,smart=False):
    return kpmdmrg.get_spismj(self,n=n,mode=mode,i=i,j=j,smart=smart)
  def get_dynamical_correlator(self,n=1000,mode="DMRG",i=0,j=0,name="XX",
                                 delta=0.02,window=[-1.0,10.0]):
    self.set_initial_wf(self.wf0) # set the initial wavefunction
    return kpmdmrg.get_dynamical_correlator(self,n=n,mode=mode,
                     i=i,j=j,name=name,delta=delta,window=window)
  def get_excited(self,n=10,mode="DMRG"):
    self.to_folder() # go to temporal folder
    if mode=="DMRG":
      self.setup_sweep()
      self.setup_task("excited",task={"nexcited":str(n)})
      self.write_hamiltonian() # write the Hamiltonian to a file
      self.run() # perform the calculation
      out = np.genfromtxt("EXCITED.OUT")
    elif mode=="ED":
      h = self.get_full_hamiltonian() # get the Hamiltonian
      import pychain.spectrum
      out = pychain.spectrum.eigenstates(h,k=n) # return energy
    else: 
      self.to_origin() # go to main folder
      raise
    self.to_origin() # go to main folder
    return out
  def get_gap(self):
    es = self.get_excited(2)
    return es[1] -es[0]
  def set_initial_wf(self,wf):
      """Use a certain wavefunction as initial guess"""
      if wf is None: return
      self.gs_from_file = True # use a wavefunction from a file
      self.starting_file_gs = wf.name # name of the wavefunction
  def get_gs(self,mode="DMRG",wf0=None):
    """Return the ground state"""
    if mode=="DMRG":
      self.gs_energy(wf0=wf0) # perform a ground state calculation
      return self.wf0 # return wavefucntion
    elif mode=="ED":
      self.to_folder() # go to temporal folder
      h = self.get_full_hamiltonian() # get the Hamiltonian 
      self.to_origin() # go to temporal folder
      return pychain.spectrum.ground_state(h)[1] # return energy
    else: raise
  def gs_energy(self,mode="DMRG",wf0=None):
    """Return the ground state energy"""
    # write the sites
    self.to_folder() # go to temporal folder
    if mode=="DMRG":
      self.set_initial_wf(wf0) # set the initial wavefunction
      self.setup_sweep()
      self.setup_task("GS")
      self.write_hamiltonian() # write the Hamiltonian to a file
      self.run() # perform the calculation
      self.wf0 = mps.MPS(self).copy() # set the ground state
      out = np.genfromtxt("GS_ENERGY.OUT") # return the ground state energy
      self.computed_gs = True
    elif mode=="ED": # use brute force
      h = self.get_full_hamiltonian() # get the Hamiltonian 
      out = pychain.spectrum.ground_state(h)[0] # return energy
    else: 
      self.to_origin() # go to main folder
      raise
    self.to_origin() # go to main folder
    self.restart = True
    return out
  def correlator(self,pairs=[[]],mode="DMRG"):
    self.to_folder() # go to temporal folder
    if mode=="DMRG": # DMRG correlation
      self.setup_sweep()
      self.setup_task("correlator")
      self.write_hamiltonian() # write the Hamiltonian to a file
      write_correlators(pairs) # write the input file
      self.run() # perform the calculation
      m = np.genfromtxt("CORRELATORS.OUT").transpose()[1] # return the correlators
    else: raise # not implemented
    self.to_origin() # go to main folder
    return m
  def get_magnetization(self):
    """Calculate the magnetization of the system"""
    mx = self.get_file("MEASURE_SX.OUT").transpose()[1]
    my = self.get_file("MEASURE_SY.OUT").transpose()[1]
    mz = self.get_file("MEASURE_SZ.OUT").transpose()[1]
    return (mx,my,mz)
  def get_file(self,name):
    """Return the electronic density"""
    if not self.computed_gs: self.get_gs() # compute gs
    self.to_folder() # go to folder
    m = np.genfromtxt(name) # read file
    self.to_origin() # go back
    return m



#from fermionchain import Fermionic_Hamiltonian
#from spinchain import Spin_Hamiltonian



from writemps import write_hoppings
from writemps import write_hubbard
from writemps import write_fields
from writemps import write_sites
from writemps import write_couplings
from writemps import write_correlators
from writemps import write_sweeps
from writemps import write_pairing






def setup_sweep(self,mode="default"):
  """Setup the sweep parameters"""
  sweep = dict() # dictionary
  sweep["cutoff"] = 1e-06
  if mode=="default": # default mode
    sweep["n"] = "20"
    sweep["maxm"] = "100" 
  elif mode=="fast": # default mode
    sweep["n"] = "3"
    sweep["maxm"] = "20" 
  elif mode=="accurate": # default mode
    sweep["n"] = "10"
    sweep["maxm"] = "50" 
  else: raise
  sweep["n"] = self.nsweeps
  sweep["maxm"] = self.maxm
  sweep["cutoff"] = self.cutoff
  self.sweep = sweep # initialize
  write_sweeps(self) # write the sweeps





