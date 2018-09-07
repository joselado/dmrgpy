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
    self.hubbard = dict() # empty dictionary
#    self.couplings.append(Coupling(0,self.ns-1,one)) # closed boundary
    # additional arguments
    self.kpmmaxm = 50 # bond dimension in KPM
    self.maxm = 30 # bond dimension in wavefunctions
    self.nsweeps = 8 # number of sweeps
    self.kpmcutoff = 1e-5 # cutoff in KPM
    self.kpmscale = 10.0
    self.restart = False # restart the calculation
    os.system("mkdir -p "+self.path) # create folder for the calculations
  def to_folder(self): os.chdir(self.path) # go to calculation folder
  def to_origin(self): 
    if os.path.isfile(self.path+"/ERROR"): raise # something wrong
    os.chdir(self.inipath) # go to original folder
  def clean(self): os.system("rm -rf "+self.path) # clean temporal folder
  def set_exchange(self,fun):
    """Set the exchange coupling between sites"""
    self.couplings = [] # empty list
    for i in range(self.ns): # loop
      for j in range(self.ns):  # loop
        g = fun(i,j) # call the function
        g = g*one # multiply by the identity
        if np.sum(np.abs(g))!=0.0: 
          c = Coupling(i,j,g) # create class
          self.couplings.append(c) # store
  def set_hoppings(self,fun):
      self.hoppings = dict()
      for i in range(self.ns): # loop
          for j in range(self.ns): # loop
              if self.sites[i]==1 and self.sites[j]==1:
                  c = fun(i,j)
                  if np.abs(c)>0.0:
                      self.hoppings[(i,j)] = Coupling(i,j,c) # store
  def set_hubbard(self,fun):
      self.hubbard = dict()
      for i in range(self.ns): # loop
          for j in range(self.ns): # loop
              if self.sites[i]==1 and self.sites[j]==1:
                  c = fun(i,j)
                  if np.abs(c)>0.0:
                      self.hubbard[(i,j)] = Coupling(i,j,c) # store
  def set_fields(self,fun):
    self.fields = [fun(i) for i in range(self.ns)] # fields
  def get_coupling(self,i,j):
    """Return the coupling between two sites"""
    for c in self.couplings:
      if i==c.i and j==c.j: return c.g
    return np.zeros((3,3)) 
  def setup_sweep(self,mode="default"):
    setup_sweep(self,mode=mode)
  def setup_task(self,mode="GS",task=dict()):
    from taskdmrg import setup_task
    setup_task(self,mode=mode,task=task)
  def write_hamiltonian(self):
    write_sites(self) # write the different sites
    write_couplings(self)  # write the couplings
    write_hoppings(self)  # write the hoppings
    write_hubbard(self)  # write hubbard terms
    write_fields(self) # write the fields
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
  def get_gs(self,mode="DMRG"):
    """Return the ground state"""
    if mode=="DMRG":
      self.gs_energy() # perform a ground state calculation
      wf = mps.MPS(self) # create an MPS
      return wf.copy() # return wavefucntion
    elif mode=="ED":
      self.to_folder() # go to temporal folder
      h = self.get_full_hamiltonian() # get the Hamiltonian 
      self.to_origin() # go to temporal folder
      return pychain.spectrum.ground_state(h)[1] # return energy
    else: raise
  def gs_energy(self,mode="DMRG"):
    """Return the ground state energy"""
    # write the sites
    self.to_folder() # go to temporal folder
    if mode=="DMRG":
      self.setup_sweep()
      self.setup_task("GS")
      self.write_hamiltonian() # write the Hamiltonian to a file
      self.run() # perform the calculation
      out = np.genfromtxt("GS_ENERGY.OUT") # return the ground state energy
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
  def magnetization(self):
    """Calculate the magnetization of the system"""
    self.gs_energy() # calculate ground state
    self.to_folder() # go to temporal folder
    mx = np.genfromtxt("MEASURE_SX.OUT").transpose()[1]
    my = np.genfromtxt("MEASURE_SY.OUT").transpose()[1]
    mz = np.genfromtxt("MEASURE_SZ.OUT").transpose()[1]
    self.to_origin() # go to main folder
    return (mx,my,mz)


#from fermionchain import Fermionic_Hamiltonian
#from spinchain import Spin_Hamiltonian



from writemps import write_hoppings
from writemps import write_hubbard
from writemps import write_fields
from writemps import write_sites
from writemps import write_couplings
from writemps import write_correlators
from writemps import write_sweeps






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
  self.sweep = sweep # initialize
  write_sweeps(self) # write the sweeps


