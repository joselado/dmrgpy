from __future__ import print_function
import numpy as np
import os
from . import kpmdmrg
from . import mps
from . import timedependent
from . import groundstate
from . import operatornames
from . import correlator
from . import densitymatrix
from . import taskdmrg
from . import cvm
from . import dcex
from . import funtk
from . import vev

#dmrgpath = os.environ["DMRGROOT"]+"/dmrgpy" # path for the program
dmrgpath = os.path.dirname(os.path.realpath(__file__)) # path for the program
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
      self.clean() # clean calculation
      self.inipath = os.getcwd() # original folder
      self.ns = len(sites) # number of sites
      self.exchange = [Coupling(i,i+1,one) for i in range(self.ns-1)] # empty list
      self.exchange = [] # empty list
      self.fields = [] # empty list
      self.hoppings = dict() # empty dictionary
      self.spinful_hoppings = dict() # empty dictionary
      self.pairing = None # pairing
      self.hubbard = dict() # empty dictionary
      self.hubbard_matrix = np.zeros((self.ns,self.ns)) # empty matrix
  #    self.exchange.append(Coupling(0,self.ns-1,one)) # closed boundary
      # additional arguments
      self.kpmmaxm = 50 # bond dimension in KPM
      self.maxm = 30 # bond dimension in wavefunctions
      self.nsweeps = 15 # number of sweeps
      self.kpmcutoff = 1e-8 # cutoff in KPM
      self.cutoff = 1e-8 # cutoff in ground state
      self.cvm_tol = 1e-5 # tolerance for CVM
      self.cvm_nit = 1e3 # iterations for CVM
      self.kpm_scale = 0.6 # scaling of the spectra for KPM
      self.kpm_accelerate = True # set to true
      self.kpm_n_scale = 3 # scaling factor for the number of polynomials
      self.restart = False # restart the calculation
      self.gs_from_file = False # start from a random wavefunction
      self.e0 = None # no ground state energy
      self.wf0 = None # no initial WF
      self.starting_file_gs = None # initial file for GS
      self.skip_dmrg_gs = False # skip the DMRG minimization
      self.computed_gs = False # computed the GS already
      self.vijkl = None # generalized interaction
      self.fit_td = False # use fitting procedure in time evolution
      self.itensor_version = 2 # ITensor version
      self.hamiltonian_multioperator = None # Multioperator for the Hamiltonian
      os.system("mkdir -p "+self.path) # create folder for the calculations
  def to_folder(self):
      """Go to a certain folder"""
#      self.inipath = os.getcwd() # record the folder
      os.chdir(self.path) # go to calculation folder
  def copy(self):
      from copy import deepcopy
      return deepcopy(self)
  def clone(self):
      """
      Clone the object and create a temporal folder
      """
      from copy import deepcopy
      name = "dmrgpy_clone_"+str(np.random.randint(100000))
      out = deepcopy(self) 
      out.path = "/tmp/"+name # new path
      out.inipath = out.path # initial path
      print("New path",out.path)
      os.system("cp -r "+self.path+"  "+out.path) # copy to the new path
      return out # return new object
  def to_origin(self): 
    if os.path.isfile(self.path+"/ERROR"): raise # something wrong
    os.chdir(self.inipath) # go to original folder
  def clean(self): 
      """
      Remove the temporal folder
      """
      os.system("rm -rf "+self.path) # clean temporal folder
  def vev_MB(self,MO,**kwargs):
      """
      Compute a vacuum expectation value
      """
      return vev.vev(self,MO,**kwargs)
  def vev(self,MO,**kwargs):
      """
      Compute a vacuum expectation value
      """
      return vev.vev(self,MO,**kwargs)
  def set_vijkl(self,f):
      """
      Create the generalized interaction
      """
      self.vijkl = f # store
  def set_spinful_hoppings(self,fun):
      """Add the spin independent hoppings"""
      raise # no longer used
#      self.computed_gs = False # say that GS has not been computed
#      if callable(fun):
#        self.spinful_hoppings = dict()
#        for i in range(self.ns): # loop
#            for j in range(self.ns): # loop
#                if self.sites[i]==1 and self.sites[j]==1:
#                    c = fun(i,j)
#                    if np.abs(c)>0.0:
#                        self.spinful_hoppings[(i,j)] = Coupling(i,j,c) # store
#      else: # assume it is a matrix
#          self.spinful_hoppings = np.matrix(fun)
  def set_pairings_MB(self,fun):
      """Add the up/down pairing"""
      self.computed_gs = False # say that GS has not been computed
      self.pairing = funtk.obj2fun(fun) # set function
#      self.pairing = dict()
#      for i in range(self.ns): # loop
#          for j in range(self.ns): # loop
#              if self.sites[i]==1 and self.sites[j]==1:
#                  c = fun(i,j)
#                  if np.abs(c)>0.0:
#                      self.pairing[(i,j)] = Coupling(i,j,c) # store
  def set_hubbard_MB(self,fun):
      self.computed_gs = False # say that GS has not been computed
      self.hubbard = dict()
      for i in range(self.ns): # loop
          for j in range(self.ns): # loop
              if self.sites[i] in [0,1] and self.sites[j] in [0,1]:
                  c = fun(i,j)
                  if np.abs(c)>0.0:
                      self.hubbard[(i,j)] = Coupling(i,j,c) # store
      m = np.zeros((self.ns,self.ns)) # initialize
      for key in self.hubbard:
          c = self.hubbard[key]
          m[c.i,c.j] = c.g # create entry
      self.hubbard_matrix = m # store matrix
  def set_fields(self,fun):
    """
    Add local magnetic fields
    """
    self.computed_gs = False # say that GS has not been computed
    self.fields = [fun(i) for i in range(self.ns)] # fields
#  def setup_sweep(self,mode="default"):
#    setup_sweep(self,mode=mode)
  def setup_task(self,mode="GS",task=dict()):
    from .taskdmrg import setup_task
    setup_task(self,mode=mode,task=task)
  def write_task(self):
      """
      Write the tasks in tasks.in
      """
      self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
  def write_hamiltonian(self):
      """
      Write the Hamiltonian in a file
      """
      from .writemps import write_hamiltonian
      self.execute(lambda: write_hamiltonian(self))
  def run(self,automatic=False): 
      """
      Run the DMRG calculation
      """
      # executable
      if self.itensor_version==2: mpscpp = dmrgpath+"/mpscpp2/mpscpp.x" 
      elif self.itensor_version==3: mpscpp = dmrgpath+"/mpscpp3/mpscpp.x" 
      else: raise
      if not os.path.isfile(mpscpp): raise
      self.execute(lambda : os.system(mpscpp+" > status.txt"))
  def entropy(self,n=1):
    """Return the entanglement entropy"""
#    self.setup_sweep()
    self.setup_task("entropy")
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.run() # perform the calculation
    return np.genfromtxt("ENTROPY.OUT")
  def get_dos(self,**kwargs):
    from .dos import get_dos
    return get_dos(self,**kwargs)
  def get_spismj(self,n=1000,mode="DMRG",i=0,j=0,smart=False):
    return kpmdmrg.get_spismj(self,n=n,mode=mode,i=i,j=j,smart=smart)
  def get_dynamical_correlator_MB(self,submode="KPM",**kwargs):
    self.set_initial_wf(self.wf0) # set the initial wavefunction
    if submode=="KPM": # KPM method
        return kpmdmrg.get_dynamical_correlator(self,**kwargs)
    elif submode=="TD": # time dependent 
        return timedependent.dynamical_correlator(self,**kwargs)
    elif submode=="CVM": # CVM mode
        return cvm.dynamical_correlator(self,**kwargs)
    elif submode=="EX": # CVM mode
        return dcex.dynamical_correlator(self,**kwargs)
    else: raise
  def get_excited(self,**kwargs):
    """Return excitation energies"""
    from . import excited
    return excited.get_excited(self,**kwargs) # return excitation energies
  def get_gap(self):
    """Return the gap"""
    es = self.get_excited(2)
    return es[1] -es[0]
  def set_initial_wf(self,wf):
      """Use a certain wavefunction as initial guess"""
      if wf is None:
        self.gs_from_file = False # use a wavefunction from a file
      else:
        self.gs_from_file = True # use a wavefunction from a file
        self.starting_file_gs = wf.name # name of the wavefunction
  def get_gs(self,best=False,n=1,**kwargs):
      """Return the ground state"""
      if best: groundstate.best_gs(self,n=n,**kwargs) # best ground state
      else: self.gs_energy(**kwargs) # perform a ground state calculation
      return self.wf0 # return wavefucntion
  def gs_energy(self,**kwargs):
      """Return the ground state energy"""
      return groundstate.gs_energy(self,**kwargs)
  def get_correlator_MB(self,**kwargs):
      """Return a correlator"""
      return correlator.get_correlator(self,**kwargs)
  def get_correlator(self,**kwargs):
      """Return a correlator, default one"""
      return correlator.get_correlator(self,**kwargs)
  def get_file(self,name):
      """Return the electronic density"""
      if not self.computed_gs: self.get_gs() # compute gs
      self.to_folder() # go to folder
      m = np.genfromtxt(name) # read file
      self.to_origin() # go back
      return m
  def execute(self,f):
      """Execute function in the folder"""
      self.to_folder() # go to folder
      out = f()
      self.to_origin() # go back
      return out # return result
  def evolution(self,**kwargs):
      """
      Perform time dependent DMRG
      """
      from . import timedependent
      return timedependent.evolution(self,**kwargs)
  def get_rdm(self,**kwargs):
      """
      Compute the reduced density matrix
      """
      return densitymatrix.reduced_dm(self,**kwargs) # return DM
  def get_kpm_scale(self):
      """
      Return an estimate of the bandwidth
      """
      return 3*self.ns # estimated bandwidth



#from fermionchain import Fermionic_Hamiltonian
#from spinchain import Spin_Hamiltonian



from .writemps import write_hoppings
from .writemps import write_hubbard
from .writemps import write_fields
from .writemps import write_sites
from .writemps import write_exchange
#from .writemps import write_sweeps
from .writemps import write_pairing






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
#  write_sweeps(self) # write the sweeps





