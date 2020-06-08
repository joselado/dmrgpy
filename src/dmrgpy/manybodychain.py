from __future__ import print_function
import numpy as np
import os
from . import mps
from . import timedependent
from . import groundstate
from . import operatornames
from . import correlator
from . import densitymatrix
from . import taskdmrg
from . import dynamics
from . import funtk
from . import vev
from . import mpsalgebra

#dmrgpath = os.environ["DMRGROOT"]+"/dmrgpy" # path for the program
dmrgpath = os.path.dirname(os.path.realpath(__file__)) # path for the program
one = np.matrix(np.identity(3))



class Coupling():
  def __init__(self,i,j,g):
    self.i = i
    self.j = j
    self.g = g




class Many_Body_Chain():
  def __init__(self,sites):
      self.sites = sites # list of the sites
      self.path = os.getcwd()+"/.mpsfolder/" # folder of the calculations
      self.clean() # clean calculation
      self.inipath = os.getcwd() # original folder
      self.ns = len(sites) # number of sites
      self.exchange = 0 # zero
      self.fields = 0 # zero
      self.pairing = 0
      self.hubbard = 0
      self.hopping = 0
      self.fermionic = False
      self.sites_from_file = False
      self.excited_gram_schmidt = False # it does not seem very effective
      self.hamiltonian = None # Hamiltonian, as a multioperator
      self.hubbard_matrix = np.zeros((self.ns,self.ns)) # empty matrix
      self.use_ampo_hamiltonian = False # use ampo Hamiltonian
  #    self.exchange.append(Coupling(0,self.ns-1,one)) # closed boundary
      # additional arguments
      self.kpmmaxm = 50 # bond dimension in KPM
      self.maxm = 30 # bond dimension in wavefunctions
      self.nsweeps = 15 # number of sweeps
      self.noise = 1e-1 # noise for dmrg
      self.kpmcutoff = 1e-8 # cutoff in KPM
      self.cutoff = 1e-8 # cutoff in ground state
      self.tevol_custom_exp = True # custom exponential function for Tevol
      self.cvm_tol = 1e-5 # tolerance for CVM
      self.cvm_nit = 1e3 # iterations for CVM
      self.kpm_scale = 0.7 # scaling of the spectra for KPM
      self.kpm_accelerate = True # set to true
      self.kpm_n_scale = 3 # scaling factor for the number of polynomials
      self.gs_from_file = False # start from a random wavefunction
      self.excited_from_file = False # read excited states
      self.e0 = None # no ground state energy
      self.wf0 = None # no initial WF
      self.skip_dmrg_gs = False # skip the DMRG minimization
      self.computed_gs = False # computed the GS already
      self.vijkl = 0 # generalized interaction
      self.fit_td = False # use fitting procedure in time evolution
      self.itensor_version = 2 # ITensor version
      self.has_ED_obj = False # ED object has been computed
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
      out.inipath = os.getcwd() # initial path
#      print("New path",out.path)
      os.system("cp -r "+self.path+"  "+out.path) # copy to the new path
      return out # return new object
  def set_hamiltonian(self,MO,restart=True): 
      """Set the Hamiltonian"""
      if restart: self.restart() # restar the calculation
      self.hamiltonian = MO
      self.use_ampo_hamiltonian = True # use ampo Hamiltonian
  def bandwidth(self,h,**kwargs):
      """Compute the bandwidth of an Hermitian operator"""
      mbc = self.clone() # clone the object
      mbc.set_hamiltonian(h) ; e0 = mbc.gs_energy(**kwargs)
      mbc.set_hamiltonian(-h) ; e1 = mbc.gs_energy(**kwargs)
      mbc.clean() # remove
      return -e0 -e1
  def to_origin(self): 
      if os.path.isfile(self.path+"/ERROR"): raise # something wrong
      os.chdir(self.inipath) # go to original folder
  def restart(self):
      """Restart the calculation"""
      self.computed_gs = False
      self.gs_from_file = False
      self.skip_dmrg_gs = False
      self.wf0 = None # initial file for GS
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
  def vev(self,MO,mode="DMRG",**kwargs): 
      if mode=="DMRG": return vev.vev(self,MO,**kwargs)
      elif mode=="ED": return self.get_ED_obj().vev(MO,**kwargs) # ED object
      else: raise
  def test_ED(self):
      """Test the ED object"""
      self.get_ED_obj().test()
  def excited_vev_MB(self,MO,**kwargs):
      """
      Compute a vacuum expectation value
      """
      return vev.excited_vev(self,MO,**kwargs)
  def excited_vev(self,MO,**kwargs): return self.excited_vev_MB(MO,**kwargs)
  def set_vijkl(self,f):
      """
      Create the generalized interaction
      """
      h = 0
      C = self.C
      Cdag = self.Cdag
      for i in range(self.ns):
        for j in range(self.ns):
          for k in range(self.ns):
            for l in range(self.ns):
                h = h + f(i,j,k,l)*Cdag[i]*C[j]*Cdag[k]*C[l]
      h = 0.5*(h+h.get_dagger())
      self.vijkl = h # store
      self.update_hamiltonian()
  def generate_bilinear(self,fun,A,B):
      """Generic bilinear term"""
      fun = funtk.obj2fun(fun) # set function
      h = 0 # initialize
      for i in range(self.ns): # loop
          for j in range(self.ns): h = h + fun(i,j)*A[i]*B[j]
      return 0.5*(h + h.get_dagger()) # Hermitian
  def update_hamiltonian(self):
      h = self.hopping + self.hubbard + self.pairing 
      h = h + self.vijkl + self.exchange
      self.set_hamiltonian(h)
  def overlap(self,wf1,wf2,**kwargs):
      """Compute the overlap"""
      return mpsalgebra.overlap(self,wf1,wf2,**kwargs)
  def exponential(self,h,wf,**kwargs):
      """Compute the overlap"""
      return mpsalgebra.exponential(self,h,wf,**kwargs)
  def applyoperator(self,A,wf,**kwargs):
      """Apply an operator"""
      return mpsalgebra.applyoperator(self,A,wf,**kwargs)
  def summps(self,wf1,wf2,**kwargs):
      """Apply an operator"""
      return mpsalgebra.summps(self,wf1,wf2,**kwargs)
  def set_pairings_MB(self,fun):
      """Generic pairing term"""
      h = self.generate_bilinear(fun,self.C,self.C)
      self.pairing = h
      self.update_hamiltonian()
  def set_hubbard_MB(self,fun):
      """Hubbard term"""
      h = self.generate_bilinear(fun,self.N,self.N)
      self.hubbard = h # store
      self.update_hamiltonian()
  def set_hoppings_MB(self,fun):
      """Hubbard term"""
      h = self.generate_bilinear(fun,self.Cdag,self.C)
      self.hopping = h # store
      self.update_hamiltonian()
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
      self.execute(lambda : taskdmrg.write_tasks(self)) # write tasks
      if self.itensor_version in [2,"2","v2","C++","cpp","c","C"]: 
          mpscpp = dmrgpath+"/mpscpp2/mpscpp.x" 
#      elif self.itensor_version==3: mpscpp = dmrgpath+"/mpscpp3/mpscpp.x" 
      elif self.itensor_version in ["julia","Julia","jl"]: 
          from . import juliarun
          juliarun.run(self)
          return
     #     mpscpp = "julia "+dmrgpath+"/mpsjulia/mpsjulia.jl" 
      else: raise
      if not os.path.isfile(mpscpp): # mpscpp.x not found, rerun with julia
          print("C++ backend not found, trying tu run with Julia version")
          self.itensor_version = "julia" # turn to Julia
          return self.run() # rerun with julia
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
  def get_dynamical_correlator_MB(self,**kwargs):
      return dynamics.get_dynamical_correlator(self,**kwargs)
  def get_dynamical_correlator(self,mode="DMRG",**kwargs):
      if mode=="DMRG": 
          return dynamics.get_dynamical_correlator(self,**kwargs)
      elif mode=="ED": 
          return self.get_ED_obj().get_dynamical_correlator(**kwargs)
      else: raise
  def get_distribution(self,mode="DMRG",**kwargs):
      if mode=="DMRG": 
          from . import distribution
          return distribution.get_distribution(self,**kwargs)
     #     raise # not implemented
       #   return dynamics.get_dynamical_correlator(self,**kwargs)
      elif mode=="ED": 
          return self.get_ED_obj().get_distribution(**kwargs)
      else: raise
  def get_excited(self,mode="DMRG",**kwargs):
      """Return excitation energies"""
      if mode=="DMRG":
        from . import excited
        return excited.get_excited(self,**kwargs) # return excitation energies
      elif mode=="ED": return self.get_ED_obj().get_excited(**kwargs) # ED
      else: raise
  def get_gap(self,**kwargs):
    """Return the gap"""
    es = self.get_excited(n=2,**kwargs)
    return es[1] -es[0]
  def get_hamiltonian(): return self.hamiltonian
  def gs_energy_fluctuation(self,**kwargs):
      """Compute the energy fluctuations"""
      h = self.get_hamiltonian()
      e = self.vev(h)
      e2 = self.vev(h*h)
      return np.sqrt(np.abs(e2-e**2))
  def set_initial_wf_guess(self,wf):
      """Set the initial guess, and perform the DMRG GS calculation"""
      self.set_initial_wf(wf,reconverge=True)
  def set_initial_wf(self,wf,reconverge=False):
      """Use a certain wavefunction as initial guess"""
      if wf is None:
        self.gs_from_file = False # use a wavefunction from a file
      else:
        self.gs_from_file = True # use a wavefunction from a file
        self.wf0 = wf.copy() # name of the wavefunction
        if reconverge: self.skip_dmrg_gs = False # reconverge the calculation
        else: self.skip_dmrg_gs = True # reconverge the calculation
  def get_gs(self,best=False,n=1,mode="DMRG",**kwargs):
      """Return the ground state"""
      if mode=="DMRG":
        if best: groundstate.best_gs(self,n=n,**kwargs) # best ground state
        else: self.gs_energy(**kwargs) # perform a ground state calculation
        return self.wf0 # return wavefunction
      elif mode=="ED": return self.get_ED_obj().get_gs()
  def gs_energy(self,mode="DMRG",**kwargs):
      """Return the ground state energy"""
      if mode=="DMRG": return groundstate.gs_energy(self,**kwargs)
      elif mode=="ED": return self.get_ED_obj().gs_energy() # ED object
      else: raise
#  def get_correlator_MB(self,**kwargs):
#      """Return a correlator"""
#      return correlator.get_correlator(self,**kwargs)
  def get_correlator(self,pairs=[],**kwargs):
      """Return a correlator, default one"""
      print("Method get_correlator is deprecated, use vev instead")
      from . import spinchain
      from . import fermionchain
      if type(self)==spinchain.Spin_Chain: 
          ops = [self.Sz[i]*self.Sz[j] for (i,j) in pairs]
      elif type(self)==fermionchain.Fermionic_Chain: 
          ops = [self.Cdag[i]*self.C[j] for (i,j) in pairs]
      else: raise
      return np.array([self.vev(o,**kwargs) for o in ops])
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
  def get_operator(self,name,i=None):
      """Return a certain multioperator"""
      from . import multioperator
      return multioperator.obj2MO([[name,i]]) # return operator



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




Many_Body_Hamiltonian = Many_Body_Chain # temporal workaround

