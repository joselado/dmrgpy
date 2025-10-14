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
from . import entanglement
from . import entropy
from . import excited
from . import effectivehamiltonian
from .writemps import write_sites
from .mode import dmrgpath
import subprocess

one = np.matrix(np.identity(3))



class Coupling():
  def __init__(self,i,j,g):
    self.i = i
    self.j = j
    self.g = g




class Many_Body_Chain():
  def __init__(self,sites,**kwargs):
      self.sites = sites # list of the sites
#      self.path = id_generator() # random ID in dmrgpy_tmp
      self.ns = len(sites) # number of sites
      self.mode = None # no mode (use the input parameter)
      self.exchange = 0 # zero
      self.fields = 0 # zero
      self.pairing = 0
      self.hubbard = 0
      self.hopping = 0
      self.resorder = False # reorder the indexes
      self.resordered_indexes = None # reordered indexes
      self.fermionic = False
      self.sites_from_file = False
      self.excited_gram_schmidt = False # it does not seem very effective
      self.hamiltonian = None # Hamiltonian, as a multioperator
      self.hubbard_matrix = np.zeros((self.ns,self.ns)) # empty matrix
      self.use_ampo_hamiltonian = False # use ampo Hamiltonian
      # additional arguments
      self.kpmmaxm = 50 # bond dimension in KPM
      self.maxm = 30 # bond dimension in wavefunctions
      self.mpomaxm = 5000 # bond dimension for operators
      self.nsweeps = 15 # number of sweeps
      self.noise = 1e-1 # noise for dmrg
      self.kpmcutoff = 1e-12 # cutoff in KPM
      self.cutoff = 1e-12 # cutoff in ground state
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
      self.ED_obj = None # no ED object
      self.kpm_extrapolate = False # use extrapolation
      self.kpm_extrapolate_factor = 2.0 # factor for the extrapolation
      self.kpm_extrapolate_mode = "plain" # mode of the extrapolation
      self.initialize(**kwargs)
      # and initialize the sites
  def initialize(self,**kwargs):
      """Initialize the sites"""
      if self.mode=="ED": return # do nothing
      if self.itensor_version in [2,"julia"]:
          from .sites import initialize
          initialize(self)
      elif self.itensor_version=="julia_live":
          from .mpsjulialive.sites import initialize
          initialize(self)
      else: raise
  def setup_julia(self):
      """Setup the Julia mode"""
      self.itensor_version = "julia_live"
      self.initialize()
  def setup_cpp(self):
      """Setup the C++ mode"""
      self.itensor_version = 2
      self.initialize()
  def get_mode(self,**kwargs):
      from .mode import get_mode
      return get_mode(self,**kwargs)
  def to_folder(self):
      """Go to a certain folder"""
#      self.inipath = os.getcwd() # record the folder
      os.chdir(self.path) # go to calculation folder
  def copy(self):
      return self.clone() # clone and create a new one
      from copy import deepcopy
      return deepcopy(self)
  def clone(self):
      """
      Clone the object and create a temporal folder
      """
      from copy import deepcopy
      name = "dmrgpy_clone_"+str(np.random.randint(10000))
      subprocess.run(["rm","-rf","/tmp/"+name]) # clean the new directory
      out = deepcopy(self) # full copy of the object 
      out.path = "/tmp/"+name # new path
      out.inipath = os.getcwd() # initial path
#      print("New path",out.path)
      subprocess.run(["cp","-r",self.path,out.path]) # copy to the new path
      return out # return new object
  def set_hamiltonian(self,MO,restart=True): 
      """Set the Hamiltonian"""
      if restart: self.restart() # restar the calculation
      self.hamiltonian = MO
      self.use_ampo_hamiltonian = True # use ampo Hamiltonian
  def get_heff(self,**kwargs):
      """Return effective Hamiltonian"""
      return effectivehamiltonian.get_effective_hamiltonian_coefficients(self,
              **kwargs)
  def bandwidth(self,h,**kwargs):
      """Compute the bandwidth of an Hermitian operator"""
      mbc = self.clone() # clone the object
      mbc.set_hamiltonian(h) ; e0 = mbc.gs_energy(**kwargs)
      mbc.set_hamiltonian(-h) ; e1 = mbc.gs_energy(**kwargs)
      mbc.clean() # remove
      return -e0 -e1
  def lowest_eigenvalue(self,X,**kwargs):
      """Given an operator X, return its smallest eigenvalue"""
      mbc = self.clone() # clone the object
      mbc.set_hamiltonian(X) 
      return mbc.gs_energy(**kwargs)
  def to_origin(self): 
      if os.path.isfile(self.path+"/ERROR"): raise # something wrong
      os.chdir(self.inipath) # go to original folder
  def restart(self):
      """Restart the calculation"""
      self.computed_gs = False
      self.gs_from_file = False
      self.has_ED_obj = False # restart ED obj
      self.skip_dmrg_gs = False
      self.wf0 = None # initial file for GS
  def is_hermitian(self,H):
      """Check if an operator is Hermitian"""
      from .mpsalgebra import is_hermitian
      return is_hermitian(self,H)
  def clean(self): 
      """
      Remove the temporal folder
      """
      subprocess.run(["rm","-rf",self.path]) # clean temporal folder
  def vev_MB(self,MO,**kwargs):
      """
      Compute a vacuum expectation value
      """
      return vev.vev(self,MO,**kwargs)
  def get_gs_degeneracy(self,**kwargs):
      from . import degeneracy
      return degeneracy.gs_degeneracy(self,**kwargs)
  def vev(self,MO,mode="DMRG",**kwargs): 
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG": 
          if self.itensor_version==2: # C++ version
              return vev.vev(self,MO,**kwargs)
          elif self.itensor_version=="julia_live": # julia live version
              from .mpsjulialive.vev import vev as vevjl
              return vevjl(self,MO,**kwargs)
          else: raise
      elif mode=="ED": return self.get_ED_obj().vev(MO,**kwargs) # ED object
      else: raise
  def test_ED(self):
      """Test the ED object"""
      self.get_ED_obj().test()
  def toMPO(self,H,**kwargs):
      """Transport an operator into a matrix-product operator"""
      return mpsalgebra.toMPO(self,H,**kwargs)
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
  def get_dagger(self,m): return m.get_dagger() # dummy method
  def overlap(self,wf1,wf2,**kwargs):
      """Compute the overlap"""
      return mpsalgebra.overlap(self,wf1,wf2,**kwargs)
  def aMb(self,wf1,M,wf2,**kwargs):
      """Compute the overlap <a|M|b>"""
      return mpsalgebra.overlap_aMb(self,wf1,M,wf2,**kwargs)
  def operator_norm(self,op,**kwargs):
      """Estimate the norm of an operator"""
      return mpsalgebra.operator_norm(self,op,**kwargs)
  def is_zero_operator(self,op,**kwargs):
      """Check if this is the zero operator"""
      out = self.operator_norm(op,**kwargs)
      return out<1e-4
  def exponential(self,h,wf,**kwargs):
      """Compute the overlap"""
      return mpsalgebra.exponential(self,h,wf,**kwargs)
  def applyoperator(self,A,wf,**kwargs):
      """Apply an operator"""
      return mpsalgebra.applyoperator(self,A,wf,**kwargs)
  def applyinverse(self,A,wf,**kwargs):
      """Apply an operator"""
      return mpsalgebra.applyinverse(self,A,wf,**kwargs)
  def summps(self,wf1,wf2,**kwargs):
      """Apply an operator"""
      return mpsalgebra.summps(self,wf1,wf2,**kwargs)
  def trace(self,A,**kwargs):
      return mpsalgebra.trace(self,A,**kwargs)
  def inverse_trace(self,A,**kwargs):
      return mpsalgebra.inverse_trace(self,A,**kwargs)
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
      from .writemps import write_sites
      self.execute(lambda: write_sites(self)) # write the different sites
      self.execute(lambda: self.hamiltonian.write("hamiltonian.in"))
  def run(self,**kwargs): 
      from .mode import run
      return run(self,**kwargs)
  def get_bond_entropy(self,wf,i,j):
      """Return the entanglement entropy of two sites"""
      return entropy.bond_entropy(self,wf,i,j)
  def get_site_entropy(self,wf,b):
      """Return the entanglement entropy of a site"""
      return entropy.site_entropy(self,wf,b)
  def get_mutual_information(self,wf,i,j):
      """Return the mutual information"""
      return entropy.mutual_information(self,wf,i,j)
  def get_pair_entropy(self,wf,i,j):
      """Return the entanglement entropy of two sites with the
      rest of the system"""
      return entropy.pair_entropy(self,wf,i,j)
  def get_correlation_matrix(self,**kwargs):
      return entanglement.get_correlation_matrix(self,**kwargs)
  def get_correlation_eigenvalues(self,**kwargs):
      return entanglement.get_correlation_eigenvalues(self,**kwargs)
  def get_correlation_entropy(self,**kwargs):
      return entanglement.get_correlation_entropy(self,**kwargs)
  def get_correlated_orbitals(self,**kwargs):
      return entanglement.get_correlated_orbitals(self,**kwargs)
  def get_correlated_density(self,**kwargs):
      return entanglement.get_correlated_density(self,**kwargs)
  def get_dynamical_correlator_MB(self,**kwargs):
      return dynamics.get_dynamical_correlator(self,**kwargs)
  def get_dynamical_correlator(self,mode="DMRG",**kwargs):
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG": 
          return dynamics.get_dynamical_correlator(self,**kwargs)
      elif mode=="ED": 
          return self.get_ED_obj().get_dynamical_correlator(**kwargs)
  def get_distribution(self,mode="DMRG",**kwargs):
      if mode=="DMRG": 
          from . import distribution
          return distribution.get_distribution(self,**kwargs)
     #     raise # not implemented
       #   return dynamics.get_dynamical_correlator(self,**kwargs)
      elif mode=="ED": 
          return self.get_ED_obj().get_distribution(**kwargs)
      else: raise
  def get_distribution_moments(self,mode="DMRG",**kwargs):
      if mode=="DMRG":
          from . import distribution
          return distribution.get_distribution_moments(self,**kwargs)
      elif mode=="ED":
          raise
      else: raise
  def get_excited(self,mode="DMRG",**kwargs):
      """Return excitation energies"""
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG":
          return excited.get_excited(self,**kwargs) # return excitation energies
      elif mode=="ED": 
          return self.get_ED_obj().get_excited(**kwargs) # ED
  def get_full_matrix(self,name):
      return self.get_ED_obj().get_operator(name) # get the full operator
  def get_excited_states(self,mode="DMRG",**kwargs):
      """Return excitation energies"""
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG":
          return excited.get_excited_states(self,**kwargs) # return es and waves
      elif mode=="ED": 
          return self.get_ED_obj().get_excited_states(**kwargs) # ED
  def get_gap(self,**kwargs):
    """Return the gap"""
    es = self.get_excited(n=2,**kwargs)
    return es[1] -es[0]
  def get_hamiltonian(): return self.hamiltonian
  def gs_energy_fluctuation(self,**kwargs):
      """Compute the energy fluctuations"""
      h = self.get_hamiltonian()
      e = self.vev(h)
      e2 = self.vev(h,npow=2)
      return np.sqrt(np.abs(e2-e**2))
  def set_initial_wf_guess(self,wf):
      """Set the initial guess, and perform the DMRG GS calculation"""
      self.set_initial_wf(wf,reconverge=True)
  def set_initial_wf(self,wf,reconverge=False):
      """Use a certain wavefunction as initial guess"""
      self.computed_gs = False
      if wf is None:
        self.gs_from_file = False # use a wavefunction from a file
      else:
        self.gs_from_file = True # use a wavefunction from a file
        self.wf0 = wf.copy() # name of the wavefunction
        if reconverge: self.skip_dmrg_gs = False # reconverge the calculation
        else: self.skip_dmrg_gs = True # reconverge the calculation
  def set_gs(self,wf):
      """Set the ground state"""
      from .groundstate import set_gs 
      groundstate.set_gs(self,wf) # set this as ground state
  def get_gs(self,best=False,n=1,mode="DMRG",**kwargs):
      """Return the ground state"""
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG": # DMRG mode
        if self.computed_gs: # if stored, rewrite and return
#            self.wf0.write(name=self.wf0.name,path=self.path)
            return self.wf0
        if best: groundstate.best_gs(self,n=n,**kwargs) # best ground state
        else: self.gs_energy(**kwargs) # perform a ground state calculation
        return self.wf0 # return wavefunction
      elif mode=="ED": return self.get_ED_obj().get_gs(**kwargs)
  def get_gs_manifold(self,**kwargs):
      return groundstate.get_gs_manifold(self,**kwargs)
  def gs_energy(self,mode="DMRG",**kwargs):
      """Return the ground state energy"""
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG": 
          if self.computed_gs: return self.e0
          return groundstate.gs_energy(self,**kwargs)
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
  def random_state(self,mode="DMRG",orthogonal=None):
      """Generate a random MPS"""
      if self.mode is not None: mode = self.mode # redefine
      if mode in ["DMRG","MPS"]:
         if self.itensor_version==2: # C++ version
             from . import mps
             if orthogonal is None: return mps.random_mps(self)
             else: return mps.orthogonal_random_mps(self,orthogonal)
         elif self.itensor_version=="julia_live": # Julia version
             from .mpsjulialive import mps
             return mps.random_mps(self)
      elif mode=="ED":
          return self.get_ED_obj().random_state()
      else: raise
  def random_mps(self,**kwargs): return self.random_state(**kwargs)
  def get_operator(self,name,i=None):
      """Return a certain multioperator"""
      from . import multioperator
      return multioperator.obj2MO([[name,i]]) # return operator



#from fermionchain import Fermionic_Hamiltonian
#from spinchain import Spin_Hamiltonian



from .writemps import write_hoppings
from .writemps import write_hubbard
from .writemps import write_fields
from .writemps import write_exchange
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


import string
import random

def id_generator(size=20, chars=string.ascii_uppercase + string.digits):
   out = ''.join(random.choice(chars) for _ in range(size))
   return "/tmp/dmrgpy_tmp/"+out # temporal folder




