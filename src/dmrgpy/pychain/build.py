
from __future__ import print_function
import os
from scipy.sparse import csc_matrix as csc
import scipy.sparse as sparse
from . import spectrum
import numpy as np
from . import read
from . import states
from ..algebra import algebra
from .. import multioperator

usecpp = False # use c++ library


maxsize = 300000 # maximum matrix size
tol = 0.0001
pathlib = "../lib" # path to library


def get_dimension(spins):
  size = 1
  for s in spins:
    size *= 2*s+1
  return int(size)

def generate_inputs(spins=[]):
  """Generate inputs for a certain spin chain,
     inputs are the different spins and couplings"""
  fs = open("spins.in","w")
  fs.write(str(len(spins))+"\n") # number of spins
  size = 1 # size
  for s in spins:
    m = int(round(2*s))+1
    size *= m # accumulate
    fs.write(str(m)+"\n") # number of m components
  if size>maxsize: 
    print("Too large matrices",size,maxsize)
    raise
  fs.close()

def run():
  """ Run the program"""
  os.system("rm -f *.op")
  os.system("./main.x")




def generate_hamiltonian(ms=[],js=[]):
  """ Build a hamiltonian"""
  h = ms[0]*0.0 # initialize hamiltonian
  for (m,j) in zip(ms,js):
    t = j*m # couplings matrix
    h = h + t +t.H # plus the transpose
  return h



def exp_val(v,A):
  """Calculate a certain expectation value"""
  return spectrum.exp_val(v,A) # use the default library




def get_manifolds(evals,evecs):
  """ Return a list with the different manifolds, splitted
  by energy """
  me = min(evals) # minimum
  wfm = [] # list for the wavefunctions in this manifold
  evalsrec = [] # list for eigenvalues left
  evecsrec = [] # list for eigenfunctions left
  for (a,v) in zip(evals,evecs):
    if abs(a-me)<tol:
      wfm.append(v) # append wavefunction
    else:
      evalsrec.append(a) # store eigval left
      evecsrec.append(v) # eigfun left
  if len(evalsrec)>0:
    return [wfm] + get_manifolds(evalsrec,evecsrec) # if still wfs, iterate
  else:
    return [wfm] # return this manifold



def disentangle_manifold(wfs,A):
  """ Disentangles the wavefunctions of a degenerate manifold
  by expressing them in terms of eigenvalues of an input operator"""
  import scipy.linalg as dlg
  ma = get_representation(wfs,A) # get the matrix form of the operator
  wfsout = [] # empty list
  evals,evecs = dlg.eigh(ma) # diagonalize

  evecs = evecs.transpose() # transpose eigenvectors
  for v in evecs: # loop over eigenvectors
    wf = wfs[0]*0.0j
    for (i,iv) in zip(range(len(v)),v): # loop over components
      wf += iv*wfs[i] # add contribution
    wfsout.append(wf.copy()) # store wavefunction
  return wfsout


def disentangle_all(evals,evecs,A):
  """Disentangles eigenvalues and eigenvectors"""
  mfs = get_manifolds(evals,evecs) # get the different manifolds
  wfout = []
  for m in mfs: # loop over manifols
    wfout += disentangle_manifold(m,A) # disentangle this manifold
  eout = [round(e,5) for e in evals]
  return eout,wfout 


def get_representation(wfs,A):
  """Gets the matrix representation of a certain operator"""
  n = len(wfs) # number of eigenfunctions
  ma = np.matrix([[0.0j for i in range(n)] for j in range(n)]) # representation of A
  from scipy.sparse import csc_matrix as csc
  sa = csc(A) # sparse matrix
  for i in range(n):
    vi = csc(np.conjugate(wfs[i])) # first wavefunction
    for j in range(n):
      vj = csc(wfs[j]).transpose() # first wavefunction
      data = (vi*sa*vj).todense()[0,0]
      ma[i,j] = data
  return ma




class Spin_chain():
  size = 0 # size of the Hamiltonian
  def __init__(self):
        self.path = os.getcwd()+"/.pychainfolder/"
        os.system("rm -rf "+self.path) # remove temporal folder
        os.system("mkdir "+self.path) # create temporal folder
        self.inipath = os.getcwd() # initial path
  def to_folder(self): os.chdir(self.path)
  def to_origin(self): os.chdir(self.inipath) # go to original folder
  def build(self,spins):
    """Creates a spin chain, using as input the spins list"""
    self.to_folder() # go to temporal folder
    if get_dimension(spins)>maxsize: 
      print("Surpased maximum allowed dimension for ED",maxsize)
      print("Dimension of the requested Hilbert space",get_dimension(spins))
      raise
    from . import chain
    ###############################
    # now read the spin operators #
    ###############################
    sobj = chain.get_chain(spins) # return a list of classes with the spins
    self.nspins = len(spins)
    self.spins = spins 
    self.sxi = sobj.sxi  # store different sx
    self.syi = sobj.syi  # store different sy
    self.szi = sobj.szi  # store different sz
    self.ski = [sobj.sxi,sobj.syi,sobj.szi]  # store by components
    self.sx = sobj.sx  # store sx
    self.sy = sobj.sy  # store sy
    self.sz = sobj.sz  # store sz
    self.wf0 = None # ground state
    self.e0 = None # ground state energy
    self.hamiltonian = None # Hamiltonian, as a multioperator
    self.sk = [sobj.sx,sobj.sy,sobj.sz]  # store by components
    self.s = self.sx*self.sx + self.sy*self.sy + self.sz*self.sz # total spin
    self.size = self.sxi[0].shape[0] # store size of the Hamiltonian
    # now read the basis
    self.basis = np.genfromtxt("basis.out").astype(int)
    self.to_origin() # go back
  def generate_hamiltonian(self,xs,ys,js,mode="operator",is_ising=False):
    """Generate a Hamiltonian"""
    h = csc(([],([],[])),shape=(self.size,self.size)) # initialize 
    if is_ising: # ising type
      for (x,y,j) in zip(xs,ys,js): # loop over couplings
        h = h + j*self.szi[x]@self.szi[y] # component
    else: # not ising model
      for (x,y,j) in zip(xs,ys,js): # loop over couplings
        h = h + j*self.sxi[x]@self.sxi[y] # component
        h = h + j*self.syi[x]@self.syi[y] # component
        h = h + j*self.szi[x]@self.szi[y] # component
    return h
  def add_heisenberg(self,xs,ys,js):
    return self.generate_hamiltonian(xs,ys,js,is_ising=False)
  def get_identity(self):
      return sparse.identity(self.sx.shape[0],dtype=np.complex)
  def gs_energy(self):
      """Compute ground state energy"""
      if self.e0 is not None: return self.e0
      if self.hamiltonian is None: raise # no Hamiltonian
      h = self.get_operator(self.hamiltonian) # generate the matrix
      (e0,wf0) = algebra.ground_state(h) # return energy
      self.wf0 = wf0 # store ground state
      self.e0 = e0 # store ground state energy
      return e0 # return energy
  def vev(self,MO):
      """Compute a certain expectation value"""
      self.gs_energy() # compute ground state energy
      op = self.get_operator(MO) # get operator
      return algebra.braket_wAw(self.wf0,op)
  def get_operator(self,name,i=0):
      """Return an operator"""
      if type(name)==multioperator.MultiOperator:
          return multioperator.MO2matrix(name,self) # return operator
      else:
        if name=="X" or name=="Sx": return self.sxi[i]
        elif name=="Y" or name=="Sy": return self.syi[i]
        elif name=="Z" or name=="Sz": return self.szi[i]
        else: raise
  def add_exchange(self,xcs,gmatrix=None): 
    """ Return matrix with the exchange field"""
    h = csc(([],([],[])),shape=(self.size,self.size)) # initialize 
    ii = 0
    for xc in xcs: # loop over couplings
      if gmatrix is None:
        h = h + xc[0]*self.sxi[ii] # add the new term
        h = h + xc[1]*self.syi[ii] # add the new term
        h = h + xc[2]*self.szi[ii] # add the new term
      else:
        for i in range(3):
          for j in range(3):
            h = h + xc[i]*gmatrix[ii][i,j]@self.ski[j][ii] # add the new term
      ii += 1 # increase counter
    return h # return zeeman term 
  def add_uniaxial_anisotropy(self,ds=None,axis=None): 
    """Add anisotropy term"""
    h = csc(([],([],[])),shape=(self.size,self.size)) # initialize 
    if ds is None: ds = [1. for i in range(self.nspins)]
    if axis is None: # assume in z direction
      axis = [[0.,0.,1.] for d in ds] # in the z direction
    axis = [np.array(r) for r in axis] # converto to array
    for i in range(len(ds)): # loop over spin sites
      r = axis[i] # get the axis
      r = r/np.sqrt(r.dot(r)) # unit vector
      h = h + ds[i]*r[0]*self.sxi[i]@self.sxi[i] # component
      h = h + ds[i]*r[1]*self.syi[i]@self.syi[i] # component
      h = h + ds[i]*r[2]*self.szi[i]@self.szi[i] # component
    return h # return anisotropy term    
  def add_dm_interaction(self,dij):
    """Adds DM interaction, favouring spin canting"""
    h = csc(([],([],[])),shape=(self.size,self.size)) # initialize 
    for i in range(len(dij)): # loop over spins
      for j in range(len(dij)): # loop over spins
        d = dij[i][j] # vector of the DM for this pair
        if d is None: continue # next iteration
        # cross product
        x = self.syi[i]*self.szi[j] - self.szi[i]*self.syi[j]
        y = self.szi[i]*self.sxi[j] - self.sxi[i]*self.szi[j]
        z = self.sxi[i]*self.syi[j] - self.syi[i]*self.sxi[j]
        v = d[0]*x + d[1]*y + d[2]*z # term in the Hamiltonian
        h = h + v # add contribution
    return h # return term in the Hamiltonian
  def add_tensor_interaction(self,t,xs=None,ys=None):
    """Generate a Hamiltonian"""
    h = csc(([],([],[])),shape=(self.size,self.size)) # initialize 
    def sij(i,j): # return a certain spin matrix i of iesim spin j
      if i==0: return self.sxi[j]
      if i==1: return self.syi[j]
      if i==2: return self.szi[j]
    if xs is None: xs = range(self.nspins)
    if ys is None: ys = range(self.nspins)
    for x in xs: # loop over couplings
      for y in ys: # loop over couplings
        if callable(t): ttmp = t(x,y) # call
        else: ttmp = t[x,y] # call
        for i in range(3): # loop over rows
          for j in range(3): # loop over columns
            h = h + ttmp[i,j]*sij(i,x)*sij(j,y) # add component
    return h
  def select_state(self,indexes):
    """Return a certain vector"""
    return states.select_state(self,indexes)
  def template(self,name="open_chain",j=1.0):
    """Create a spin chain using different templates"""
    from templates import sc_template
    return sc_template(self,name=name,j=j)


def traceless(h):
  """Turn the Hamiltonian traceless"""
  return h - np.sum(h.diagonal())/h.shape[0]*sparse.identity(h.shape[0])




