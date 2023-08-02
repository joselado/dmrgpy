from __future__ import print_function
from . import build
import numpy as np
import scipy.linalg as lg
from . import tensorial
from . import dmrg



def biladder(s1,s2,j1=1.0,j2=1.0,j12=1.0):
  """Create a ladder"""
  spins = [s1,s2]
  sc = build.Spin_chain()
  sc.build(spins)
  c1 = np.sqrt(j1+0j) # coupling
  c2 = np.sqrt(j2+0j) # coupling
  c12 = np.sqrt(j12+0j) # coupling
  def left(i):
    return [c1*sc.sxi[0],c1*sc.syi[0],c1*sc.szi[0],c2*sc.sxi[1],c2*sc.syi[1],c2*sc.szi[1]]
  def right(i):
    return [c1*sc.sxi[0],c1*sc.syi[0],c1*sc.szi[0],c2*sc.sxi[1],c2*sc.syi[1],c2*sc.szi[1]]
  def onsite(i):
    return j12*(sc.sxi[0]*sc.sxi[1] + sc.syi[0]*sc.syi[1] + sc.szi[0]*sc.szi[1])
  odict = dict() # dictionary
  odict["right"] = right
  odict["left"] = left
  odict["onsite"] = onsite
  odict["site_operator_generator"] = right
  dmrgdict = dmrg.dmrgdict() # rest of parameters
  for key in dmrgdict: odict[key] = dmrgdict[key] # copy dictionary
  return odict


def monochain(s1,d=0.0,b=[0.,0.,0.],fun=None,J=[1.,1.,1.]):
  """
  Create a chain
   - b is the magnetic field
   - J is the different components of the exchange
   - d is the single ion anisotropy
  """
  if fun is None: fun = lambda i: 1.
  spins = [s1]
  sc = build.Spin_chain()
  sc.build(spins)
  def left(i):
    return [sc.sxi[0],sc.syi[0],sc.szi[0]]
  def right(i):
    return [sc.sxi[0],sc.syi[0],sc.szi[0]]
  def onsite(i):
    if callable(d): dtmp = d(i)
    else: dtmp = d
    if callable(b): raise # bi = b(i) this does not work
    else: bi = b
    return dtmp*(sc.szi[0]*sc.szi[0]) + (bi[0]*sc.sxi[0] +
         bi[1]*sc.syi[0] + bi[2]*sc.szi[0])
  odict = dict() # dictionary
  odict["right"] = right
  odict["left"] = left
  odict["onsite"] = onsite
  odict["site_operator_generator"] = left
  def f(i,j): return J # isotropic by default
  odict["coupling"] = f # assign function
  dmrgdict = dmrg.dmrgdict() # rest of parameters
  for key in dmrgdict: odict[key] = dmrgdict[key] # copy dictionary
  return odict



def frusladder(s1,j=1.0,jf=1.0):
  """Create a frustrated ladder"""
  spins = [s1,s1]
  sc = build.Spin_chain()
  sc.build(spins)
  cf = np.sqrt(jf+0j) # coupling
  c = np.sqrt(j+0j) # coupling
  def left(i):
    out1 = [c*sc.sxi[0],c*sc.syi[0],c*sc.szi[0]]
    out1 += [c*sc.sxi[1],c*sc.syi[1],c*sc.szi[1]]
    out1 += [cf*sc.sxi[1],cf*sc.syi[1],cf*sc.szi[1]]
    out1 += [cf*sc.sxi[0],cf*sc.syi[0],cf*sc.szi[0]]
    return out1
  def right(i):
    out2 = [c*sc.sxi[0],c*sc.syi[0],c*sc.szi[0]]
    out2 += [c*sc.sxi[1],c*sc.syi[1],c*sc.szi[1]]
    out2 += [cf*sc.sxi[0],cf*sc.syi[0],cf*sc.szi[0]]
    out2 += [cf*sc.sxi[1],cf*sc.syi[1],cf*sc.szi[1]]
    return out2
  def onsite(i):
    return c*(sc.sxi[0]*sc.sxi[1] + sc.syi[0]*sc.syi[1] + sc.szi[0]*sc.szi[1])
  odict = dict() # dictionary
  odict["right"] = right
  odict["left"] = left
  odict["onsite"] = onsite
  odict["site_operator_generator"] = right
  odict["left_coupling"] = 1.0
  odict["right_coupling"] = 1.0
  dmrgdict = dmrg.dmrgdict() # rest of parameters
  for key in dmrgdict: odict[key] = dmrgdict[key] # copy dictionary
  return odict


def get_representation(A,wfs):
  """Get representation of an operator"""
  nw = len(wfs) # number of waves
  m = np.matrix(np.zeros((nw,nw),dtype=np.complex_))
  for i in range(nw):
    for j in range(nw):
      m[j,i] = np.conjugate(wfs[i]).dot(A*wfs[j])
  return m # return the matrix


def target_total_s(s=None,sz=None,inds=None):
  """Returns a function that will target a state a a certain total S,
  and a certain Sz"""
  def tfun(indict,wfs):
    """Function that decides which wavefunction to target"""
    # create the total spin operator
    sis = [None,None,None] # total spin operators
    for i in range(3): # loop over components
      si = indict["left_site_operators_full"][i]
      si = si + indict["right_site_operators_full"][i]
      si = si + indict["sum_right_block_operators_full"][i]
      sis[i] = si + indict["sum_left_block_operators_full"][i] # store
    # matrix with spin operators has been built
    ssop = sis[0]*sis[0] + sis[1]*sis[1] + sis[2]*sis[2] # total S
    szop = sis[2] # Sz
    # get the eigenvectors
    if s is not None: # target a certain s
      ssm = get_representation(ssop,wfs) # get the matrix
      try: wfs = retain_sector(ssm,s*(s+1),wfs) # filter by total S
      except: return []
      print("Eigenvalues SS",np.round(lg.eigvalsh(ssm),2))
    if sz is not None: # target a certain sz
      szm = get_representation(szop,wfs) # get the matrix
      try: wfs = retain_sector(szm,sz,wfs) # filter by total S
      except: return []
      print("Eigenvalues SZ",np.round(lg.eigvalsh(szm),2))
    # sort wavefunctions by energy
#    h = indict["hamiltonian"]
#    es = [tensorial.exp_val(wfi,h) for wfi in wfs] # energies
#    wfs = [wfi for (e,wfi) in sorted(zip(es,wfs))]
    # and retain only certain indexes
#    if inds is not None: wfs = [wfs[i] for i in inds]
    return wfs
#    exit()
  return tfun # return the function


def retain_sector(M,val,wfs,tol=0.1):
  """Return wavefunctions that share a certain eigenvalue in a
    subspace"""
  (ls,R) = lg.eigh(M) # eigenvectors and eigenvalues
  R = R.T # transpose()
  inds = [] # indexes
  for (i,l) in zip(range(len(ls)),ls):
    if np.abs(l-val)<tol: inds.append(i) # store index 
  if len(inds)==0:
#    print("WARNING, sector not found")
    raise
    return []
  Ro = np.matrix([R[i] for i in inds]).H.T # build new matrix
  wfsout = [] # empty list
  for i in range(Ro.shape[0]): # loop over out waves
    wfi = wfs[0]*0. # initialize
    for j in range(Ro.shape[1]): # loop
      wfi += Ro[i,j]*wfs[j]
    wfsout.append(wfi) # store
  return wfsout # return waves




