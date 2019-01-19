from __future__ import print_function,division
from .tensorial import tensorial_LO
import numpy as np
from . import traceoverf90
from . import tensorial
from scipy.sparse import csc_matrix,kron
from scipy.sparse import coo_matrix
from . import tensorialf90
import time
from scipy.sparse import linalg as slg # linear algebra library
from scipy import linalg as lg # linear algebra library
from scipy.sparse.linalg import LinearOperator
import scipy.sparse as sp
from . import dmrg

diagmax = 100000


def tensorialoperator(op1,op2,sparse=True,LO=False):
  """Return the tendorial operator from two input operators"""
  return kron(csc_matrix(op1),csc_matrix(op2))
#  if LO: return tensorial_LO(op1,op2)
#  n1 = op1.shape[0] # dimension of the first
#  n2 = op2.shape[0] # dimension of the second
#  if sparse: # use sparse routine
#    op1 = coo_matrix(op1) # convert to coo_matrix
#    op2 = coo_matrix(op2) # convert to coo_matrix
#    # fix in case one of the operators is only zeros
#    if len(op1.col) ==0:
#      op1.col = np.array([0])
#      op1.row = np.array([0])
#      op1.data = np.array([0.])
#    if len(op2.col) ==0:
#      op2.col = np.array([0])
#      op2.row = np.array([0])
#      op2.data = np.array([0.])
#    (data3,row3,col3) = tensorialf90.sparse_tensorial_operator(op1.data,
#                     op1.row+1,op1.col+1,op2.data,op2.row+1,op2.col+1,n1,n2)
#    op3 = sp = csc_matrix((data3,(col3-1,row3-1)),shape=(n1*n2,n1*n2))
##    print("Time in tensorial operator",tnew-told)
#    op3 = csc_matrix(op3)
#    op3.eliminate_zeros() # remove zero entries
#    return op3.T # clear zero entries
#  else:
#    try: op1 = op1.todense() # dense matrices
#    except: pass
#    try: op2 = op2.todense() # dense matrices
#    except: pass
#    op3 = tensorialf90.tensorial_operator(op1,op2) # operator as tensorial  
#  return csc_matrix(op3) # return the tensorial operator



def traceover(wf,n1,n2):
  """For an input wavefunction, trace over the n2 space"""
  if len(wf) != n1*n2: raise # check that is has the correct dimensions
  dmat = traceoverf90.traceover(wf,n1,n2) # get the density matrix
  return np.matrix(dmat) # return density matrix



def getentropy(es):
  es = es[es>0.000000001] # only retain values that will not give an error
  s = -np.sum(es*np.log(es))
  return s


class GSout():
    """Class for the output of the ground state"""
    pass




def groundstate(h,dim1,dim2,v0=None,target=0,diag_states=20,
                        indict=None):
  """Get the ground state"""
  gsout = GSout() # create output
  target = indict["target"] # state to retain
  # if dictionary has been provided
  retain_states = indict["retain_states"]
  if h.shape[0]>diagmax:
    if not dmrg.silent: print("Too large dimension, ",h.shape[0])
    raise
  if v0 is not None:
    if h.shape[0] != v0.shape[0]:
      v0=None
      print("Wrong dimension of v0")
  told = time.clock() # old time
  if not dmrg.silent: print("Dimension of the Hamiltonian",h.shape)
#  if v0 is not None: diag_states = 2
  diag_states = indict["diag_states"] # state to retain
  if target>diag_states-1: diag_states = target+1 # if not enough states
  if retain_states>diag_states: diag_states = retain_states # if not enough states
  if indict["diag_mode"]=="gradient": # ground state
    X = np.random.random((h.shape[0],diag_states)) 
#+ 1j*np.random.random((h.shape[0],diag_states)) # initial vector
    (es,wfs) = slg.lobpcg(h,X=X,largest=False,tol=indict["tol"],maxiter=20000)
  if indict["diag_mode"]=="arpack": # ground state
    (es,wfs) = slg.eigsh(h,k=diag_states,which="SA",
                  v0=None,tol=indict["tol"]) # ground state
  wfs = np.transpose(wfs) # transpose
  tnew = time.clock() # old time
  if not dmrg.silent: print("Time in ground state",tnew-told)
  if len(es)>1:
    eout = sorted(es)[target] # which eigenvalue to retain
#    print("Ground state energy is",mine)
    for (e,wf) in zip(es,wfs):
      if e==eout: wfout = wf ; break
  else: # legth 1
    eout = es[0]
    wfout = wfs[0]
  ####################################
  ####################################
  if callable(indict["target_function"]):
    target_function = indict["target_function"]   
    wfsout = []
    for wfsi in getmanifolds(es,wfs): 
      wfsout += target_function(indict,wfsi) # so it is in that manifold
    wfs = wfsout # update list of waves retained
    print("Found ",len(wfsout),"useful waves") 
    wfout = wfs[indict["target_state"]] # State that will be stored
    es = [tensorial.exp_val(wfi,h) for wfi in wfs] # energies
    eout = tensorial.exp_val(wfout,h) 
    print("Energies",es-min(es))
  ####################################
  ####################################
  # calculate the density matrices
  retain_states = min(retain_states,len(wfs)) # number of waves
  dmats = [traceover(wfs[i],dim1,dim2) for i in range(retain_states)] # trace over right site + block
  # compute the "distance" between excited and GS dmat
  dmatdis = [] # empty list
  for i in range(1,len(dmats)): # loop
      A = dmats[i]*dmats[0] # define product difference
      B = dmats[0]*dmats[i] # define product difference
      C = A-B # difference
      dmatdis.append((C*C.H).trace()[0,0].real) # distance
  dmat = 0. # initialize
  for d in dmats: dmat += d # add contribution
  dmat /= retain_states # normalize
#  dmat = dmat.T # transpose
  gsout.energy = eout
  gsout.energies = sorted(es)
  gsout.wf = wfout
  gsout.dm = dmat
  gsout.dmdis = dmatdis # distance
#  print(dmatdis)
#  print(es)
  return gsout # return output



def coupledhamiltonian(ops1,ops2,cs=None,LO=False):
  """Calculate the Hamiltonian for two coupled operators"""
  id1 = sp.eye(ops1[0].shape[0],dtype=np.complex) # identity operator
  id2 = sp.eye(ops2[0].shape[0],dtype=np.complex) # identity operator
  nout = ops1[0].shape[0]*ops2[0].shape[0] # output dimension 
#  LO = False # do not use linear operators
#  if LO: # initialize a trivial linear operator 
#    def fun(v): return np.zeros(v.shape[0])
#    h = LinearOperator((nout,nout),matvec=fun) # zero linear operator
  h = csc_matrix((nout,nout),dtype=np.complex) # zero sparse matrix
#(nout,nout,dtype=np.complex) # zero sparse matrix
  nc = len(ops1) # number of operators
  if cs is None: cs = [1. for ic in range(nc)]
  for ic in range(nc): # loop over couplings, to right part
    op1 = tensorialoperator(ops1[ic],id2,LO=LO)
    op2 = tensorialoperator(id1,ops2[ic],LO=LO)
    h = h + cs[ic]*op1*op2 # coupling 1 in tensorial
  return h



def getmanifolds(es,wfs,tol=0.0001):
  """Returns a list with the wavefunction, sorted in different
  manifolds with the same energy"""
  mout = [] # list with list
  se = [es[0]] # energies
  for e in es: # loop over energies
    if np.min(np.abs(np.array(se)-e))>tol: se.append(e) # store energy
  se =sorted(se)  # sort energies
  for e in se: # loop over non degenerate energies
    mi = []
    for i in range(len(es)): # loop over waves
      if np.abs(es[i]-e)<tol:
        mi.append(wfs[i]) # store
    mout.append(mi) # store the group
  return mout # return groups





