from __future__ import print_function
import scipy.sparse.linalg as slg 
import scipy.sparse as sp
from scipy.sparse import csc_matrix as csc
import scipy.linalg as lg 
import numpy as np
#from numba import jit

nfull = 2000 # dimension for using full diagonalizetion
maxfull = 10000 # hard limit using full diagonalizetion
nprec = 9 # number of digits

def degeneracy(h):
  """Get degeneracy of ground state"""
  eig = lg.eigvalsh(h.todense())
#  eig,eigvec = slg.eigsh(h,k=1,which="SA",maxiter=100000)
#  h0 = h - eig*sp.eye(h.shape[0]) # shift to zero the matrix
#  print(eig)
  deg = np.sum(np.abs(eig-min(eig))<0.0001)
  print(deg)
  return deg


def gap(h):
  """Return the gap"""
  es = eigenstates(h,k=3)
  return es[1]-es[0]


def eigenstates(h,operator=None,k=None,evals=False):
  """Calculate spectrum of the system"""
  if k is not None:
    if k<1: k = 1
    if k>h.shape[0]-1: k = None # if states wanted bigger than size
#  if h.shape[0]<4000: k = None # small matrices
#  print(k)
  if k is None:
    if h.shape[0]>maxfull: raise
    if not evals:
      eig = lg.eigvalsh(h.todense())
      return np.sort(eig)
    else:
      eig,eigvec = lg.eigh(h.todense())
      eigvec = np.conjugate(np.transpose(eigvec)) # arrange
      return eig,eigvec
  else:
    eig,eigvec = slg.eigsh(h,k=k,which="SA",maxiter=100000)
    if not evals: return np.sort(eig)
    else: 
      eigvec = np.conjugate(np.transpose(eigvec)) # arrange
      return eig,eigvec

def ground_state(h,nmax=nfull):
  """Get a ground state"""
  if h.shape[0]>nmax:  
    print("Calling ARPACK")
    eig,eigvec = slg.eigsh(h,k=10,which="SA",maxiter=100000)
  else:  
    print("Full diagonalization")
    eig,eigvec = lg.eigh(h.todense())
  return eig[0],eigvec.transpose()[0]



def ground_states(h,k=None,de=0.00001):
  """Return the manifold of ground states"""
  if k is None: k = 10 # first iteration
  if k == 0 : k=10 # at least ten
#  k = None
#  eig,evecs = slg.eigsh(h,k=k,which="SA",maxiter=100000)
  eig,evecs = eigenstates(h,k=k,evals=True) # get eigenvalues and eigenvectors
#  print(eig)
  mine = np.min(eig) # minimum eigenvalue 
  gs = [] # empty list
  for i in range(len(eig)):
    if np.abs(eig[i]-mine)<de: gs.append(evecs[i]) # store
  if len(eig)/len(gs)<2:
    print("Not enought states, calling again",k,len(gs),mine)
    return ground_states(h,k=2*k,de=de) # recursive
  return mine,gs # return ground state waves

def correlation_wave(sc,wave,dm=0.001):
  """Calculates the correlation matrix of a spin chain"""
  ns = len(sc.sxi) # number of sites
  m = np.matrix(np.array([[0. for i in range(ns)] for j in range(ns)])) 
  for i in range(ns): # loop over site i
    for j in range(ns): # loop over site j
      # build the operator
      sisj = sc.sxi[i]*sc.sxi[j] + sc.syi[i]*sc.syi[j] + sc.szi[i]*sc.szi[j]
      # calculate the expectation value      
      m[i,j] = exp_val(wave,sisj)    
  return m

def exp_val(v,A):
  """Calculate a certain expectation value"""
  vi = csc(np.conjugate(v)) # wavefunction
  vj = csc(v).transpose() # wavefunction
  data = (vi*A*vj).todense()[0,0]
  data = round(data.real,nprec)
  return data # return the value


def viAvj(vi,A,vj):
  """Calculate a certain expectation value"""
  vi = csc(np.conjugate(vi)) # wavefunction
  vj = csc(vj).transpose() # wavefunction
  data = (vi*A*vj).trace()[0,0]
  return data # return the value









def write_correlation(sc,waves):
  """Writtes the different eigenvector of the correlation matrix"""
  try: waves[0][0] # check that it is a list
  except: waves = [waves] ; print("Converting to list in write_correlation")
  ns = len(sc.sxi) # number of sites
  m = np.matrix(np.array([[0. for i in range(ns)] for j in range(ns)])) 
  for wave in waves: # loop over waves
    m = m + correlation_wave(sc,wave) # get correlation matrix
  m = m/len(waves) # normalize by the number of waves
  # now diagonalize
  eig,evecs = lg.eigh(m) # eigenvectors and eigenvalues
  evecs = np.conjugate(np.transpose(evecs)) # transpose
  ie = 0
  for (e,v) in zip(eig,evecs): # loop over waves
    fo = open("CORRELATION_"+str(ie)+"_"+str(e)+".OUT","w")
    for i in range(len(v)):  fo.write(str(i)+"    "+str(v.real[i])+"\n")
    ie += 1 # increase counter
    fo.close()
  print(m,eig)
  return eig # return eigenvalues


def write_splitting(waves,ops):
  """Writtes in different files the different splittings created
  by splitting operators"""
  fo = open("SPLITTINGS.OUT","w") 
  iop = 0 # initialize
  for op in ops: # loop over operators
    ns = len(waves) # number of sites
    m = np.matrix(np.array([[0. for i in range(ns)] for j in range(ns)])) 
    for i in range(ns): # loop over waves
      for j in range(ns):  # loop over waves
        m[i,j] = viAvj(waves[i],op,waves[j]) # matrix element
    eig = lg.eigvalsh(m) # eigenvalues       
    fo.write(str(iop)+"  ") 
    for e in eig: fo.write(str(e.real)+"  ") # write ener
    fo.write("\n") # next line
    iop += 1 # increase counter
  fo .close() # close file  



def exp_val(v,A):
  """Calculate a certain expectation value"""
  vi = csc(np.conjugate(v)) # wavefunction
  vj = csc(v).transpose() # wavefunction
  data = (vi*A*vj).todense()[0,0]
#  data = round(data.real,5)
  return data.real # return the value




def variational(h):
  from scipy.optimize import minimize
  v0 = np.random.random(h.shape[0]) # random vector
  def fun(v):
    return exp_val(v,h)
  res = minimize(fun,v0,method='CG')
  raise
  return 0.,res.x


def full_spectra(h):
  if h.shape[0]>2000: raise
  eig = lg.eigvalsh(h.todense())
  np.savetxt("FULL_SPECTRA.OUT",np.matrix([range(len(eig)),eig]).T)
  return eig

