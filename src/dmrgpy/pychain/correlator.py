from __future__ import print_function
import numpy as np
from .spectrum import viAvj
from . import spectrum
from scipy.sparse import csc_matrix as csc
from scipy.sparse import identity
import scipy.sparse.linalg as lg
import scipy.sparse.linalg as slg
from ..algebra import algebra
from ..algebra import kpm

# calculate dynamical correlators

def spismj(sc,h0,es=None,i=0,j=0,delta=0.1):
  """Calculate a correlation function S+S- in a frequency window"""
  e0,wf0 = spectrum.ground_state(h0) # get the ground state
  if es is None:
    es = np.linspace(-1.0,7.0,int(40/delta))
  sm = sc.sxi[i] - 1j*sc.syi[i] # S-
  sp = sc.sxi[j] + 1j*sc.syi[j] # S+
  iden = np.identity(h0.shape[0],dtype=np.complex) # identity
  out = []
  for e in es: # loop over energies
    g = ((iden*(e+e0+1j*delta)-h0).I - (iden*(e+e0-1j*delta)-h0).I)/2.
    op = sp*g*sm # operator
    o = viAvj(wf0,op,wf0) # correlator
    out.append(o)
  return es,-np.array(out).imag/np.pi # return result



def dynamical_correlator(sc,h0,es=None,i=0,j=0,
        delta=0.1,namei="X",namej="X",mode="full"):
  """Calculate a correlation function SiSj in a frequency window"""
  e0,wf0 = spectrum.ground_state(h0) # get the ground state
  if es is None:
    es = np.linspace(-1.0,7.0,int(40/delta))
  if namei=="X": sm = sc.sxi[i]
  elif namei=="Y": sm = sc.syi[i]
  elif namei=="Z": sm = sc.szi[i]
  else: raise
  if namej=="X": sp = sc.sxi[j]
  elif namej=="Y": sp = sc.syi[j]
  elif namej=="Z": sp = sc.szi[j]
  else: raise
  if mode=="full" and h0.shape[0]>100: mode = "cv"
  iden = identity(h0.shape[0],dtype=np.complex) # identity
  if mode=="full": iden = iden.todense() # dense matrix
  out = []
  for e in es: # loop over energies
      if mode=="full": # using exact inversion
        g = ((iden*(e+e0+1j*delta)-h0).I - (iden*(e+e0-1j*delta)-h0).I)/2.
        op = sp*g*sm # operator
        o = viAvj(wf0,op,wf0) # correlator
      elif mode=="cv": # correction vector algorithm
          o1 = solve_cv(h0,wf0,sp,sm,e+e0,delta=delta) # conjugate gradient
          o2 = solve_cv(h0,wf0,sp,sm,e+e0,delta=-delta) # conjugate gradient
          o = (o1 - o2)/2. # substract
      else: raise # not recognised
      out.append(o)
  return es,1j*np.array(out)/np.pi # return result




def solve_cv(h0,wf0,si,sj,w,delta=0.0):
     # use algorithm to solve A*x = b
     iden = identity(h0.shape[0],dtype=np.complex) # identity
     b = -delta*sj*np.matrix(wf0).T # create the b vector
     A = (h0 - w*iden)*(h0-w*iden) + iden*delta*delta # define A matrix
     b = np.array(b).reshape((b.shape[0],)) # array
     x,info = slg.cg(A,b,tol=1e-10) # solve the equation
#     print(np.max(np.abs(A*np.matrix(x).T-np.matrix(b).T)))
     x = np.matrix(x).T # column vector
     x = 1j*x + (h0 - w*iden)*x/delta # full correction vector
#     v1 = (iden*(w+1j*delta)-h0).todense().I*sj*np.matrix(wf0).T
#     print(np.max(np.abs(x-v1))) # check the difference
#     x = v1.copy() # copy
     x = si*x # apply second operator
     o = (np.matrix(wf0).H.T*x).trace()[0,0] # compute the braket
     return o





def spismj_kpm(sc,h0,es=np.linspace(-1.0,4.0,300),i=0,j=0,delta=0.1,n=500):
  """Calculate a correlation function S+S- in a frequency window"""
  e0,wf0 = lg.eigsh(-h0,k=1,ncv=20,which="LA")
  emax,wfmax = lg.eigsh(h0,k=1,ncv=20,which="LA")
  e0,wf0 = -e0[0],np.transpose(wf0)[0]
  emax = emax[0]
  sm = sc.sxi[i] - 1j*sc.syi[i] # S-
  sp = sc.sxi[j] - 1j*sc.syi[j] # S+
  vi = sm*csc(wf0).transpose() 
  vj = sp*csc(wf0).transpose() 
  h = -identity(h0.shape[0])*e0+h0
  x = es # energies
  scale = (emax-e0)*1.3
  (xs,ys) = kpm.dm_vivj_energy(h,vi,vj,scale=scale,
                                    npol=n,ne=n*10,x=x)
  return xs,ys




def dynamical_correlator_kpm(sc,h0,es=np.linspace(-1.0,4.0,300),i=0,j=0,
             n=500,namei="X",namej="X"):
  """Calculate a correlation function S+S- in a frequency window"""
  e0,wf0 = slg.eigsh(-h0,k=1,ncv=20,which="LA")
  emax,wfmax = slg.eigsh(h0,k=1,ncv=20,which="LA")
  e0,wf0 = -e0[0],np.transpose(wf0)[0]
  emax = emax[0]
  if namei=="X": sm = sc.sxi[i]
  elif namei=="Y": sm = sc.syi[i]
  elif namei=="Z": sm = sc.szi[i]
  else: raise
  if namej=="X": sp = sc.sxi[j]
  elif namej=="Y": sp = sc.syi[j]
  elif namej=="Z": sp = sc.szi[j]
  else: raise
  vi = sm*csc(wf0).transpose() 
  vj = sp*csc(wf0).transpose() 
  h = -identity(h0.shape[0])*e0+h0
  scale = (emax-e0)*1.2 # scale of the KPM method
  if es is None: es = np.linspace(-scale,scale,n*4)*0.99
  x = es[abs(es)<scale] # restrict the interval
  (xs,ys) = kpm.dm_vivj_energy(h,vi,vj,scale=scale,
                                    npol=n*4,ne=n*10,x=x)
  # now interpolate the function
  ys = ys/scale # normalize
  from scipy.interpolate import interp1d
  fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
  fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
  return es,fr(es)+1j*fi(es) # return





def static(sc,h,namei="X",namej="X",i=0,j=0,apply_hamiltonian=False):
    """Return a certain correlator"""
    wf = spectrum.ground_state(h)[1] # get wavefunction
    def getop(name,i):
        if name=="X" or name=="Sx": return sc.sxi[i]
        elif name=="Y" or name=="Sy": return sc.syi[i]
        elif name=="Z"or name=="Sz": return sc.szi[i]
        else: raise
    if apply_hamiltonian: m = getop(namei,i)@h@getop(namej,j) # get operator
    else: m = getop(namei,i)@getop(namej,j) # get operator
    return algebra.braket_wAw(wf,m,wf)
