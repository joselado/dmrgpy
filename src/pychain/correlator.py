from __future__ import print_function
import numpy as np
from spectrum import viAvj
import spectrum
from scipy.sparse import csc_matrix as csc
from scipy.sparse import identity
import scipy.sparse.linalg as lg

# calculate dynamical correlators

def spismj(sc,h0,es=None,i=0,j=0,delta=0.1):
  """Calculate a correlation function S+S- in a frequency window"""
  e0,wf0 = spectrum.ground_state(h0) # get the ground state
  if es is None:
    es = np.linspace(-1.0,7.0,int(40/delta))
  sm = sc.sxi[i] + 1j*sc.syi[i] # S+
  sp = sc.sxi[j] - 1j*sc.syi[j] # S-
  iden = np.identity(h0.shape[0],dtype=np.complex) # identity
  out = []
  for e in es: # loop over energies
    op = sp*(iden*(e+e0+1j*delta)-h0).I*sm # operator
#    op = (iden*(e+e0+1j*delta)-h0).I # operator
    o = viAvj(wf0,op,wf0) # correlator
#    o = op.trace()[0,0]
    out.append(o)
  return es,-np.array(out).imag/np.pi # return result


def spismj_kpm(sc,h0,es=np.linspace(-1.0,4.0,300),i=0,j=0,delta=0.1,n=500):
  """Calculate a correlation function S+S- in a frequency window"""
  e0,wf0 = lg.eigsh(-h0,k=1,ncv=20,which="LA")
  emax,wfmax = lg.eigsh(h0,k=1,ncv=20,which="LA")
  e0,wf0 = -e0[0],np.transpose(wf0)[0]
  emax = emax[0]
  sm = sc.sxi[i] + 1j*sc.syi[i] # S+
  sp = sc.sxi[j] - 1j*sc.syi[j] # S-
  vi = sm*csc(wf0).transpose() 
  h = -identity(h0.shape[0])*e0+h0
  import kpm
  x = es # energies
  scale = (emax-e0)*1.3
  (xs,ys) = kpm.dm_vivj_energy(h,vi,vi,scale=scale,
                                    npol=n,ne=n*10,x=x)
  return xs,ys/scale




