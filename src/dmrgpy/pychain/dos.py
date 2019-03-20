
import numpy as np
import scipy.sparse.linalg as lg
import scipy.sparse.linalg as slg

def dos_kpm(h0,delta=1e-1):
  """Compute full DOS"""
  e0,wf0 = slg.eigsh(-h0,k=1,ncv=20,which="LA")
  emax,wfmax = slg.eigsh(h0,k=1,ncv=20,which="LA")
  e0,wf0 = -e0[0],np.transpose(wf0)[0]
  emax = emax[0]
  scale = (emax-e0)*1.2 # scale of the KPM method
  from ..algebra import kpm
  npol = int(scale/delta) # number of polynomials
  return kpm.tdos(h0,scale)


