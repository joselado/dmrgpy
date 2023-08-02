from __future__ import print_function
import numpy as np

import scipy.linalg as lg
import scipy.optimize as optimize

try:
  import density_matrixf90 as dm90
  use_fortran = True
except:
  print("PROBLEM WITH FORTRAN COMPILATION, ENTROPY IS NOT CALCULATED")
  use_fortran = False


def reduced_density_matrix(wave,basis,site=0,use_fortran=True):
  """Calculates the reduced density matrix of a certain wave, projecting
  onto the subspace site"""
  dimdm = np.max(basis[:,site])+1 # dimension of density matrix
  if dimdm==len(wave): raise # if no partition
  # the basis has to be transposed
  if use_fortran:
    dmat = dm90.reduced_density_matrix(wave,basis.T+1,site+1,dimdm)
  else:
    print("PROBLEM WITH FORTRAN COMPILATION, ENTROPY IS NOT CALCULATED")
    dmat = np.zeros((dimdm,dimdm),dtype=np.complex_)
    dmat[0,0] = 1.
  return dmat



def entropy_matrix(m):
  """Return the entropy of a matrix"""
  if np.abs(m.trace()-1.0)>0.0001: 
    print(m.trace())
    raise
  eig = lg.eigvalsh(m)
  s = 0.0
  for e in eig: # loop over eigenvalues
    if e>0.00001: s += e*np.log(e) # add contribution to the entropy
  return -s




def entropy(wave,basis,site=0):
  """Calculate the entanglement entropy of a certain site"""
  dmat = reduced_density_matrix(wave,basis,site=site)
  dmat = np.matrix(dmat) # convert to matrix
  if (np.sum(np.abs(dmat-dmat.H)))> 0.0001: raise # check is hermitian
  return entropy_matrix(dmat)




def thermal_entropy(waves,energies,basis,site=0,
                                    temps=np.linspace(0.001,1.,100)):
  """Thermal density matrix"""
  dmats = [] # list with all the density matrix
  for wave in waves: # loop over states
    dmat = reduced_density_matrix(wave,basis,site=site)
    dmats.append(dmat) # store density matrix
  # now do the loop over temepratures
  ss = []
  for it in temps: # loop over t
    dmat = dmats[0]*0. # initialize
    z = 0.0 # partition function
    for (d,energy) in zip(dmats,energies): # loop over energies
      zfac = np.exp(-energy/it) # boltzman factor
      z += zfac # add to partition function
      dmat += zfac*d # add to density matrix
    dmat /= z # normalize
    s = entropy_matrix(dmat) # store entropy
#    print(s,zfac)
    ss.append(s) # store entropy
  return ss


def minimize_entropy(waves,basis,site=0):
  """Minimize the entanglement entropy, and return both the entropy and
  wavefunction of the state"""
  def get_wf(cs):
    w0 = waves[0]*cs[0] # initialize wavefunction
    for i in range(len(waves)-1): # loop over wavefunctions
      w0 += (cs[2*i+1] + 1j*cs[2*i+2])*waves[i+1]
    w0 = w0/np.sqrt(w0.dot(np.conjugate(w0))) # normalize wavefunction
    return w0 # return wavefunction
  def fun(cs):
    """Function to optimize"""
    return entropy(get_wf(cs),basis) # return entropy
  x0 = np.random.random(2*len(waves)-1) # initialize
  bounds = [(-1.,1.) for x in x0]
  optresult = optimize.minimize(fun,x0,method="SLSQP",bounds=bounds)
  smin = fun(optresult.x) # minimum entropy
  wfmin = get_wf(optresult.x) # minimum wavefunction
  return smin








def basis_group(basis,sites=[0]):
  """Generates a new basis, where sites are grouped in a single site"""
  bt = basis.T # transpose basis




