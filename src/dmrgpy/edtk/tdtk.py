import numpy as np
import scipy.linalg as lg
import scipy.sparse.linalg as slg
from scipy.sparse import csc_matrix
from scipy.sparse import identity
#from numba import jit




from scipy.integrate import solve_ivp


def evolve(w,h,t=0.0,mode="scipy",dt=0.01,de=0.0,dp=0.0):
  """Evolve the wavefunctions using the Schrodinger equation"""
  if de != 0.0 or dp != 0.0:
    return discrete_relaxed_evolution(w,h,t,dt,de=de,dp=dp)
  elif mode=="caley":
    return discrete_evolution(w,h,t,dt)
  elif mode=="taylor":
    return discrete_evolution_taylor2(w,h,t,dt)
  elif mode=="scipy":
    return scipy_evolution(w,h,t)
  else: raise

def scipy_evolution(psi,h,t):
  """e^{+iht}|psi>, exactly (to rounding) and unitarily.

  This was solve_ivp RK45 at scipy's default rtol=1e-3/atol=1e-6, whose
  error is set by t times the absolute energy of the state rather than by
  anything physical: a constant +20 added to H moved evolve_and_measure(
  mode="ED") from 2.1e-7 to 5.6e-4 off exact at dt=0.1, and at dt=0.2 on an
  8-site chain ED, the reference every ED-versus-DMRG real-time test is
  held to, was 1.2 per cent off where "python" TDVP was at 7e-9
  (2026-09-24c audit, finding 17). The sign stays e^{+iht}: evolution_DC
  is built on it, and evolution_ABC passes -h."""
  return slg.expm_multiply(1j*t*h,psi)




from scipy.sparse import issparse

# @jit
def discrete_evolution(wave,h,t=0.0,dt=0.001,order=1):
  """ Evolves a wavefunction using Caley's form"""
  if issparse(wave):
    w = np.array(wave.todense())
  else:
    w = np.array(wave) # convert into sparse matrix
  nt = np.round(np.int(t/dt)) # number of steps
  if nt == 0: # no steps
    nt = 1 # one step
  dt = t/nt # renormalize time interval steps
  iden = identity(h.shape[0],dtype=np.complex128)
  if order==1: # order of pade approximant
    u1 = iden + 1j*h*dt/2.
    u2 = iden - 1j*h*dt/2.
  for i in range(nt): # loop over steps
    b = u1*w # right side of the equation
    w = slg.spsolve(u2,b) # solve the equation
  return w




def discrete_evolution_taylor2(wave,h,t=0.0,dt=0.001):
  """ Evolves a wavefunction using Caley's form"""
  w = np.array(wave) # convert into sparse matrix
  nt = np.round(np.int(t/dt)) # number of steps
  if nt == 0: # no steps
    nt = 1 # one step
  dt = t/nt # renormalize time interval steps
  iden = identity(h.shape[0],dtype=np.complex128)
  for i in range(nt): # loop over steps
    wtmp = h*dt*w # temporal vector
    w = (iden - h*dt/2.)*wtmp + 1j*wtmp # second order formula
    w /= np.sqrt(w.dot(np.conjugate(w))) # normalize
  return w









def discrete_relaxed_evolution(wave,h,t=0.0,dt=0.001,order=1,de=0.0,dp=0.0):
  """ Evolves a wavefunction using Caley's form"""
  w = np.array(wave) # convert into sparse matrix
  nt = np.round(np.int(t/dt)) # number of steps
  if nt == 0: # no steps
    nt = 1 # one step
  dt = t/nt # renormalize time interval steps
  iden = identity(h.shape[0],dtype=np.complex128)
  if order==1: # order of pade approximant
    u1 = iden + (1j*h - de*h)*dt
    u2 = iden - (1j*h - de*h)*dt
  for i in range(nt): # loop over steps
    b = u1*w # right side of the equation
    w = slg.spsolve(u2,b) # solve the equation
    phi = dp*np.random.random(len(w))*2.*np.pi # random phases
    w *= np.exp(phi) # add random phase
    w /= np.sqrt(w.dot(np.conjugate(w))) # normalize
  return w


