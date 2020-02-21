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
  def f(t,psi):
      return 1j*h@psi # evolution
  tspan = [0.,t] # times
  t_eval = [tspan[1]]
  v0 = psi # initial wavefunction
  sol = solve_ivp(f,tspan,v0,method="RK45",t_eval=t_eval)
  return sol.y[:,0]




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
  iden = identity(h.shape[0],dtype=np.complex)
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
  iden = identity(h.shape[0],dtype=np.complex)
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
  iden = identity(h.shape[0],dtype=np.complex)
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


