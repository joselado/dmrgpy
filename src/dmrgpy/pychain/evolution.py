import numpy as np
import scipy.linalg as lg
import scipy.sparse.linalg as slg
from scipy.sparse import csc_matrix
from scipy.sparse import identity
#from numba import jit


name_sx = "SX.OUT"
name_sy = "SY.OUT"
name_sz = "SZ.OUT"

def evolve(waves,h,t=0.0,mode="caley",dt=0.01,de=0.0,dp=0.0):
  """Evolve the wavefunctions using the Schrodinger equation"""
  if de != 0.0 or dp != 0.0:
    return [discrete_relaxed_evolution(w,h,t,dt,de=de,dp=dp) for w in waves]
  elif mode=="caley":
    return [discrete_evolution(w,h,t,dt) for w in waves]
  elif mode=="taylor":
    return [discrete_evolution_taylor2(w,h,t,dt) for w in waves]
#  elif mode=="full":
#    u = lg.expm(1j*h*t) # evolution operator
#    wout = []
#    for wave in waves: # loop over wavefunctions
#      w = u*csc_matrix([wave]).T # evolve
#      w = np.array(w.todense().reshape(w.shape[0]))[0] # convert to array
#      wout.append(w) # store wavefunction
#    return wout # return waves
#  else: raise


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




def record_evolution(sc,wave,h,t=0.0,nt=100,dt=0.01,de=0.0,dp=0.0,
                      name_e="ENERGY.OUT",name_s="ENTROPY.OUT",
                      name_m = "MAGNETIZATION.OUT"):
  """Saves in a file the evolution of the expectation values
  of the spin chain"""
  import entanglement
  import build
  ff = "{:10.4f}".format
  fe = open(name_e,"w") # file to save the evolution
  fs = open(name_s,"w") # file to save the evolution
  fm = open(name_m,"w") # file to save the evolution
  fx = open(name_sx,"w") # file to save the evolution
  fy = open(name_sy,"w") # file to save the evolution
  fz = open(name_sz,"w") # file to save the evolution
  w = wave.copy() # copy the wavefunction
  tij = float(t)/nt # time to wait to write
  for it in range(nt): # loop over time steps
    w = discrete_relaxed_evolution(w,h,tij,dt,de=de,dp=dp) # new wavefunction
    # calculate expectation values
    e = build.exp_val(w,h) # energy
    fe.write(ff(tij*it)+"    "+ff(e)+"\n") # save energy
    fs.write(ff(tij*it)+"    ") # save time
    fm.write(ff(tij*it)+"    ") # save time
    fx.write(ff(tij*it)+"    ") # save time
    fy.write(ff(tij*it)+"    ") # save time
    fz.write(ff(tij*it)+"    ") # save time
    # loop over spins
    for i in range(sc.nspins): # loop over spins
      entropy = entanglement.entropy(w,sc.basis,site=i)
      sx = build.exp_val(w,sc.sxi[i])
      sy = build.exp_val(w,sc.syi[i])
      sz = build.exp_val(w,sc.szi[i])
      fs.write(ff(entropy)+"   ")
      fm.write(ff(sx)+"   ")
      fx.write(ff(sx)+"   ")
      fm.write(ff(sy)+"   ")
      fy.write(ff(sy)+"   ")
      fm.write(ff(sz)+"   ")
      fz.write(ff(sz)+"   ")
    fs.write("\n") # next line
    fm.write("\n") # next line
    fx.write("\n") # next line
    fy.write("\n") # next line
    fz.write("\n") # next line
  fe.close()
  fm.close()
  fx.close()
  fy.close()
  fz.close()
  fs.close()
  return w # return wavefunction

