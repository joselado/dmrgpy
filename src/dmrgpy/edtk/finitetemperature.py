# routines for exact diagonalization at finite temperature
from ..algebra import algebra
import numpy as np

try:
  from numba import jit
except:
  print("Numba is not present")
  def jit(*args,**kwargs): # dummy
    return *args

def thermal_rho(h,beta=1.0):
    """Return the thermal density matrix"""
    out = algebra.expm(-beta*h) # return density matrix
    out = out/np.trace(out) # normlize
    return out


def gs_energy(h,beta=1.0):
    """Compute the ground state energy"""
#    return np.trace(h@thermal_rho(h,beta=beta)).real
    es = algebra.eigvalsh(h) # eigenvalues
    ws = np.exp(-beta*es) # statistical weights
    ws = ws/np.sum(ws) # normalize statistical weights
    return np.sum(es*ws) # return weighted average
    
def measure(h,ops,beta=1.0):
    """Measure certain operators"""
    rho = thermal_rho(h,beta=beta) # density matrix
    return np.array([np.trace(op@rho) for op in ops]) # return values


def dynamical_correlator(h,a0,b0,delta=2e-2,
        T=1.0,es=np.linspace(-1.0,3.0,300)):
    """Compute a dynamical correlator at finite temperature"""
    emu,vs = algebra.eigh(h) 
    U = np.array(vs) # matrix
    Uh = np.conjugate(np.transpose(U)) # Hermitian
    A = Uh@a0@U # get the matrix elements
    B = Uh@b0@U # get the matrix elements
    out = 0.0+es*0.0*1j # initialize
    out = dynamical_sum(emu,1./T,es+1j*delta,A,B,out) # perform the summation
    return (es,out.imag/np.pi) # return correlator

@jit(nopython=True)
def dynamical_sum(es,beta,ws,A,B,out):
    """Return the sum giving the dynamical correlator"""
    out = out*0.0 # initialize
    es = es-np.min(es) # remove minimum
    expes = np.exp(-beta*es) # boltzman factors
    n = len(es) # number of energies
    for iw in range(len(ws)): # loop over frequencies
      for i in range(n): # loop over energies
        for j in range(n): # loop over energies
            tmp = expes[i]*A[i,j]*B[j,i]
            tmp *= 1./(ws[iw]+es[i] - es[j])
            out[iw] = out[iw] + tmp
    return out/np.sum(expes) # return dynamical correlator



