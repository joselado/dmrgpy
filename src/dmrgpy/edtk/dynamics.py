from ..algebra import algebra
import numpy as np
#from numba import jit


def dynamical_correlator(h,a0,b0,delta=2e-2,
        es=np.linspace(-1.0,3.0,300)):
    """Compute a dynamical correlator at finite temperature"""
    emu,vs = algebra.eigh(h)
    U = np.array(vs) # matrix
    Uh = np.conjugate(np.transpose(U)) # Hermitian
    A = Uh@a0@U # get the matrix elements
    B = Uh@b0@U # get the matrix elements
    out = 0.0+es*0.0*1j # initialize
    out = dynamical_sum(emu,es+1j*delta,A,B,out) # perform the summation
    return (es,out.imag/np.pi) # return correlator

#@jit(nopython=True)
def dynamical_sum(es,ws,A,B,out):
    """Return the sum giving the dynamical correlator"""
    out = out*0.0 # initialize
    es = es-np.min(es) # remove minimum
    n = len(es) # number of energies
    for iw in range(len(ws)): # loop over frequencies
        i = 0
        for j in range(n): # loop over energies
            tmp = A[i,j]*B[j,i]
            tmp *= 1./(ws[iw]+es[i] - es[j])
            out[iw] = out[iw] + tmp
    return out # return dynamical correlator


