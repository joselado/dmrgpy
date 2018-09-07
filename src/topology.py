import numpy as np


def berry_phase(f,nk=10):
    """Compute the Berry phase, requires a function
    that returns a manybody Hamiltonian"""
    ks = np.linspace(0.,1.,nk,endpoint=False) # kpoints
    wfs = [f(k).get_gs() for k in ks] # get the GS wavefunctions
    berry = 1.0 + 0.0j # initialize
    for i in range(nk-1): # loop
      fac = wfs[i].dot(wfs[i+1])
      berry *= fac
    fac = wfs[nk-1].dot(wfs[0]) # last one
    berry *= fac # last factor
    berry = np.arctan2(berry.imag,berry.real)/np.pi # berry phase
    return berry


