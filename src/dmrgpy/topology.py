import numpy as np


def berry_phase(f,nk=100):
    """Compute the Berry phase, requires a function
    that returns a manybody Hamiltonian"""
    ks = np.linspace(0.,1.,nk,endpoint=False) # kpoints
    f(0).get_gs() # first iteration
    wfs = [] # empty list
    for k in ks:
        sc = f(k) # get object
        sc.sites_from_file = True # use the old sites
        wfs.append(sc.get_gs()) # get wavefunction
    berry = 1.0 + 0.0j # initialize
    for i in range(nk-1): # loop
      fac = wfs[i].dot(wfs[i+1])
      print(np.abs(fac))
#      if np.abs(fac)>0.9: 
      berry *= fac/np.abs(fac)
    fac = wfs[nk-1].dot(wfs[0]) # last one
    berry *= fac # last factor
    berry = np.arctan2(berry.imag,berry.real)/np.pi # berry phase
    return berry


