### functions to compute the fermionic parity


from .fermionchain import Fermionic_Chain

def explicit_parity(wf):
    """Given a wavefunction, compute the parity using
    an explicit algorithm"""
#    if not isinstance(wf.MBO,Fermionic_Chain): # check
#        print("Only implemented for fermionic systems")
#        raise
    N = wf.MBO.N # density operators
    wf0 = wf.copy() # make a copy
    for Ni in N: # loop over operators
        wf = -2*Ni*wf + wf # iterate over the wavefunction
    return wf0.dot(wf) # return the parity



def fermi_string_parity(wf):
    """Given a wavefunction, compute the parity using
    an explicit algorithm"""
#    if not isinstance(wf.MBO,Fermionic_Chain): # check
#        print("Only implemented for fermionic systems")
#        raise
    F = wf.MBO.F # Fermi string operators
    wf0 = wf.copy() # make a copy
    for Fi in F: # loop over operators
        wf = Fi*wf # iterate over the wavefunction
    return wf0.dot(wf) # return the parity



