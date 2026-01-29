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




def fermi_string_parity_iterative_collect(wf):
    """Given a wavefunction, compute the parity using
    an explicit algorithm using Fermi strings"""
    ncut = 40 # perform as many products
    MBO = wf.MBO
    F = wf.MBO.F # Fermi string operators
    wf0 = wf.copy() # make a copy
    ii = 0 # initialize counter
    O = MBO.Id # initialize
    Ot = MBO.Id # initialize
    Ot = MBO.toMPO(Ot) # transform to MPO
    for Fi in F: # loop over operators
        O = Fi*O # iterate
        ii += 1 # increase counter
        if ii==ncut: # maximum number reached
            O = MBO.toMPO(O) # transform to MPO
            Ot = Ot*O # iterate
#            wf = O*wf # apply to wavefunction
            ii = 0 # start over
            O = MBO.Id # start over
#        wf = Fi*wf # iterate over the wavefunction
    return wf.aMb(Ot,wf) # return the parity




def fermi_string_parity_iterative(wf):
    """Given a wavefunction, compute the parity using
    an explicit algorithm using Fermi strings"""
    F = wf.MBO.F # Fermi string operators
    wf0 = wf.copy() # make a copy
    for Fi in F: # loop over operators
        wf = Fi*wf # iterate over the wavefunction
    return wf0.dot(wf) # return the parity



def fermi_string_parity(wf):
    """Create the full parity operator"""
    F = wf.MBO.F # Fermi string operators
    O = 1. # initialize
    for Fi in F: # loop over operators
        O = Fi*O # iterate
    return wf.aMb(O,wf) # return the parity



def get_fermionic_parity(wf,fpmode="full",**kwargs):
    if fpmode=="full":
        return fermi_string_parity(wf) # parity of the state
    elif fpmode=="iterative":
        return fermi_string_parity_iterative(wf) # parity of the state
    else: raise



