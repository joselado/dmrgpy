from . import states
from ..pychain.spectrum import ground_state

def gs_energy(m0):
    """Compute ground state energy"""
    h = states.one2many(m0) # single to many body Hamiltonian
    return ground_state(h)[0] # return GS energy
