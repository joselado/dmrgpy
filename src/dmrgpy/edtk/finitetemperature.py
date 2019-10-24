# routines to spin chains at finite temperature
from ..algebra import algebra

def thermal_rho(h,beta=1.0):
    """Return the thermal density matrix"""
    return algebra.expm(-beta*m) # return density matrix


def gs_energy(h,beta=1.0):
    """Compute the ground state energy"""
    es = algebra.eigvalsh(h) # eigenvalues
    ws = np.exp(-beta*es) # statistical weights
    ws = ws/np.sum(ws) # normalize statistical weights
    return np.sum(es*ws) # return weighted average
    


