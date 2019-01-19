
from .algebra import Av
from .evolution import discrete_evolution


def otoc(w0,h,A,B,t,dt=1e-4):
    """Compute the out of time ordered correlator"""
    def evolt(w,t): # evolve a wavefunction
        return discrete_evolution(w,h,t,dt)
    # compute all the expectation values


