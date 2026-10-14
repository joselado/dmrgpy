import numpy as np
from .excited import gram_smith

def sector_gap(self,A,d=1.,lagrange=1.0,info=False,mode="DMRG"):
    """Compute a gap in a sector that commutes with the Hamiltonian:
    E0(<A>+d) - E0(<A>), the lowest energy of h0 + lagrange*(A-<A>-d)^2
    minus the ground-state energy of h0. mode="DMRG"|"ED" picks the solver
    of both solves (it used to be refused with a TypeError)."""
    mbo = self.copy() # copy the object
    h0 = mbo.hamiltonian # get the Hamiltonian
    # now try to diagonalize in a new sector
    e0 = mbo.gs_energy(mode=mode) # ground state energy
    n0 = mbo.vev(A,mode=mode).real # expectation value in the ground state
    B = A - n0 - d # operator for the Lagrange multiplier
    while True:
      h1 = h0 + lagrange*B*B # define a Hamiltonian that has a minimum at +d VEV
      mbo.set_hamiltonian(h1) # set the new Hamiltonian
      e1 = mbo.gs_energy(mode=mode) # lowest energy with the penalty
      n1 = mbo.vev(A,mode=mode).real # expectation value in that state
      if abs(n1-n0-d)<1e-1: break # break the loop
      else:
          if info: print("Not found yet, recomputing",n0,n1,d)
      lagrange *= 2 # twice as big
    return float(np.real(e1-e0)) # return the energy difference
