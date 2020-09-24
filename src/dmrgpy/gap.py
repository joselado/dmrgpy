import numpy as np
from .excited import gram_smith

def sector_gap(self,A,d=1.,lagrange=1.0,info=False):
    """Compute a gap in a sector that commutes with the Hamiltonian"""
    mbo = self.copy() # copy the object
    h0 = mbo.hamiltonian # get the Hamiltonian
    diff = h0*A - A*h0 # check that it commutes
#    if not mbo.is_zero_operator(diff,simplify=False,ntries=1): 
#        raise # not implemented
    # now try to diagonalize in a new sector
    wf0 = mbo.get_gs() # get the true ground state
    n0 = mbo.vev(A) # compute this expectation value
    B = A - n0 - d # operator for the Lagrange multiplier
    while True:
      h1 = h0 + lagrange*B*B # define a Hamiltonian that has a minimum at +d VEV
      mbo.set_hamiltonian(h1) # set the new Hamiltonian
      wf1 = mbo.get_gs() # get the new ground state
#      ws = [wf0,wf1] # set of wavefunctions
#      ws = gram_smith(ws) # orthogonalize
#      wf1 = ws[1]
      n1 = mbo.vev(A) # compute this expectation value
      info = True
      if abs(n1-n0-d)<1e-1: break # break the loop
      else:
          if info: print("Not found yet, recomputing",n0,n1,d)
      lagrange *= 2 # twice as big
    e0 = wf0.overlap(h0*wf0) # ground state energy
    e1 = wf1.overlap(h1*wf1) # ground state energy
    return (e1-e0).real # return the energy difference


