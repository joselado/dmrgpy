from . import mps
import numpy as np


def best_gs(sc,n=1):
    """Compute many ground states, and retain only the best one"""
    wfs = [] # list with GS
    es = [] # list with energies
    emin = 1e8 # grund state energy
    wf0 = None
    for i in range(n): # loop
        sc.computed_gs = False # initialize
        e0 = sc.gs_energy() # gounrd state energy
        if e0<emin: wf0 = sc.wf0 # copy wavefunction
    sc.set_initial_wf(wf0) # set the wavefunction


def gs_energy(self,wf0=None):
      self.to_folder() # go to temporal folder
      self.set_initial_wf(wf0) # set the initial wavefunction
      self.setup_sweep()
      self.setup_task("GS")
      self.write_hamiltonian() # write the Hamiltonian to a file
      self.run() # perform the calculation
      self.wf0 = mps.MPS(self).copy() # set the ground state
      out = np.genfromtxt("GS_ENERGY.OUT") # return the ground state energy
      self.e0 = out # store ground state energy
      self.computed_gs = True
      self.to_origin() # go to temporal folder
      return out # return energy


