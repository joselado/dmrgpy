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
        e0 = sc.gs_energy() # ground state energy
        if e0<emin: wf0 = sc.wf0 # copy wavefunction
    sc.set_initial_wf(wf0) # set the wavefunction

def gs_energy_bestof(self,**kwargs):
    maxm = self.maxm
    maxm2 = maxm
    self.maxm = maxm2
    gs_energy_many(self,**kwargs)
    self.maxm = maxm
    return gs_energy_single(self,reconverge=True)

def gs_energy_many(self,n=20,**kwargs):
    """Compute many ground states, and retain only the best one"""
    wfs = [] # list with GS
    es = [] # list with energies
    emin = 1e8 # grund state energy
    wf0 = None
    for i in range(n): # loop
        self.computed_gs = False # initialize
        self.gs_from_file = False
        self.skip_dmrg_gs = False
        e0 = gs_energy_single(self,reconverge=False,
                **kwargs) # ground state energy
        if e0<emin: 
            wf0 = self.wf0.copy() # copy wavefunction
            emin = e0
        else: self.wf0.clean() # remove wavefunction
        print("Best of",n,i,e0,emin)
    wf0.rename("psi_GS.mps")
    self.set_initial_wf(wf0) # set the wavefunction
    print("Final energy",self.vev(self.hamiltonian).real)
    return emin

def gs_energy_single(self,wf0=None,reconverge=None):
    """
    Return the ground state energy
    """
    if wf0 is not None: 
        self.set_initial_wf(wf0) # set the initial wavefunction
    if reconverge is not None: # overwrite skip_dmrg_gs
        self.skip_dmrg_gs = not reconverge # if the computation should be rerun
    self.execute(lambda: self.setup_task("GS"))
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.run() # perform the calculation
    self.wf0 = mps.MPS(self,name="psi_GS.mps")#.copy() # set the ground state
    # get the ground state energy
    out = self.execute(lambda: np.genfromtxt("GS_ENERGY.OUT"))
    self.e0 = out # store ground state energy
    self.computed_gs = True
    self.sites_from_file = True
    self.gs_from_file = True
    self.skip_dmrg_gs = True
    self.set_initial_wf(self.wf0) # set the initial wavefunction
    return out # return energy

def gs_energy(self,policy="single",**kwargs):
    if policy=="single":
        return gs_energy_single(self,**kwargs)
    if policy=="many":
        return gs_energy_many(self,**kwargs)
    else: raise


