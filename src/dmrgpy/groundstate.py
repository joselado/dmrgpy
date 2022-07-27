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

def gs_energy_single(self,wf0=None,reconverge=None,maxde=None,maxdepth=5):
    """
    Return the ground state energy
    """
    if wf0 is not None: 
        self.execute(wf0.write) # write wavefunction
        self.set_initial_wf(wf0,reconverge=True) # set the initial wavefunction
    if reconverge is not None: # overwrite skip_dmrg_gs
        self.skip_dmrg_gs = not reconverge # if the computation should be rerun
    self.execute(lambda: self.setup_task("GS"))
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.run() # perform the calculation
    # get the ground state energy
    out = self.execute(lambda: np.genfromtxt("GS_ENERGY.OUT"))
    self.e0 = out # store ground state energy
    self.computed_gs = True
    self.sites_from_file = True
    self.gs_from_file = True
    self.skip_dmrg_gs = True
    wf0 = mps.MPS(MBO=self,name="psi_GS.mps").copy() # set GS
    self.set_initial_wf(wf0) # set the initial wavefunction
    if maxde is not None: # enforce a maximum fluctuation in the energy
      e = self.vev(self.hamiltonian)  
      e2 = self.vev(self.hamiltonian,npow=2) 
      de = np.sqrt(abs(e2-e**2)) # fluctuation in the energy
      de = de/self.ns # normalize by the number of sites
      if de>maxde and maxdepth>0: # if a maximum energy fluctuation 
          maxm,nsweeps = self.maxm,self.nsweeps
          noise = self.noise
          print("Energy fluctuation = ",de,maxm)
          self.maxm = maxm*2
          self.nsweeps = 2 # just two sweeps
          self.noise = 0.0
          gs_energy_single(self,maxde=maxde,reconverge=True,
                  maxdepth=maxdepth-1) # execute again
          self.maxm = maxm
          self.nsweeps = nsweeps # restore
          self.noise = noise
    self.computed_gs = True # ground state has been computed
    return out # return energy

def gs_energy(self,policy="single",**kwargs):
    if policy=="single":
        return gs_energy_single(self,**kwargs)
    if policy=="many":
        return gs_energy_many(self,**kwargs)
    else: raise





def lowest_energy(self,h):
    """Return the lowest energy of the Hamiltonian"""
    raise # not finished yet
    self.execute(lambda: h.write("hamiltonian.in")) # write Hamiltonian
    task = {"GS":"true",
            }
    self.task = task
    self.write_task()
    self.run() # perform the calculation
    out = self.execute(lambda: np.genfromtxt("GS_ENERGY.OUT"))
    return out



def get_gs_manifold(MBO,n=2,tol=1e-3,**kwargs):
    """Return the ground state manifold, i.e. all the states with the
    lowest energy"""
    (es,wfs) = MBO.get_excited_states(n=n,**kwargs)
    e0 = es[0] # ground state
    ngs = len(es[np.abs(es-e0)<tol]) # number of ground states
    if ngs<n: # all the GS found
        wfo = []
        for (e,w) in zip(es,wfs):
            if np.abs(e-e0)<tol: wfo.append(w)
        return wfo
    else: 
        print("Recalling with ",n+1,"states")
        return get_gs_manifold(MBO,n=n+1,tol=tol,**kwargs)




