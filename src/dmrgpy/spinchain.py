from .manybodychain import Many_Body_Hamiltonian
import numpy as np
from .dmrgpy2pychain import correlator as correlatorpychain

Spin_Hamiltonian = Many_Body_Hamiltonian

class Spin_Hamiltonian(Many_Body_Hamiltonian):
    """Class for spin Hamiltonians"""
    def __init__(self,sites):
        Many_Body_Hamiltonian.__init__(self,sites)
        # default exchange constants
        self.set_exchange(lambda i,j: abs(i-j)==1*1.0)
    def sisj_edge(self):
        if not self.computed_gs: self.get_gs() # compute ground state
        pairs = [(0,i) for i in range(self.ns)] # create pairs
        cs = self.correlator(pairs=pairs,mode="DMRG") # compute
        ns = range(len(cs)) # number of correlators
        return np.array([ns,cs])
    def get_correlator(self,pairs=[[]],mode="DMRG",**kwargs):
        """Return the correlator"""
        if mode=="DMRG": # using DMRG
            return Many_Body_Hamiltonian.get_correlator(self,pairs=pairs,
                    **kwargs)
        elif mode=="ED": # using exact diagonalization
            self.to_folder()
            m = correlatorpychain.correlator(self,pairs=pairs,**kwargs)
            self.to_origin() # go to main folder
            return m
        else: raise
