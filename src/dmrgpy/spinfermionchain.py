from . import fermionchain



class Spin_Fermion_Hamiltonian(fermionchain.Spinful_Fermionic_Hamiltonian):
    def __init__(self,sites):
        """Create the sites"""
        super().__init__(len(sites)) # initialize the Hamiltonian
