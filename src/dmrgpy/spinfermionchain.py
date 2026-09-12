from . import fermionchain



class Spin_Fermion_Hamiltonian(fermionchain.Spinful_Fermionic_Hamiltonian):
    def __init__(self,sites,**kwargs):
        """Create the sites"""
        # **kwargs (itensor_version=, ...) is forwarded the same way the
        # fermionic subclasses in fermionchain.py forward it
        super().__init__(len(sites),**kwargs) # initialize the Hamiltonian
