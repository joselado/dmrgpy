from manybodychain import Many_Body_Hamiltonian


class Kondo_Chain_Hamiltonian(Many_Body_Hamiltonian):
    """Class for fermionic Hamiltonians"""
    def __init__(self,n):
        sites = []
        for i in range(n): # loop over sites
            sites.append(1) # append a fermion
        for i in range(n): # loop over sites
            sites.append(2) # append a spin
        sites = [1 for i in range(nf)] + spins # fermions and spins
        Many_Body_Hamiltonian.__init__(self,sites)



