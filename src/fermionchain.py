from manybodychain import Many_Body_Hamiltonian


class Fermionic_Hamiltonian(Many_Body_Hamiltonian):
    """Class for fermionic Hamiltonians"""
    def __init__(self,n):
        Many_Body_Hamiltonian.__init__(self,[1 for i in range(n)])



