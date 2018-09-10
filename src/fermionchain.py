from manybodychain import Many_Body_Hamiltonian
import numpy as np


class Fermionic_Hamiltonian(Many_Body_Hamiltonian):
    """Class for fermionic Hamiltonians"""
    def __init__(self,n):
        Many_Body_Hamiltonian.__init__(self,[1 for i in range(n)])
    def get_density(self):
        """Return the electronic density"""
        m = self.get_file("MEASURE_N.OUT") # get the file
        return m.transpose()[1] # return density
    def get_density_fluctuation(self):
        """Return the electronic density"""
        d = self.get_file("MEASURE_N.OUT").transpose()[1] # get the file
        d2 = self.get_file("MEASURE_N2.OUT").transpose()[1] # get the file
        return d2-d**2 # return density fluctuations
    def get_delta(self):
        """Return the electronic density"""
        m = self.get_file("MEASURE_DELTA.OUT") # get the file
        return m.transpose()[1] # return delta



