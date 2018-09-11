from manybodychain import Many_Body_Hamiltonian
import numpy as np
import scipy.linalg as lg


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
    def hamiltonian_free(self,pairs=[[]]):
        """Compute the free correlator"""
        if len(self.hubbard)!=0: raise
        else:
          m = np.zeros((self.ns,self.ns)) # matrix
          for key in self.hoppings:
              t = self.hoppings[key]
              m[t.i,t.j] = t.g
        return m
    def correlator_free(self,pairs=[[]]):
          m = self.hamiltonian_free()
          (es,vs) = lg.eigh(m) # diagonalize
          vs = vs.transpose()
          out = []
          for p in pairs:
              o = 0.0 # initialize
              for (e,v) in zip(es,vs):
                  if e<=0.0: o += v[p[0]]*np.conjugate(v[p[1]]) # add
              out.append(o)
          return np.array(out)*2.0 # return




