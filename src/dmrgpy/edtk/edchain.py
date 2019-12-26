from ..algebra import algebra
from .. import multioperator


class EDchain():
    """Generic class for an ED chain"""
    def __init__(self):
        self.computed_gs = False
    def get_operator(self,name,i=0):
        """Return an operator"""
        return self.operators[(name,i)] # return the operator
    def get_hamiltonian(self):
        """Return the Hamiltonian"""
        return multioperator.MO2matrix(self.hamiltonian,self) # return operator
    def gs_energy(self):
        """Return ground state energy"""
        h = self.get_hamiltonian()
        return algebra.ground_state(h)[0]
    def get_gs(self):
        """Get ground state wavefunction"""
        if self.computed_gs: return self.wf0
        else: 
          e0,wf0 = algebra.ground_state(self.get_hamiltonian())
          self.wf0 = wf0
          self.computed_gs = True
          return self.wf0
    def vev(self,op):
        """Return a vacuum expectation value"""
        wf0 = self.get_gs()
        op = multioperator.MO2matrix(op,self) # return operator
        return algebra.braket_wAw(wf0,op)
    def get_excited(self,**kwargs):
        """Excited states"""
        h = self.get_hamiltonian()
        return algebra.lowest_eigenvalues(h,**kwargs)

