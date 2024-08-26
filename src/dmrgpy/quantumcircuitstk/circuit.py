from ..spinchain import Spin_Chain


class Quantum_Circuit():
    def __init__(self,n):
        ED = True
        sc = Spin_Chain(["S=1/2" for i in range(n)])
        if ED=True:
            sc.MBO = sc.get_ED_obj() # get the object
        self.X = 2*sc.Sx # Pauli operators
        self.Y = 2*sc.Sy
        self.Z = 2*sc.Sz





def generate_ansatz(self):
    """Return a callable function"""
    raise




