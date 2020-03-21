from ..edtk import edchain
import numpy as np

class Spin_Chain(edchain.EDchain):
    """Z3 parafermion chain"""
    def __init__(self,MBO):
        n = MBO.ns # number of sites
        super().__init__() # super initializetion
        self.localdim = [2 for i in range(n)]
        self.hamiltonian = MBO.hamiltonian


