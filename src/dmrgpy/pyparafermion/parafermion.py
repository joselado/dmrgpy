from ..edtk import edchain
import numpy as np

class Parafermion_Chain(edchain.EDchain):
    """Z3 parafermion chain"""
    def __init__(self,MBO):
        n = MBO.ns # number of sites
        super().__init__() # super initializetion
        self.localdim = [3 for i in range(n)]
        # create all operators
        zero = np.zeros((3,3),dtype=np.complex) # get zero
        N = zero.copy()
        N[1,1] = 1.0
        N[2,2] = 2.0
        for i in range(n): self.create_operator(N,i=i,name="N")
        S = zero.copy()
        S[0,2] = 1.0
        S[1,0] = 1.0
        S[2,1] = 1.0
        S = S.T # be consistent with the itensor basis
        for i in range(n): self.create_operator(S,i=i,name="Sig")
        for i in range(n): self.create_operator(S.T,i=i,name="SigDag")
        T = zero.copy()
        T[0,0] = 1.0
        T[1,1] = np.exp(1j*2*np.pi/3.)
        T[2,2] = np.exp(1j*4*np.pi/3.)
        for i in range(n): self.create_operator(T,i=i,name="Tau")
        Td = np.conjugate(T.T)
        for i in range(n): self.create_operator(Td,i=i,name="TauDag")
        h = MBO.hamiltonian
        self.hamiltonian = h


