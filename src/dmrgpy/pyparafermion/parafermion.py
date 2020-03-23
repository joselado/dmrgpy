from ..edtk import edchain
import numpy as np

omega = np.exp(1j*2*np.pi/3)


class Parafermion_Chain(edchain.EDchain):
    """Z3 parafermion chain"""
    def __init__(self,MBO):
        n = MBO.ns # number of sites
        super().__init__() # super initializetion
        # create all operators
        Z = MBO.Z # type of parafermion
        self.Z = Z
        self.localdim = [Z for i in range(n)]
        zero = np.zeros((Z,Z),dtype=np.complex) # get zero
        N = zero.copy()
        for i in range(Z): N[i,i] = i
        for i in range(n): self.create_operator(N,i=i,name="N")
        S = zero.copy()
        for i in range(Z-1): S[i,i+1] = 1.0
        S[Z-1,0] = 1.0
        for i in range(n): self.create_operator(S,i=i,name="Sig")
        for i in range(n): self.create_operator(S.T,i=i,name="SigDag")
        T = zero.copy()
        for i in range(Z): T[i,i] = np.exp(i*1j*2*np.pi/Z)
        for i in range(n): self.create_operator(T,i=i,name="Tau")
        Td = np.conjugate(T.T)
        for i in range(n): self.create_operator(Td,i=i,name="TauDag")
        h = MBO.hamiltonian
        self.hamiltonian = h
        # assign the multioperators
        self.Chi = MBO.Chi
        self.Chid = MBO.Chid
        self.Tau = MBO.Tau
        self.Id = MBO.Id
        self.Taud = MBO.Taud
        self.Psi = MBO.Psi
        self.Psid = MBO.Psid
    def test(self):
        test_commutation(self)


def test_commutation(self):
    """Perform a test of the commutation relations"""
    Chi = [self.get_operator(o) for o in self.Chi]
    Chid = [self.get_operator(o) for o in self.Chid]
    Tau = [self.get_operator(o) for o in self.Tau]
    Taud = [self.get_operator(o) for o in self.Taud]
    Psi = [self.get_operator(o) for o in self.Psi]
    Psid = [self.get_operator(o) for o in self.Psid]
#    Id = self.get_operator(self.Id)
    n = len(Chi) # number of sites
    ntries = 8 # number of tries
    omega = np.exp(1j*2*np.pi/self.Z)
    for ii in range(ntries):
        i = np.random.randint(n)
        j = np.random.randint(n)
        if i>=j: continue
        d = Chi[i]@Chi[j] - omega*Chi[j]@Chi[i]
        if np.max(np.abs(d))>1e-7:
            print("Chi,Chi failed",i,j)
            raise
        d = Psi[i]@Psi[j] - omega*Psi[j]@Psi[i]
        if np.max(np.abs(d))>1e-7:
            print("Psi,Psi failed",i,j)
            raise
        d = Chi[i]@Psi[j] - omega*Psi[j]@Chi[i]
        if np.max(np.abs(d))>1e-7:
            print("Psi,Chi failed",i,j)
            raise
    print("Commutation test passed")





