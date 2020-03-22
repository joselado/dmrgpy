from ..edtk import edchain
import numpy as np

omega = np.exp(1j*2*np.pi/3)


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
#        d = C[i]@C[j] + C[j]@C[i]
#        if np.max(np.abs(d))>1e-7:
#            print("C,C failed",i,j)
#            raise
#        if i==j: d = Cdag[i]@C[j] + C[j]@Cdag[i] - Id
#        else: d = Cdag[i]@C[j] + C[j]@Cdag[i]
#        if np.max(np.abs(d))>1e-7:
#            print("Cdag,C failed",i,j)
#            raise
    print("Commutation test passed")





