from .manybodychain import Many_Body_Chain
import numpy as np


class Parafermionic_Chain(Many_Body_Chain):
    """Class for a parafermionic chain"""
    def __init__(self,n):
        self.N = [self.get_operator("N",i) for i in range(n)]
        self.Sig = [self.get_operator("Sig",i) for i in range(n)]
        self.Sigd = [self.get_operator("SigDag",i) for i in range(n)]
        self.Tau = [self.get_operator("Tau",i) for i in range(n)]
        self.Taud = [self.get_operator("TauDag",i) for i in range(n)]
        self.Id = self.get_operator("Id",1)
        Many_Body_Chain.__init__(self,[-2 for i in range(n)])
        self.use_ampo_hamiltonian = True # use ampo
        self.Chi = []
        self.Psi = []
        for i in range(n): 
            t = 1
            for j in range(i): t = t*self.Tau[j]
            self.Chi.append(t*self.Sig[i])
            self.Psi.append(t*self.Sig[i]*self.Tau[i])
    def get_ED_obj(self):
        from .pyparafermion import parafermion
        obj = parafermion.Parafermion_Chain(self)
        return obj
    def get_dynamical_correlator(self,mode="DMRG",**kwargs):
        if mode=="DMRG":
            return super().get_dynamical_correlator_MB(**kwargs)
        elif mode=="ED":
            return self.get_ED_obj().get_dynamical_correlator(**kwargs)





