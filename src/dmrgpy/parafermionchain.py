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
    def get_ED_obj(self):
        from .pyparafermion import parafermion
        obj = parafermion.Parafermion_Chain(self)
        return obj





