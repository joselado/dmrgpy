from . import spinchain
import numpy as np
from dmrgpy.pychain import rlc
from dmrgpy.pychain import dmrg as classicdmrg

# this is a wrapper to compute simple spin chains

class SSC():
    def __init__(self,s=0.5,n=10,b=[0.,0.,0.],J=[1.,1.,1.]):
        spins = [int(2*s)+1 for i in range(n)]
        self.sc = spinchain.Spin_Hamiltonian(spins) # create the object
        def fj(i,j):
            if abs(i-j)==1: return np.diag(J)
            else: return 0.
        self.sc.set_exchange(fj) # set exchange coupling
        if callable(b): bf = b
        else: bf = lambda i: b
        self.sc.set_fields(bf)
        self.maxm = 20
        # parameters for classic DMRG
        hdict = rlc.monochain(s,b=b,J=J) 
        params = classicdmrg.dmrgdict()
        for key in hdict: params[key] = hdict[key] # copy dictionary
        self.dmrgp = params # store
        self.dmrgp["finite"] = True # finite system
        self.dmrgp["finite_num_ite"] = 3 # iterations finite system
        self.dmrgp["target_length"] = n # system size
    def gs_energy(self,mode="DMRG"):
        """Compute ground state energy using different methods"""
        if mode=="MPS" or mode=="DMRG":
            self.sc.maxm = self.maxm # set bond dimension
            return self.sc.gs_energy(mode="DMRG") # return energy
        elif mode=="ED": return self.sc.gs_energy(mode="ED")
        elif mode=="classicDMRG":
            return classicdmrg.infinite_dmrg(self.dmrgp).energy # return energy
    def get_excited(self,n=2,mode="DMRG"):
        """Compute ground state energy using DMRG"""
        if mode=="MPS" or mode=="DMRG":
            self.sc.maxm = self.maxm # set bond dimension
            return self.sc.get_excited(mode="DMRG",n=n) # return energy
        elif mode=="ED": return self.sc.get_excited(mode="ED",n=n)
        elif mode=="classicDMRG":
            self.dmrgp["retain_states"] = n # states to retain in the DM
            es = classicdmrg.infinite_dmrg(self.dmrgp).energies # return energy
            return np.round(es[0:n],5) # return excitation energies
    def get_gap(self,mode="DMRG"):
        """
        Compute the gap
        """
        es = self.get_excited(mode=mode)
        return es[1]-es[0]
    def get_dm_dis(self,n=4,mode="classicDMRG",f=None):
        """
        Return the distances between the GS DM and excited ones
        """
        if mode=="classicDMRG": # classic DMRG mode
            self.dmrgp["retain_states"] = n # states to retain in the DM
            self.dmrgp["DM_target"] = 0 # states to use in the DM
            self.dmrgp["function_DM_distance"] = f # function to compute DM
                                                   # distance
            out = classicdmrg.infinite_dmrg(self.dmrgp) # perform calculation
            return out.iteration_dmdis # return list
        else: raise # not implemented








