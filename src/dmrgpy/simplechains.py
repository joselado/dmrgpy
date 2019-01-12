from . import spinchain
from dmrgpy.pychain import rlc
from dmrgpy.pychain import dmrg as classicdmrg

# this is a wrapper to compute simple spin chains

class SSC():
    def __init__(self,s=0.5,n=10,b=[0.,0.,0.]):
        spins = [int(2*s)+1 for i in range(n)]
        self.sc = spinchain.Spin_Hamiltonian(spins) # create the object
        if callable(b): bf = b
        else: bf = lambda i: b
        self.sc.set_fields(bf)
        self.maxm = 20
        # parameters for classic DMRG
        hdict = rlc.monochain(s,b=b) 
        params = classicdmrg.dmrgdict()
        for key in hdict: params[key] = hdict[key] # copy dictionary
        self.dmrgp = params # store
        self.dmrgp["finite"] = True # finite system
        self.dmrgp["target_length"] = n # system size
    def gs_energy(self,mode="DMRG"):
        """Compute ground state energy using DMRG"""
        if mode=="MPS" or mode=="DMRG":
            self.sc.maxm = self.maxm # set bond dimension
            return self.sc.gs_energy(mode="DMRG") # return energy
        elif mode=="ED": return self.sc.gs_energy(mode="ED")
        elif mode=="classicDMRG":
            return classicdmrg.infinite_dmrg(self.dmrgp).energy # return energy




