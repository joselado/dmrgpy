import numpy as np
import os
from copy import deepcopy
import subprocess


class Wavefunction():
    def __init__(self,mode="DMRG",path=None,name=None,v=None):
        self.mode = mode
        if mode=="DMRG":
            if path is None: self.path = os.getcwd()+"/.mpsfolder/"
            else: self.path = path
            if name is None: self.name = str(np.random())+".mps"
            else: self.name = name
        if mode=="ED": self.v = v # store vector

    def copy(self): # copy wavefunction
        if self.mode=="ED": return np.array(self.v) # return vector
        elif self.mode=="DMRG":
            out = deepcopy(self) # copy object
            name = str(np.random())+".mps" # new name
            subprocess.run(["cp",self.path+"/"+self.name,self.path+"/"+name])
            out.name = name # store the name
            return out # return new object





