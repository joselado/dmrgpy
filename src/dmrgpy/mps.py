from __future__ import print_function
from copy import deepcopy
import os
import numpy as np

class MPS():
    """Object for an MPS"""
    def __init__(self,sc,name="psi_GS.mps"):
        self.sc = sc # many body object
        self.sc.wf0 = None # no wavefunction
        self.path = sc.path # path to the spin chain folder
        self.name = name # initial name
        self.factor = 1.0 # factor of the mps
    def dot(self,x):
        return dot(self,x) # dot function
    def copy(self,name=None):
        """Copy this wavefunction"""
        out = deepcopy(self) # copy everything
        if name is None:
          name = id_generator()+".mps" # create a new name
        out.name = name
        self.execute(lambda: os.system("cp "+self.name+"  "+out.name))
        return out
    write = copy
    def rename(self,name):
        self.execute(lambda: os.system("mv "+self.name+"  "+name))
        self.name = name
    def execute(self,f):
        pwd = os.getcwd() # path
        os.chdir(self.path) # go
        f()
        os.chdir(pwd) # go back
    def clean(self):
        self.execute(lambda: os.system("rm "+self.name))
        del self
    def __mul__(self,x):
        """Multiply by an scalar"""
        out = self.copy()
        out.factor *= x # multiply
        return out





import string
import random

def id_generator(size=20, chars=string.ascii_uppercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))





