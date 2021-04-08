from __future__ import print_function
from copy import deepcopy
import os
import numpy as np
from . import multioperator

class MPS():
    """Object for an MPS"""
    def __init__(self,MBO=None,name="psi_GS.mps"):
#        self.sc = sc # many body object
#        self.sc.wf0 = None # no wavefunction
        if MBO is None:
            self.path = os.getcwd() # current directory
            self.MBO = None
        else:
            self.path = MBO.path # path to the many body object folder
            self.MBO = MBO
        self.name = name # initial name
#        self.factor = 1.0 # factor of the mps
        self.mps = open(self.path+"/"+name,"rb").read() # read the MPS
        self.sites = open(self.path+"/sites.sites","rb").read() # read sites
    def set_MBO(self,MBO):
        """Set the MBO"""
        self.path = MBO.path # path to the many body object folder
        self.MBO = MBO # set the object
    def dot(self,x):
        if self.MBO is not None: return self.MBO.overlap(self,x)
        else: raise
    def overlap(self,x):
        if self.MBO is not None: return self.MBO.overlap(self,x)
        else: raise
    def aMb(self,M,b):
        if self.MBO is not None: return self.MBO.aMb(self,M,b)
        return self.dot(M*b) # workaround
    def __radd__(self,x): return self + x
    def __add__(self,x):
        if x==0: return self # do nothing
        if self.MBO is not None: return self.MBO.summps(self,x)
        else: raise
    def __sub__(self,x):
        return self + (-1)*x
    def __neg__(self,x):
        return (-1)*x
    def __truediv__(self,x): return self*(1./x)
    def copy(self,name=None):
        """Copy this wavefunction"""
        out = deepcopy(self) # copy everything
        if name is None:
          name = id_generator()+".mps" # create a new name
        out.name = name
#        self.execute(lambda: os.system("cp "+self.name+"  "+out.name))
        return out
    def write(self,name=None,path=None):
        """Write the MPS in a folder"""
        if name is None: name = self.name
        if path is None: path = self.path
        open(path+"/"+name,"wb").write(self.mps) # write the MPS
        open(path+"/sites.sites","wb").write(self.sites) # write the sites
    def get_entropy(self,b=None):
        """Compute entanglement entropy in a bond"""
        if b is None: # compute all 
            return np.mean([self.get_entropy(i) for i in range(1,self.MBO.ns)])
        if self.MBO is not None: return self.MBO.get_entropy(self,b=b)
        else: raise
    def get_site_entropy(self,i):
        if self.MBO is not None: return self.MBO.get_site_entropy(self,i)
        else: raise # not implemented
    def get_bond_entropy(self,i,j):
        if self.MBO is not None: return self.MBO.get_bond_entropy(self,i,j)
        else: raise # not implemented
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
    def normalize(self):
        """Normalize a wavefunction"""
        norm = np.sqrt(self.dot(self).real) # norm
        if norm>1e-8: return self*(1./norm)
        else: return None
    def __rmul__(self,A):
        """Multiply by an operator"""
        if self.MBO is not None:
            if type(A)==multioperator.MultiOperator: # MO type
                return self.MBO.applyoperator(A,self) # apply the operator
            elif multioperator.isnumber(A):
                return A*multioperator.identity()*self
            else: raise
        else: raise
    def __mul__(self,x):
        if multioperator.isnumber(x):
            return x*multioperator.identity()*self
        else: raise





import string
import random

def id_generator(size=20, chars=string.ascii_uppercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))



from .randommps import random_mps
from .randommps import orthogonal_random_mps
from .randommps import random_product_state









