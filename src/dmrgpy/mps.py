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
  def copy(self):
    """Copy this wavefunction"""
    out = deepcopy(self) # copy everything
    name = id_generator()+".mps" # create a new name
    out.name = name
    os.system("cp "+self.path+"/"+self.name+"  "+out.path+out.name)
    return out
  def __mul__(self,x):
    """Multiply by an scalar"""
    out = self.copy()
    out.factor *= x # mutiply
    return out





import string
import random

def id_generator(size=20, chars=string.ascii_uppercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))




def dot(wf1,wf2):
  """Compute scalar product"""
  raise
  sc = wf1.sc
  sc.setup_task(mode="overlap") # compute overlap mode
  os.system("cp "+wf1.path+wf1.name+"  "+wf1.path+"/overlap_wf1.mps") # copy
  os.system("cp "+wf2.path+wf2.name+"  "+wf2.path+"/overlap_wf2.mps") # copy
  os.system("rm -f "+sc.path+"/OVERLAP.OUT") # remove the file
  sc.run(automatic=True) # run the calculation
  out = np.genfromtxt(sc.path+"/OVERLAP.OUT") # get the data
  return (out[0]+1j*out[1])*wf1.factor*wf2.factor # return result




