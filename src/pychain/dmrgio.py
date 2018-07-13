from __future__ import print_function,division
import numpy as np
import os

def savedict(indict,name="chaindict"):
  """Save dictionary in a folder"""
  os.system("rm -rf "+name) # remove folder
  for key in indict: # loop over keys
    obj = indict[key] # get the object
    if type(obj) is np.matrix:
      np.save(name++"/"+key+"_TYPE_.npz"+obj) # save stuff



