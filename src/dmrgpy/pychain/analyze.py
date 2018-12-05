from __future__ import print_function
import os
import numpy as np
import spectrum


fform = "{0:.5f}".format

def spins(h,sc,k=1):
  """Calculate the expectation value of spin operators in the ground state""" 
  (es,vs) = spectrum.eigenstates(h,evals=True,k=k)
  es = es-min(es)
  fo = open("SPIN_VALUE.OUT","w")
  for (e,v) in zip(es,vs):
    fo.write(fform(e)+"   ") # write energy
    for si in range(len(sc.sxi)):
      sx = spectrum.exp_val(v,sc.sxi[si])
      sy = spectrum.exp_val(v,sc.syi[si])
      sz = spectrum.exp_val(v,sc.szi[si])
      fo.write(fform(sx)+"   ")
      fo.write(fform(sy)+"   ")
      fo.write(fform(sz)+"   ")
    fo.write("\n")
  fo.close()



