from __future__ import print_function
import sys
import os

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy.pychain import dmrg
from dmrgpy.pychain import rlc
import scipy.linalg as lg

hdict = rlc.frusladder(.5,j=-4,jf=1.)
hdict = rlc.biladder(.5,.5)
hdict = rlc.monochain(0.5,d=0.,b=[0.,0.,0.],fun = lambda i: 1.+0.2*(-1)**i)

def f1(d,idmrg=0):
  if idmrg>40: return 4
  else: return 1

def f2(d,idmrg=0):
  if idmrg>40: return 40
  else: return 20

# input
params = dmrg.dmrgdict()
for key in hdict: params[key] = hdict[key] # copy dictionary
params["number_of_states"] = 40
params["diag_states"] = 1
params["diag_mode"] = "arpack"
params["target_function"] = None

#params["dynamic_number_of_states"] = f2
params["retain_states"] = 4
#params["dynamic_retain_states"] = f1
params["target_state"] = 0
params["tol"] = 0.00001
params["ensure_symmetry"] = True
params["ndmrg"] = 1000
params["target_length"] = 500
params["store_operators"] = False
params["finite_num_ite"] = 3
params["avoid_edge"] = 0
params["finite"] = True
#params["target_function"] = rlc.target_total_s(sz=0)


# perform calculation

dmrg.infinite_dmrg(params) # perform calculation
