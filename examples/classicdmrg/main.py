# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy.pychain import dmrg
from dmrgpy.pychain import rlc
import scipy.linalg as lg
ns = 40 # number of spins
maxm = 20 # bond dimension
b = np.random.random(3)*0.
b = [0.,0.,0.1]
hdict = rlc.monochain(1.,d=0.,b=b,fun = lambda i: 1.)
def f1(d,idmrg=0):
  if idmrg>40: return 4
  else: return 1
def f2(d,idmrg=0):
  if idmrg>40: return 40
  else: return 20
# input
params = dmrg.dmrgdict()
for key in hdict: params[key] = hdict[key] # copy dictionary
params["number_of_states"] = maxm
params["diag_states"] = 1
params["diag_mode"] = "arpack"
params["target_function"] = None
#params["dynamic_number_of_states"] = f2
params["retain_states"] = 10
#params["dynamic_retain_states"] = f1
params["target_state"] = 0
params["tol"] = 0.00001
#params["ensure_symmetry"] = True
params["ndmrg"] = 1000
params["target_length"] = ns
params["store_operators"] = False
params["finite_num_ite"] = 1
params["avoid_edge"] = 0
params["finite"] = True
#params["target_function"] = rlc.target_total_s(sz=0)
# perform calculation
import time
t0 = time.clock()
import profile
#profile.run('dmrg.infinite_dmrg(params)',sort=1)
out = dmrg.infinite_dmrg(params) # perform calculation
print("Energy with classic DMRG",out.energy)
#print(out.iteration_dmdis)
t1 = time.clock()
e0 = out.energy
# Now use ITensor
from dmrgpy import spinchain
sc = spinchain.Spin_Hamiltonian([2 for i in range(ns)])
sc.set_fields(lambda i: b)
sc.maxm = maxm # bond dimension
e1 = sc.gs_energy()
print("Energy with MPS",e1)
t2 = time.clock()
print("Time in classic DMRG",t1-t0)
print("Time in MPS",t2-t1)
if abs(e0-e1)>0.01: raise


