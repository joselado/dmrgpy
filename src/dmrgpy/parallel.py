# routines to call a function in parallel
from __future__ import print_function
import scipy.linalg as lg
from . import algebra
import os
try:
  from multiprocess import Pool
  import multiprocess
  maxcpu = multiprocess.cpu_count()
except:
    print("Multiprocess not working")
    def Pool(n=1): # workaround
            class mpool():
                def map(self,f,xs):
                  return [f(x) for x in xs]
                def terminate(self): return None # dummy function
            return mpool()
    maxcpu = 1 # maximum number of cpus

cores = 1 # call in a single by default


def set_cores(n=1):
    global cores
    cores = n



def multicall(fs):
    """
    Execute several independent functions
    """
    def ff(i):
        return fs[i]() # execute that function
    return pcall(ff,range(len(fs)))



def pcall_serial(fun,args):
  """Function to call in serial"""
  return [fun(a) for a in args]



def pcall_mp(fun,args,cores=cores):
    """Calls a function for every input in args"""
    mainpool = Pool(cores) # create pool
    out = mainpool.map(fun,args) # return list
    mainpool.terminate()
    del mainpool # delete pool
    return out



def pcall(fun,args): # define the function
  path = os.getcwd() # get the current directory
  global cores
  if cores==1: 
      out = pcall_serial(fun,args) # one core, simply iterate
  else: 
      out = pcall_mp(fun,args,cores=cores) # call in parallel
  os.chdir(path) # go to the original path
  return out # return outputs

