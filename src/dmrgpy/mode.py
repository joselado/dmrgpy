import os
from . import taskdmrg

dmrgpath = os.path.dirname(os.path.realpath(__file__)) # path for the program
mpscpp = dmrgpath+"/mpscpp2/mpscpp.x" # C++ executable

def run(self,automatic=False):
    """
    Run the DMRG calculation
    """
    if get_mode(self)=="ED": 
        return # do nothing
    # executable
    self.execute(lambda : taskdmrg.write_tasks(self)) # write tasks
    if self.itensor_version in [2,"2","v2","C++","cpp","c","C"]:
        if os.path.isfile(mpscpp): 
            from . import cpprun
            cpprun.run(self) # run the C++ version
            return
#    elif self.itensor_version==3: mpscpp = dmrgpath+"/mpscpp3/mpscpp.x" 
    elif self.itensor_version in ["julia","Julia","jl"]:
        from . import juliarun
        juliarun.run(self)
        return
    else: raise
#      if not os.path.isfile(mpscpp): # mpscpp.x not found, rerun with julia
#          print("C++ backend not found, trying to run with Julia version")
#          self.setup_julia() # turn to Julia
#          return self.run() # rerun with julia



def get_mode(self,mode="DMRG"):
    """Return the mode of the calculation"""
    # if there is no C++, then use ED
    if not os.path.isfile(mpscpp): 
        print("C++ not compiled, using default ED routines")
        return "ED" # use exact diagonalization
    # if there is an enforced mode, then use that one
    if self.mode is not None: return self.mode # use the enforced mode
    else: return mode # use default


