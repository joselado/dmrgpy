from . import taskdmrg
from . import mps

def applyoperator(self,A,wf):
    """Apply operator to a many body wavefunction"""
    task = {"applyoperator":"true",
            "applyoperator_wf0":wf.name,
            "applyoperator_multioperator":"applyoperator_multioperator.in",
            "applyoperator_wf1":"applyoperator_wf1.mps",
            }
    self.execute(lambda: A.write(name="applyoperator_multioperator.in"))
    self.task = task
    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
    self.execute( lambda : self.run()) # run calculation
    return mps.MPS(self,name="applyoperator_wf1.mps").copy() # copy

random_mps = mps.random_mps
