from . import taskdmrg
from . import mps

def applyoperator(self,A,wf):
    """Apply operator to a many body wavefunction"""
    task = {"applyoperator":"true",
            "applyoperator_wf0":wf.name,
            "applyoperator_wf1":"applyoperator_wf1.mps",
            }
    self.task = task
    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
    self.execute( lambda : self.run()) # run calculation
    return mps.MPS(self,name="applyoperator_wf1.mps").copy() # copy



