from . import taskdmrg

def applyoperator(mb,A,wf0="psi_GS.mps",wf1="wf1.mps"):
    """Apply operator to a many body wavefunction"""
    task = {"applyoperator":"true",
            "applyoperator_wf0.mps":wf0,
            "applyoperator_wf1.mps":wf1,
            }
    self.task = task
    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
    self.execute( lambda : self.run()) # run calculation



