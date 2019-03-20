import numpy as np


def get_excited(self,n=2,noise=0.0):
    """Return excited state energies"""
    self.to_folder() # go to temporal folder
    task = {"nexcited":str(n),
            "noise":str(noise),
            }
    self.setup_task("excited",task=task)
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.run() # perform the calculation
    out = self.execute(lambda: np.genfromtxt("EXCITED.OUT"))
    self.to_origin() # go to main folder
    return out

