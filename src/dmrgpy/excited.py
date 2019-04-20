import numpy as np


def get_excited(self,n=2,noise=0.0,scale=1.0):
    """Return excited state energies"""
    task = {"excited":"true",
            "nexcited":str(n),
            "noise":str(noise),
            "scale_lagrange_excited":str(scale),
            }
    self.task = task
    self.write_task()
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.run() # perform the calculation
    out = self.execute(lambda: np.genfromtxt("EXCITED.OUT"))
    return out

