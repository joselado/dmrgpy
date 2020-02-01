import numpy as np


def get_excited(self,n=2,noise=0.0,scale=10.0,fluctuations=False):
    """Return excited state energies"""
    self.get_gs()
    if self.excited_gram_schmidt: sm = "true"
    else: sm = "false"
    task = {"excited":"true",
            "nexcited":str(n),
            "noise":str(noise),
            "excited_gram_schmidt":sm,
            "scale_lagrange_excited":str(scale),
            }
    self.task = task
    self.write_task()
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.run() # perform the calculation
    out = self.execute(lambda: np.genfromtxt("EXCITED.OUT").T)
    if fluctuations: return out
    else: return out[0]

