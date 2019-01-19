from __future__ import print_function
import sys
import os

import numpy as np

from dmrgpy import simplechains

b = np.random.random(3)
J = np.random.random(3)

ssc = simplechains.SSC(s=.5,n=14,b=b,J=J)


print("Energy with ED",ssc.gs_energy(mode="ED"))
print("Energy with MPS",ssc.gs_energy(mode="MPS"))
print("Energy with classic DMRG",ssc.gs_energy(mode="classicDMRG"))

ne = 4 # number of excited states

print("Energies with ED",ssc.get_excited(mode="ED",n=ne))
print("Energies with MPS",ssc.get_excited(mode="MPS",n=ne))
print("Energies with classic DMRG",ssc.get_excited(mode="classicDMRG",n=ne))
