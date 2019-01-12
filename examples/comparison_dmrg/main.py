from __future__ import print_function
import sys
import os

import numpy as np

from dmrgpy import simplechains

b = np.random.random(3)

ssc = simplechains.SSC(s=1,n=10,b=b)


print("Energy with ED",ssc.gs_energy(mode="ED"))
print("Energy with MPS",ssc.gs_energy(mode="MPS"))
print("Energy with classic DMRG",ssc.gs_energy(mode="classicDMRG"))


