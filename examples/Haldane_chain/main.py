# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# you may need to uncomment next line to add the library
from dmrgpy import simplechains # basic library for spin chains
s = 1 # spin in the chain
n = 10 # length of the chain
ssc = simplechains.SSC(s=s,n=n) # get the spin chain object
# solve with exact diagonalization
print("Energy with ED",ssc.gs_energy(mode="ED"))
# solve with ITensor MPS algorith
print("Energy with MPS",ssc.gs_energy(mode="MPS"))
# solve with classical DMRG algorithm (slow)
print("Energy with classic DMRG",ssc.gs_energy(mode="classicDMRG"))
# compute excitations of the model
ne = 4 # number of states to compute
print("Energies with ED",ssc.get_excited(mode="ED",n=ne))
print("Energies with MPS",ssc.get_excited(mode="MPS",n=ne))
print("Energies with classic DMRG",ssc.get_excited(mode="classicDMRG",n=ne))











