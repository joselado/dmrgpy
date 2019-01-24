# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import simplechains
# This example shows the quantum phase transition in the trnasverse Ising model
# generate a 1D Ising chain
def get_dm(bx):
    """
    Compute the distance between density matrices
    """
    b = [bx,0.,0.] # magnetic field
    J = [0.,0.,1.] # Ising coupling
    ssc = simplechains.SSC(s=.5,n=20,b=b,J=J) # generate the Ising chain
    #### You can define your function that computes
    ### the distance between density matrices, input is a list of DM,
    ### it should return a list###
    def fDM(ds): 
        """This example uses the commutator"""
        def f(d1,d2): 
            A=d1*d2;B=d2*d1;C=A-B; C=C*C.H ; return C.trace()[0,0].real
        return [f(ds[0],ds[i]) for i in range(1,len(ds))] # return distances
    o = ssc.get_dm_dis(n=10,f=fDM) # return distance between DM
    o = np.mean(o) # perform average of the distances
    print("Computing B_x = ",bx,"result",o)
    return o
bs = np.linspace(0.2,2.,5) # list of fields
gs = [get_dm(b) for b in bs] # list of distances
import matplotlib.pyplot as plt
plt.plot(bs,gs,marker="o")
plt.xlabel("Bx")
plt.ylabel("Distance between DM matrix")
plt.show()
