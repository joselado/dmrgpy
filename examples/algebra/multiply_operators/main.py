# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

def get_traces(n):
    """Build a Heisenberg chain of length n plus non-Hermitian terms,
    convert to an MPO, and return the traces of a few operator products."""
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    Sx,Sy,Sz = sc.Sx,sc.Sy,sc.Sz # lists with the respective operators

    # generate a Heisenberg operator as example
    A0 = 0
    for i in range(n-1):
        A0 = A0 + Sx[i]*Sx[i+1]
        A0 = A0 + Sy[i]*Sy[i+1]
        A0 = A0 + Sz[i]*Sz[i+1]

    for i in range(n): # add some terms that make the operator non-Hermitian
        A0 = A0 + 0.3*1j*Sy[i]*Sx[i]
    A0 = A0 + 0.2*1j # add a diagonal imaginary
    # this controls the accuracy of the tensor-network calculation
    # sc.mpomaxm = 80
    # sc.maxm = 80
    # (the higher, the more accurate) typical minimal values are 50-100

    # in order to make products in a scalable way, convert to a new object
    A = sc.toMPO(A0) # transforms into an MPO (that can be multiplied)
    # the object A is a "StaticOperator" you can do A*B  but not A + B
    #A = sc.toMPO(A0,mode="ED") # uncomment to test it with ED using this instead

    Ad = A.get_dagger() # compute the Hermitian

    # define a few operators via multiplying
    AA = A*A
    AdA = Ad*A
    AAd = A*Ad
    AAdA = A*Ad*A

    return A.trace(),Ad.trace(),AA.trace(),AAd.trace(),AAdA.trace()

ns = [2,4,6,8] # sweep the chain length (kept small, traces grow with Hilbert space)
traces = {"A":[],"A^dagger":[],"A*A":[],"A*A^dagger":[],"A*A^dagger*A":[]}
labels = list(traces.keys())
for n in ns:
    trA,trAd,trAA,trAAd,trAAdA = get_traces(n)
    print("n =",n)
    print("Trace of A",trA)
    print("Trace of A^\dagger",trAd)
    print("Trace of A*A",trAA)
    print("Trace of A*A^\dagger",trAAd)
    print("Trace of A*A^\dagger*A",trAAdA)
    for key,val in zip(labels,[trA,trAd,trAA,trAAd,trAAdA]):
        traces[key].append(abs(val))

import matplotlib.pyplot as plt
for key in labels:
    plt.plot(ns,traces[key],marker="o",label=key)
plt.xlabel("Chain length n")
plt.ylabel("|Trace|")
plt.yscale("log")
plt.legend()
plt.show()
