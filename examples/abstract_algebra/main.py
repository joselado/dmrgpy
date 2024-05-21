# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
SC = spinchain.Spin_Chain(spins) # create the spin chain object

# Let us now create some random states

def random_state():
    """Function returning a random state of the many-body system"""
    # choose one of the two modes
    mode = "MPS" # create an state as a matrix-product state
#    mode = "ED" # create an exact many-body state
    return SC.random_state(mode=mode) # return a random vector


# now generate different wavefunction and operators

# generate two random wavefunctions |1> and |2>

Psi1 = random_state() # generate a random state
Psi2 = random_state() # generate a random state

# now multiply a wavefunction by a complex number z as |3> = z|2>
z = 0.3 + 0.1j # complex number
Psi3 = z*Psi2 # multiply by an scalar

# now define an operator H to act on the wavefunctions

H = 0 # initialize operator
for i in range(n-1): # loop over NN links
    H = H + SC.Sx[i]*SC.Sx[i+1] # XX coupling
    H = H + SC.Sy[i]*SC.Sy[i+1] # YY coupling
    H = H + SC.Sz[i]*SC.Sz[i+1] # ZZ coupling

# now apply the Operator to the vector |2> as |4> = H*|2>

Psi4 = H*Psi2 # apply operator to the state

# if you wish to define a method to do this, you can do it in the following way
Hp = H.copy()
import types
dot = lambda self,v: self*v # apply operator to vector
Hp.dot = types.MethodType(dot, Hp) # add method to the object
Psi4p = Hp.dot(Psi2) # matrix times vector with a method


# states can be summed as if it were conventional vector |5> = |2> + |4>

Psi5 = Psi2 + Psi4 # sum two states

# overlap between vector can be computed with the .dot bound method <2|4>
# return a complex number

O24 = Psi2.dot(Psi4)




