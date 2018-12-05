from __future__ import print_function
import numpy as np

import build

def sc_template(sc=None,spins=None,name="open_chain",j=1.0):
  """Generate the Hamiltonian of a certain spin model, using
  the input spin chain"""
  if sc is None:
    if spins is None: raise
    sc = build.Spin_chain() # create class
    sc.build(spins) # create the object
  if name=="open_chain": # open chain
    xs = [i for i in range(sc.nspins-1)]
    ys = [i for i in range(1,sc.nspins)]
    js = [j for i in range(1,sc.nspins)]
    print(js)
    return sc.generate_hamiltonian(xs,ys,js) 
  elif name=="closed_chain": # open chain
    xs = [i for i in range(sc.nspins)]
    ys = [i for i in range(1,sc.nspins)]+[0]
    js = [j for i in range(sc.nspins)]
    return sc.generate_hamiltonian(xs,ys,js) 


def get_chain(spins,jcall):
  """Return a spin chain with the following js"""
  xs,ys,js = [],[],[]
  for ii in range(len(spins)): # loop over spins
    for jj in range(len(spins)): # loop over spins
      j = jcall(ii,jj) # get this coupling
      if j != 0.0: # if non-zero
        xs.append(ii)
        ys.append(jj)
        js.append(j)
  np.savetxt("COUPLINGS.OUT",np.matrix([xs,ys,js]).T)
  sc = build.Spin_chain() # create class
  sc.build(spins) # create the object
  return (sc,sc.generate_hamiltonian(xs,ys,js)) # return the chain

  







