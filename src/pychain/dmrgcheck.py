import build
import spectrum
import rlc
import numpy as np
import dmrg

# library to compare CI with DMRG

def heisenberg_ci(spins=[.5,.5]):
  """Model for a Heisenberg chain"""
  sc = build.Spin_chain() # create class
  sc.build(spins) # create the object
  h = sc.template(name="open_chain")
  (e,w) = spectrum.ground_state(h) 
  return e


def heisenberg_dmrg(spins=[.5,.5]):
  """Model for a Heisenberg chain"""
  hdict = rlc.monochain(spins[0],d=0.,b=[0.,0.,0.],fun = lambda i: 1.)
  datadict = dmrg.dmrgdict()
  for key in hdict: datadict[key] = hdict[key] # copy dictionary
  datadict["target_length"] = len(spins)
  datadict["finite_num_ite"] = 3
  datadict["finite"] = True
  dmrg.infinite_dmrg(datadict)
  es = np.genfromtxt("ENERGY.OUT")[-1]
  try: return es[1]
  except: return es

