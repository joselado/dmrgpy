from __future__ import print_function
import numpy as np
import scipy.sparse.linalg as slg
import scipy.sparse as sp
import tensorial

# functions to perform posprocessing on the DMRG results
def correlator(indict,outfile="CORRELATORS.OUT"):
  """Calculates the correlation function between different
  operators in different sites"""
  bops = indict["block_operators"] # list with operators
  ns = len(bops) # number of sites
  dmat = indict["density_matrix"] # density matrix
  # now calculate the different correlators
  corrs = [] # empty list
  sops = indict["site_operators"]
  for i in range(len(bops)):
    bopsi = bops[i] # get this operators
    Os = []
    for (o1,o2) in zip(sops,bopsi): Os.append(o1*o2) # store
    corrs.append([float(ns-i)]+[np.sum((dmat*O).diagonal()).real for O in Os]) # compute
  corrs = np.matrix(corrs)
  # now write in file
  np.savetxt(outfile,corrs)
  print("Correlators written in ",outfile)



def dynamical_correlator(indict,keyfile="DYNAMICAL_CORRELATOR_",
                          energies=np.linspace(0.,3.,50)):
  """Calculates the correlation function between different
  operators in different sites, with energy dependency"""
  if not indict["has_full_operators"]: return
  lBops = indict["left_block_operators_full"] # list with operators
  ns = len(lBops) # number of sites in the block
  lSops = indict["left_site_operators_full"] # left site operators
  wf = indict["wf"] # get the wavefunction
  e0 = indict["energy"] # get the wavefunction
  h = indict["hamiltonian"] # get the wavefunction
  for i in range(ns): # loop over sites in the block
    lBs = lBops[i] # get the block operators
    name = keyfile + str(ns-i)+"_.OUT" # name of the key file
    fo = open(name,"w") # open the file
    for e in energies: # loop over frequencies
      fo.write(str(e)+"  ")
      for (o1,o2) in zip(lSops,lBs): # loop over pairs of operators
        cf = dyncorr(o1,o2,e0,wf,h,omega=e,eta=0.01)
        fo.write(str(cf)+"   ") # write in file
      fo.write("\n")
    fo.close() # close file
    print("DYNCORR written in ",name)











def dyncorr(a,b,e0,v0,h,omega=0.0,eta=0.01):
  """Calculates a dynamical correlator. This funcion impements
  <0|A(E+omega+ieta - H)^(-1)B|0> """
  bv = b*v0 # right term
  iden = tensorial.identity(len(v0)) # identity matrix
  g = iden*(omega+e0+1j*eta) - h # matrix to invert
  g = tensorial.slo2lo(g)
#  print(type(g)) ; exit()
  (gbv,info) = slg.lgmres(g,bv) # solve the equation  
  cf = np.conjugate(v0).dot(a*gbv) # correlation function
  return cf.imag # return the expectation value
