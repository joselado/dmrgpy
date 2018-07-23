from __future__ import print_function
import numpy as np
#from numba import jit
from scipy.sparse import csc_matrix
from .read import write_matrix

def generate_basis(spins):
   """Generate the basis for a spin chain"""
   basis = [] # list with the basis
   ns = len(spins) # number of spins
   v = np.array([0 for s in spins]) # initialize
   ms = [int(2*s+1) for s in spins] # number of projections
   while True: # infinite loop
     for i in range(ns-1): # do the special modulus
       if v[i]==(ms[i]): # if reached the maximum
         v[i+1] += 1 # increase the next one
         v[i] = 0 # put in zero
     if v[ns-1]==ms[ns-1]: # if maximum has been reached
       return basis # return the basis
     basis.append(v.copy()) # append this vector
     v[0] += 1 # increase the first value 



def get_dictionary(basis):
   """Get the dictionary relating each vector with the site"""
   d = dict() # create dictionary
   i = 0 # initialize
   for b in basis:
     d[tuple(b)] = i # add entry to the dictionary
     i += 1 # increase
   return d


  


def get_szi(spins,basis,ispin):
  """Get the Sz operator"""
  ii = np.zeros(len(basis),dtype=np.int) # i index 
  jj = np.zeros(len(basis),dtype=np.int) # j index
  vals = np.zeros(len(basis),dtype=np.complex) # value
#  l = (spins[ispin] - 1.)/2. # total spin in this site
  l = spins[ispin]  # total spin in this site
  for ib in range(len(basis)): # loop over basis vectors
    ii[ib] = ib # store index
    jj[ib] = ib # store index
    vals[ib] = basis[ib][ispin] - l # value of the entry
  sz = csc_matrix((vals,(ii,jj)),shape=(len(basis),len(basis))) # get Sz
  return sz



def get_spi(spins,basis,ispin,bdict=None):
  """Get the Sz operators operator"""
  ii = np.zeros(len(basis),dtype=np.int) # i index 
  jj = np.zeros(len(basis),dtype=np.int) # j index
  vals = np.zeros(len(basis),dtype=np.complex) # value
  l = spins[ispin]  # total spin in this site
  if bdict is None: bdict = get_dictionary(basis) # generate the dictionary
  for ib in range(len(basis)): # loop over basis vectors
    ii[ib] = ib # store index
    m = basis[ib][ispin] - l# value of z component
    if basis[ib][ispin] == 2*spins[ispin]: # if reached maximum value
      jj[ib] = ib # put as zero a diagonal element, stupid but harmless
      vals[ib] = 0. # put to zero
      continue # next iteration, this ladder does not go anywhere
    else:  # zero value if you cannot go up
      value = np.sqrt(l*(l+1) - m*(m+1)) # prefactor
      vout = basis[ib].copy() # copy the initial vector
      vout[ispin] += 1 # increase the counter
      jj[ib] =  bdict[tuple(vout)] # get the output vector
      vals[ib] = value # value of the entry
  sp = csc_matrix((vals,(jj,ii)),shape=(len(basis),len(basis))) # get Sz
  return sp


def get_smi(spins,basis,ispin,bdict=None):
  """Get the Sz operators operator"""
  ii = np.zeros(len(basis),dtype=np.int) # i index 
  jj = np.zeros(len(basis),dtype=np.int) # j index
  vals = np.zeros(len(basis),dtype=np.complex) # value
  l = spins[ispin]  # total spin in this site
  if bdict is None: bdict = get_dictionary(basis) # generate the dictionary
  for ib in range(len(basis)): # loop over basis vectors
    ii[ib] = ib # store index
    m = (basis[ib][ispin] - l)# value of z component
    if basis[ib][ispin] == 0: # if reached minimum value
      jj[ib] = ib # put as zero a diagonal element, stupid but harmless
      vals[ib] = 0. # put to zero
      continue # next iteration, this ladder does not go anywhere
    else:  # zero value if you cannot go up
      value = np.sqrt(l*(l+1) - m*(m-1)) # prefactor
      vout = basis[ib].copy() # copy the initial vector
      vout[ispin] -= 1 # decrease the counter
      jj[ib] =  bdict[tuple(vout)] # get the output vector
      vals[ib] = value # value of the entry
  sp = csc_matrix((vals,(jj,ii)),shape=(len(basis),len(basis))) # get Sz
  return sp


def write_basis(basis):
  """Write the basis in basis.out"""
  fo = open("basis.out","w")
  fo.write("# SIZE = "+str(len(basis))+"\n")
  for b in basis: # loop over states
    for ib in b: fo.write(str(ib)+"   ")
    fo.write("\n")
  fo.close()



def generate_chain(spins):
  """Generate the terms for a certain Hamiltonian from a certian basis"""
  basis = generate_basis(spins) # get the basis vectors
  write_basis(basis) # write the basis
  bdict = get_dictionary(basis) # generate dictionary
  sxs,sys,szs = [],[],[]
  for ispin in range(len(spins)): # loop over different spins
    szi = get_szi(spins,basis,ispin) # get sz operator
    smi = get_smi(spins,basis,ispin,bdict=bdict) # get s- operator
    spi = get_spi(spins,basis,ispin,bdict=bdict) # get s+ operator
    sxi = (spi + smi)/2. # sxi
    syi = (spi - smi)/2j # syi
    sxs.append(sxi)
    sys.append(syi)
    szs.append(szi)
  return sxs,sys,szs

def write_chain(spins):
  """Write all the files for the chains"""
  (sxs,sys,szs) = generate_chain(spins) # return the chain
  ispin = 0
  for sxi,syi,szi in zip(sxs,sys,szs):
    write_matrix("sx_"+str(ispin)+"_.op",sxi)
    write_matrix("sy_"+str(ispin)+"_.op",syi)
    write_matrix("sz_"+str(ispin)+"_.op",szi)
    ispin += 1





if __name__=="__main__":
#  basis = generate_basis([4,3]) # return the basis
#  d = get_dictionary(basis)
  chain = generate_chain([4,3]) # return the basis
  print(chain)
