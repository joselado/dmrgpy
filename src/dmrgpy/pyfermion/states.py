from __future__ import print_function
import numpy as np
from scipy.sparse import csc_matrix

dimmax = 50000 # maximum dimension of the matrix


def constrain_nelectrons(ne=1):
  """Return a function that will constrain the total number of
  electrons to ne"""
  def f(v):
    if np.sum(v)==ne: return True
    else: return False
  return f



def get_basis(nsites,nelectrons):
  """Get the basis for this many electrons in that many orbitals"""
  f = constrain_nelectrons(ne=nelectrons) # function to constain
  return generate_basis(nsites,constrain=f)



def generate_basis(ne,constrain=lambda x: True):
  """Generate basis for ne fermions"""
  basis = [] # list with the basis
  v = np.array([0 for s in range(ne)]) # initialize
  ms = [2 for s in range(ne)] # number of projections
  while True: # infinite loop
    for i in range(ne-1): # loop over orbitals
      if v[i]>=ms[i]: # if reached the maximum number of states in this site
        v[i+1] += 1 # increase the next one
        v[i] = 0 # put in zero
    if v[ne-1]==ms[ne-1]: # if maximum has been reached
      return basis # return the basis 
    if constrain(v): # check if this state is accepted
      basis.append(v.copy()) # append this vector
    v[0] += 1 # increase the first value for the next iteration
  

def get_dictionary(basis):
   """Get the dictionary relating each vector with the site"""
   d = dict() # create dictionary
   i = 0 # initialize
   for b in basis:
     d[tuple(b)] = i # add entry to the dictionary
     i += 1 # increase
   return d



def write_basis(basis):
  """Write the basis in basis.out"""
  fo = open("basis.out","w")
  fo.write("# SIZE = "+str(len(basis))+"\n")
  for b in basis: # loop over states
    for ib in b: fo.write(str(ib)+"   ")
    fo.write("\n")
  fo.close()


def fermi_sign(v,i,j):
  """Return the sign coming from statistics"""
  if v[i]==0: return 0. # next iteration, destructor empty
  if (i!=j) and v[j]==1: return 0. # next iteration, creator filled
  fac = 1 # start in 1
  for ii in range(i): # how many before destroying 
    if v[ii]==1: fac *= -1 # multiply by -1
  v[i] = 0 # destroy electron
  for jj in range(j): # how many before creating 
    if v[jj]==1: fac *= -1 # multiply by -1
  return fac


def one2many_basis(m,basis,bdict=None):
  """Convert a one body operator into a manybody operator"""
  if bdict is None: bdict = get_dictionary(basis) # get the dictionary
  nm = m.shape[0] # dimension of one body
  if nm!=len(basis[0]): raise # error if wrong dimensions
  on = len(basis) # dimension of output matrix
  nout = len(basis) # length of the output matrix
  mout = csc_matrix(([],([],[])),shape=(on,on)) # output matrix
  for i in range(nm): # loop over states
    for j in range(nm): # loop over states
      a = m[j,i] # matrix element
      ii = np.zeros(on,dtype=np.int) # indexes
      jj = np.zeros(on,dtype=np.int) # indexes
      vals = np.zeros(on,dtype=np.complex) # values
      for ib in range(len(basis)): # loop over basis elements
        b = basis[ib] # get the vector
        if b[i]==0: continue # next iteration, destructor empty
        if (i!=j) and b[j]==1: continue # next iteration, creator filled
        bo = b.copy() # copy vector
        bo[i] = 0 # empty the level
        bo[j] = 1 # fill the level
        jb = bdict[tuple(bo)] # get index of the out vector
     #   except: continue # state is not in basis
        ii[ib] = ib # store index 
        jj[ib] = jb # store index 
        vals[ib] = a*fermi_sign(b.copy(),i,j) # store value
      mi = csc_matrix((vals,(jj,ii)),shape=(on,on)) #  matrix
      mout = mout + mi # add contribution
  return mout # return matrix
  

def four2many(m,basis,bdict=None):
  """Convert a four fermion operator into a many body operator"""
  if bdict is None: bscit = get_dictionary(basis) # get the dictionary
  nm = m.shape[0] # dimension of one body
  if nm!=len(basis[0]): raise # error if wrong dimensions
  on = len(basis) # dimension of output matrix
  mout = csc_matrix(shape=(on,on)) # output matrix
  for i in range(nm): # loop over states
    for j in range(nm): # loop over states
      for k in range(nm): # loop over states
        for l in range(nm): # loop over states
          a = m[l][k][j][i] # matrix element
          ii = np.zeros(on,dtype=np.int) # indexes
          jj = np.zeros(on,dtype=np.int) # indexes
          vals = np.zeros(on,dtype=np.complex) # values
          for ib in range(len(basis)): # loop over basis elements
            b = basis[ib] # get the vector
            if b[i]==0: continue # next iteration
            if b[j]==0: continue # next iteration
            if ((i!=k) and (j!=k)) and b[k]==1: continue # next iteration
            if ((i!=l) and (j!=l)) and b[l]==1: continue # next iteration
            bo = b # copy vector
            bo[i] = 0 # empty the level
            bo[j] = 0 # empty the level
            bo[k] = 1 # fill the level
            bo[l] = 1 # fill the level
            jb = bdict[tuple(bo)] # get index of the out vector
            ii[ib] = ib # store index 
            jj[ib] = jb # store index 
            s = 1 # sign of the matrix element
            if (k+l+i+j)%2==1: s*= -1 # change the sign
            vals[ib] = a*s # store value
          mi = csc_matrix((vals,(jj,ii)),shape=(on,on)) #  matrix
          mout = mout + mi # add contribution
  return mout # return matrix


def write_basis(basis):
  """Write the basis in basis.out"""
  fo = open("basis.out","w")
  fo.write("# SIZE = "+str(len(basis))+"\n")
  for b in basis: # loop over states
    for ib in b: fo.write(str(ib)+"   ")
    fo.write("\n")
  fo.close()




def one2many(m,ne=None):
  """Return a matrix in the many body basis, with ne electrons"""
  if 2**m.shape[0]>dimmax: raise # too big system
  if ne is None: # no constrain
    f = lambda x: True
  else: # with constrain
    f = constrain_nelectrons(ne) # function to constrain the number of electrons
  nsites = m.shape[0] # shape of the matrix
  basis = generate_basis(nsites,constrain=f) # generate basis
  return one2many_basis(m,basis)




if __name__=="__main__":
  nsites = 3 # 3 electrons
  basis = generate_basis(nsites) # return the basis
  m = np.random.random((nsites,nsites))
  m = m + m.transpose()
  bdict = get_dictionary(basis) 
  one2many(m,basis,bdict)
  print(basis)




