import numpy as np
from scipy.sparse import csc_matrix

def one2many(ids,op=None,i=-1):
   """Function to transform to many body basis given identity operators"""
   tmp = np.zeros((1,1),dtype=np.complex) # initialize
   tmp[0,0] = 1.0
   for j in range(len(ids)): # loop over sites
       if i!=j: op2 = ids[j] # identity
       else: op2 = op # operator
       tmp = np.kron(tmp,op2) # tensor product
   tmp = csc_matrix(tmp) # return operator
   tmp.eliminate_zeros()
   return tmp


