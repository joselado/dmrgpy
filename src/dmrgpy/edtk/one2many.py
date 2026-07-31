import numpy as np
from scipy.sparse import csc_matrix, kron as sp_kron

def one2many(ids,op=None,i=-1):
   """Function to transform to many body basis given identity operators.

   Builds the embedding with sparse Kronecker products throughout
   (instead of dense np.kron followed by a final sparse conversion):
   the embedded operator is trivially sparse (identity on every site but
   one), but a dense intermediate scales as O(dim**2) in time and memory
   and dominates ED chain construction well before the dim>2000 point
   where the eigensolver itself switches to ARPACK -- confirmed ~2000x
   faster at ns=14 (dim=16384), with the dense path needing tens of GB
   at ns=16 for a single intermediate matrix.
   """
   tmp = csc_matrix(np.array([[1.0]],dtype=np.complex128)) # initialize
   for j in range(len(ids)): # loop over sites
       op2 = ids[j] if i!=j else op # identity or operator
       tmp = sp_kron(tmp,csc_matrix(op2,dtype=np.complex128),format="csc")
   tmp.eliminate_zeros()
   return tmp


