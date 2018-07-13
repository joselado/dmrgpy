from __future__ import print_function
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import coo_matrix
import tensorialf90
import numpy as np

def tensorial_LO(op1,op2,sparse=True,fortran=True,adapted=True):
  """Perform the tensorial product, returning a LinearOperator"""
  op1 = coo_matrix(op1)
  op2 = coo_matrix(op2)
  dim1 = op1.shape[0] # size of the matrix
  dim2 = op2.shape[0] # size of the matrix
  nn1 = len(op1.col)
  nn2 = len(op2.col)
  col1 = op1.col+1
  row1 = op1.row+1
  data1 = op1.data
  col2 = op2.col+1
  row2 = op2.row+1
  data2 = op2.data
  if fortran:
    if np.sum(np.abs(op1-np.identity(op1.shape[0])))<0.000001 and adapted:
      def fun(v):
        """Function that returns a vector"""
        return tensorialf90.tensorial_idbv(dim1,
                  data2,row2,col2,v,dim2)
    elif np.sum(np.abs(op2-np.identity(op2.shape[0])))<0.000001 and adapted:
      def fun(v):
        """Function that returns a vector"""
        return tensorialf90.tensorial_aidv(dim2,data1,row1,col1,
                v,dim1)
    else:
      def fun(v):
        """Function that returns a vector"""
        return tensorialf90.tensorial_abv(data1,row1,col1,
                  data2,row2,col2,v,dim1,dim2)
  else:
    def fun(v):
      """Function that returns a vector"""
      out = v*0.0j # zero vector
      for ii in range(nn1): # loop over first subspace
        for jj in range(nn2): # loop over first subspace
          i1 = col1[ii]-1 
          j1 = row1[ii]-1
          i2 = col2[jj]-1 
          j2 = row2[jj]-1 
          i3 = dim2*i1 + i2
          j3 = dim2*j1 + j2
       #   print(dim1,dim2,i1,i2,j1,j2,i3,j3)
          out[i3] += data1[ii]*data2[jj]*v[j3]
      return out # return hte vector
  return LinearOperator( (dim1*dim2,dim1*dim2), matvec=fun )


def identity(n):
  def f(v):
    return v
  return LinearOperator( (n,n), matvec=f )


def slo2lo(op):
  """Create a linear operator from a sum of linear operators"""
  def f(v):
    return op*v
  n = op.shape[0]
  return LinearOperator( (n,n), matvec=f )



def exp_val(wf,A):
  """Calculate the expectation value of a wavefunction"""
  return np.conjugate(wf).dot(A*wf).real


