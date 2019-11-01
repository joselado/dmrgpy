import numpy as np

def one2many(ids,op):
   """Function to transform to many body basis given identity operators"""
   tmp = np.zeros((1,1),dtype=np.complex) # initialize
   for j in range(len(ids)): # loop over sites
       if i!=j: op2 = ids[j] # identity
       else: op2 = op # operator
       op = np.tensordot(op,ids[j]) # tensor product
   return op # return operator


