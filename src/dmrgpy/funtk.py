# toolkit to deal with functions and matrices

import numpy as np

def fun2list(f,n):
    """
    Transforma function into a list, where the first two entries
    are the indexes and the next one is the coupling
    """
    if f is None: return [] # empty list
    out = [] # empty list
    for i in range(n):
      for j in range(n):
        o = f(i,j)
        if np.abs(o)>1e-8:
            out.append([i,j,o]) # store
    return out # return list



def obj2fun(a):
    """
    Transforma a certain object in a callable function
    """
    if callable(a): return a
    elif type(a)==np.array: return lambda i,j: a[i,j]
    else: raise


def obj2mat(a):
    """
    Transform an object into a matrix
    """
    if type(a)==np.matrix: return a
    elif callable(a): raise

