from __future__ import print_function
import numpy as np

def angular(x,y,z):
  xy = x*y - y*x
  xy = xy - 1j*z
  data = xy.data
  if len(data)>0:
    if np.max(np.abs(data))>0.00001:
      raise



def zero(x,y):
  xy = x*y - y*x
  if np.max(np.abs(xy.data))>0.00001:
    raise

