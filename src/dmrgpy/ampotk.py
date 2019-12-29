# routines to use autompo features
import numpy as np
from . import multioperator

def names(i):
    if i==0: return "Sx"
    if i==1: return "Sy"
    if i==2: return "Sz"

def exchange2ampo(self):
  """Return ampo list for exchange"""
  try: self.Sx
  except: return 0
  h = 0 # initialize
  S = [self.Sx,self.Sy,self.Sz]
  for o in self.exchange: # loop over couplings
      for i in range(3):
        for j in range(3):
            c = o.g[i,j]
            h = h + c*S[i][o.i]*S[j][o.j]
  return h # return list

def field2ampo(self):
  """Return ampo list for exchange"""
  try: self.Sx
  except: return 0
  h = 0 # initialize
  for ii in range(len(self.fields)): # loop over couplings
      o = self.fields[ii] # this site
      h = h + o[0]*self.Sx[ii]
      h = h + o[1]*self.Sy[ii]
      h = h + o[2]*self.Sz[ii]
  return h # return list


def hoppings2ampo(self):
    h = 0
    for key in self.hoppings:
        i,j = key[0],key[1]
        c = self.hoppings[key].g
        h = h + c*self.Cdag[i]*self.C[j]
    return h

def hubbard2ampo(self):
    h = 0
    for key in self.hubbard:
        i,j = key[0],key[1]
        c = self.hubbard[key].g
        h = h + c*self.N[i]*self.N[j]
    return h

def pairing2ampo(self):
    if self.pairing is None: return 0
    ns = self.ns
    h = 0
    for i in range(ns):
      for j in range(ns):
          c = self.pairing(i,j)
          out = c*self.C[i]*self.C[j]
          h = h + out + out.get_dagger()
    return h


def vijkl2ampo(self):
    if self.vijkl is None: return 0
    h = 0
    n = self.ns
    C = self.C
    Cdag = self.Cdag
    for i in range(n):
      for j in range(n):
        for k in range(n):
          for l in range(n):
              h = h + self.vijkl(i,j,k,l)*Cdag[i]*C[j]*Cdag[k]*C[l]
    return h



def write_all(self):
    """Write everything in a file"""
    h = 0
    h = h + exchange2ampo(self) # get exchange
    h = h + field2ampo(self) # get magnetic field
    h = h + hoppings2ampo(self) # get hoppings
    h = h + hubbard2ampo(self) # get hubbard
    h = h + pairing2ampo(self) # get pairing
    h = h + vijkl2ampo(self) # get Vijkl
    self.execute(lambda: h.write(name="hamiltonian.in"))


