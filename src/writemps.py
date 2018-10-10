import numpy as np


def write_hamiltonian(self):
    write_sites(self) # write the different sites
    write_couplings(self)  # write the couplings
    write_hoppings(self)  # write the hoppings
    write_spinful_hoppings(self)  # write the hoppings
    write_pairing(self)  # write the pairing
    write_hubbard(self)  # write hubbard terms
    write_fields(self) # write the fields


def write_hubbard(self):
  """Write couplings in a file"""
  fo = open("hubbard.in","w")
  cs = self.hubbard
  fo.write(str(len(cs))+"\n")
  for key in self.hubbard: # loop
    c = self.hubbard[key] # loop
    if self.sites[c.i]!=1: raise
    if self.sites[c.j]!=1: raise
    fo.write(str(c.i)+"  ")
    fo.write(str(c.j)+"  ")
    fo.write(str(c.g.real)+"\n")
  fo.close()



def write_hoppings(self):
  """Write couplings in a file"""
  fo = open("hoppings.in","w")
  cs = self.hoppings
  fo.write(str(len(cs))+"\n")
  for key in self.hoppings: # loop
    c = self.hoppings[key] # loop
    if self.sites[c.i]!=1: raise
    if self.sites[c.j]!=1: raise
    fo.write(str(c.i)+"  ")
    fo.write(str(c.j)+"  ")
    fo.write(str(c.g.real)+"  ")
    fo.write(str(c.g.imag)+"\n")
  fo.close()





def write_spinful_hoppings(self):
  """Write couplings in a file"""
  fo = open("spinful_hoppings.in","w")
  cs = self.spinful_hoppings
  if type(cs)==type(dict()): # dictionary type
    fo.write(str(4*len(cs))+"\n")
    for key in self.spinful_hoppings: # loop
      c = self.hoppings[key] # loop
      if self.sites[c.i]!=1: raise
      if self.sites[c.j]!=1: raise
      # loop over elements
      g = np.matrix(c.g) # convert to matrix
      for i in range(2):
        for j in range(2):
          fo.write(str(2*c.i+i)+"  ")
          fo.write(str(2*c.j+j)+"  ")
          fo.write(str(g[i,j].real)+"  ")
          fo.write(str(g[i,j].imag)+"\n")
  else: # assume matrix type
      from scipy.sparse import coo_matrix
      cs = coo_matrix(cs) # transform to coo matrix
      row = cs.row
      col = cs.col
      data = cs.data
      fo.write(str(len(data))+"\n")
      for (r,c,d) in zip(row,col,data):
          fo.write(str(r)+"   ")
          fo.write(str(c)+"   ")
          fo.write(str(d.real)+"   ")
          fo.write(str(d.imag)+"\n")
  fo.close()







def write_pairing(self):
  """Write couplings in a file"""
  fo = open("pairing.in","w")
  cs = self.pairing
  fo.write(str(len(cs))+"\n")
  for key in self.pairing: # loop
    c = self.pairing[key] # loop
    if self.sites[c.i]!=1: raise
    if self.sites[c.j]!=1: raise
    fo.write(str(c.i)+"  ")
    fo.write(str(c.j)+"  ")
    fo.write(str(c.g.real)+"  ")
    fo.write(str(c.g.imag)+"\n")
  fo.close()









def write_fields(self):
  """Write fields in a file"""
  fo = open("fields.in","w")
  fo.write(str(len(self.fields))+"\n")
  for i in range(len(self.fields)): # loop
    fo.write(str(i)+"  ")
    fo.write(array2string(self.fields[i])+"\n")
  fo.close()



def write_sites(self):
  fo = open("sites.in","w")
  fo.write(str(self.ns)+"\n") # write number of sites
  for si in self.sites:
    if si<7: fo.write(str(si)+"\n")
    elif si==1: print("Warning, some sites are not spin operators")
    else: raise
  fo.close()




def write_couplings(self):
  """Write couplings in a file"""
  fo = open("couplings.in","w")
  cs = self.couplings
  out = [] # empty list
  stored = [] # empty list
  for c in cs:
    for ii in range(3):
      for jj in range(3):
        l = c.g[ii,jj]
        if np.abs(l)!=0.0: # if nonzero
            out += [(c.i,c.j,ii,jj,l.real,l.imag)] # store
            if (c.i,c.j,ii,jj) in stored: 
                print("Repeated coupling",c.i,c.j,ii,jj)
                raise
            stored += [(c.i,c.j,ii,jj)] # store
  fo.write(str(len(out))+"\n") # write that number
  for o in out: # loop over elements
        (i,j,ii,jj,p,q) = o
        fo.write(str(i)+"  ")
        fo.write(str(j)+"  ")
        fo.write(str(ii)+"  ")
        fo.write(str(jj)+"  ")
        fo.write(str(p)+"  ")
        fo.write(str(q)+"\n")
  fo.close()



def write_correlators(pairs):
  """Write the pairs of correlators in a file"""
  fo = open("correlators.in","w") # open the file
  fo.write(str(len(pairs))+"\n") # write in a file
  for p in pairs:
    # the first has to be smaller than the second
    fo.write(str(p[0])+"  "+str(p[1])+"\n")
  fo.close()




def array2string(m):
  out = ""
  for j in m:
    out += "  "+str(j)
  return out

def matrix2string(m):
  """Convert a 3x3 matrix into an string"""
  out = " "
  m = np.array(m)
  for i in m:
    for j in i:
      out += "  "+str(j)
  return out # return





def write_sweeps(self):
  """Write sweep info"""
  fo = open("sweeps.in","w")
  fo.write("sweeps\n{\n")
  fo.write("nsweeps = "+str(self.sweep["n"])+"\n")
  fo.write("maxm = "+str(self.sweep["maxm"])+"\n")
  fo.write("cutoff = "+str(self.sweep["cutoff"])+"\n}\n")
  fo.close()

