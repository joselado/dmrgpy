import numpy as np


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
  fo.write(str(len(cs))+"\n")
  for c in cs:
    fo.write(str(c.i)+"  ")
    fo.write(str(c.j)+"  ")
    fo.write(matrix2string(c.g)+"\n")
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
  fo.write("nsweeps = "+self.sweep["n"]+"\n")
  fo.write("maxm = "+self.sweep["maxm"]+"\n}\n")
  fo.write("cutoff = "+str(self.sweep["cutoff"])+"\n}\n")
  fo.close()

