# routines to use autompo features


def names(i):
    if i==0: return "Sx"
    if i==1: return "Sy"
    if i==2: return "Sz"

def exchange2ampo(self):
  """Return ampo list for exchange"""
  out = [] # empty list
  for o in self.exchange: # loop over couplings
      for i in range(3):
        for j in range(3):
            c = o.g[i,j]
            if abs(c)>1e-5:
                out.append([c.real,c.imag,names(i),o.i+1,names(j),o.j+1])
  return out # return list

def field2ampo(self):
  """Return ampo list for exchange"""
  out = [] # empty list
  for ii in range(len(self.fields)): # loop over couplings
      o = self.fields[ii] # this site
      for i in range(3): # components
            c = o[i]
            if abs(o[i])>1e-5: # non-zero
                out.append([c.real,c.imag,names(i),ii+1])
  return out # return list






def write_all(self):
    """Write everything in a file"""
    out = [] # empty list
    out += exchange2ampo(self) # get exchange
    out += field2ampo(self) # get magnetic field
    write_ampo(out,"hamiltonian.in")


def write_ampo(out,name):
    f = open(name,"w")
    f.write(str(len(out))+"\n") # number of lines
    for o in out:
      f.write(str((len(o)-2)//2)+"\n") # number of terms
      for io in o:
          f.write(str(io)+"  ")
      f.write("\n")
    f.close()




