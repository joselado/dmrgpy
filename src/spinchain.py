from __future__ import print_function
import numpy as np
import os
import kpmdmrg
import pychainwrapper

dmrgpath = os.environ["DMRGROOT"] # path for the program
one = np.identity(3)



class Coupling():
  def __init__(self,i,j,g):
    self.i = i
    self.j = j
    self.g = g


class Spin_Hamiltonian():
  def __init__(self,spins):
    self.spins = spins # list of the spins
    self.path = os.getcwd()+"/.mpsfolder" # folder of the calculations
    self.inipath = os.getcwd() # original folder
    self.ns = len(spins) # number of spins
    self.couplings = [Coupling(i,i+1,one) for i in range(self.ns-1)] # empty list
    self.fields = [] # empty list
#    self.couplings.append(Coupling(0,self.ns-1,one)) # closed boundary
    # additional arguments
    self.kpmmaxm = 10 # bond dimension in KPM
    self.kpmscale = 10.0
    os.system("mkdir -p "+self.path) # create folder for the calculations
  def to_folder(self): os.chdir(self.path) # go to calculation folder
  def to_origin(self): os.chdir(self.inipath) # go to original folder
  def clean(self): os.system("rm -r "+self.path) # clean temporal folder
  def set_exchange(self,fun):
    """Set the exchange coupling between spins"""
    self.couplings = [] # empty list
    for i in range(self.ns): # loop
      for j in range(self.ns):  # loop
        g = fun(i,j) # call the function
        g = g*one # multiply by the identity
        if np.sum(np.abs(g))!=0.0: 
          c = Coupling(i,j,g) # create class
          self.couplings.append(c) # store
  def set_fields(self,fun):
    self.fields = [fun(i) for i in range(self.ns)] # fields
  def get_coupling(self,i,j):
    """Return the coupling between two sites"""
    for c in self.couplings:
      if i==c.i and j==c.j: return c.g
    return np.zeros((3,3)) 
  def setup_sweep(self,mode="default"):
    setup_sweep(self,mode=mode)
    write_sweeps(self) # write the sweeps
  def setup_task(self,mode="GS",task=dict()):
    setup_task(self,mode=mode)
    for key in task: self.task[key] = task[key] # additional arguments
    write_tasks(self) # write the tasks
  def write_hamiltonian(self):
    write_sites(self) # write the different sites
    write_couplings(self)  # write the couplings
    write_fields(self) # write the fields
  def run(self): 
    os.system(dmrgpath+"/mpscpp/mpscpp > status.txt") # run the DMRG calculation
  def entropy(self,n=1):
    """Return the entanglement entropy"""
    self.setup_sweep()
    self.setup_task("entropy")
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.run() # perform the calculation
    return np.genfromtxt("ENTROPY.OUT")
  def get_full_hamiltonian(self):
    return pychainwrapper.get_full_hamiltonian(self)
  def get_pychain(self):
    return pychainwrapper.get_pychain(self)
  def get_dos(self,n=1000,mode="DMRG"):
    return kpmdmrg.get_dos(self,n=n,mode=mode)
  def get_spismj(self,n=1000,mode="DMRG",i=0,j=0,smart=False):
    return kpmdmrg.get_spismj(self,n=n,mode=mode,i=i,j=j,smart=smart)
  def get_excited(self,n=10,mode="DMRG"):
    self.to_folder() # go to temporal folder
    if mode=="DMRG":
      self.setup_sweep()
      self.setup_task("excited",task={"nexcited":str(n)})
      self.write_hamiltonian() # write the Hamiltonian to a file
      self.run() # perform the calculation
      out = np.genfromtxt("EXCITED.OUT")
    elif mode=="ED":
      h = self.get_full_hamiltonian() # get the Hamiltonian
      import pychain.spectrum
      out = pychain.spectrum.eigenstates(h,k=n) # return energy
    else: 
      self.to_origin() # go to main folder
      raise
    self.to_origin() # go to main folder
    return out
  def gs_energy(self,mode="DMRG"):
    """Return the ground state energy"""
    # write the spins
    self.to_folder() # go to temporal folder
    if mode=="DMRG":
      self.setup_sweep()
      self.setup_task("GS")
      self.write_hamiltonian() # write the Hamiltonian to a file
      self.run() # perform the calculation
      out = np.genfromtxt("GS_ENERGY.OUT") # return the ground state energy
    elif mode=="ED": # use brute force
      h = self.get_full_hamiltonian() # get the Hamiltonian 
      import pychain.spectrum
      out = pychain.spectrum.ground_state(h)[0] # return energy
    else: 
      self.to_origin() # go to main folder
      raise
    self.to_origin() # go to main folder
    return out
  def correlator(self,pairs=[[]],mode="DMRG"):
    if mode=="DMRG": # DMRG correlation
      self.setup_sweep()
      self.setup_task("correlator")
      self.write_hamiltonian() # write the Hamiltonian to a file
      write_correlators(pairs) # write the input file
      self.run() # perform the calculation
      return np.genfromtxt("CORRELATORS.OUT").transpose()[1] # return the correlators
    else: raise # not implemented
  def magnetization(self):
    """Calculate the magnetization of the system"""
    self.gs_energy() # calculate ground state
    mx = np.genfromtxt("MEASURE_SX.OUT").transpose()[1]
    my = np.genfromtxt("MEASURE_SY.OUT").transpose()[1]
    mz = np.genfromtxt("MEASURE_SZ.OUT").transpose()[1]
    return (mx,my,mz)




def array2string(m):
  out = ""
  for j in m:
    out += "  "+str(j)
  return out

def matrix2string(m):
  """Convert a 3x3 matrix into an string"""
  out = " "
  for i in m:
    for j in i:
      out += "  "+str(j)
  return out # return



def write_tasks(self):
  fo = open("tasks.in","w")
  fo.write("tasks\n{\n")
  for key in self.task:
    fo.write(key+" = "+self.task[key]+"\n")
  fo.write("}\n")
#("GS = true\ngap = false\ncorrelator = false\n}\n")
  fo.close()



def write_sweeps(self):
  """Write sweep info"""
  fo = open("sweeps.in","w")
  fo.write("sweeps\n{\nnsweeps = ")
  fo.write("nsweeps = "+self.sweep["n"]+"\n")
  fo.write("maxm = "+self.sweep["maxm"]+"\n}\n")
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
    if p[0]<p[1]:
      fo.write(str(p[0])+"  "+str(p[1])+"\n")
    elif p[0]>p[1]:
      fo.write(str(p[1])+"  "+str(p[0])+"\n")
    else: raise # error
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
  for si in self.spins: 
    if 1<si<7: fo.write(str(si)+"\n")
    else: raise
  fo.close()


def setup_sweep(self,mode="default"):
  """Setup the sweep parameters"""
  sweep = dict() # dictionary
  if mode=="default": # default mode
    sweep["n"] = "6"
    sweep["maxm"] = "40" 
  elif mode=="fast": # default mode
    sweep["n"] = "3"
    sweep["maxm"] = "40" 
  else: raise
  self.sweep = sweep # initialize




def setup_task(self,mode="GS"):
  """Setup the sweep parameters"""
  task = dict() # dictionary
  if mode=="GS": # default mode
    task["GS"] = "true"
  elif mode=="excited": # default mode
    task["excited"] = "true"
  elif mode=="correlator": # default mode
    task["GS"] = "true"
    task["correlator"] = "true" 
  elif mode=="entropy": # default mode
    task["GS"] = "true"
    task["entropy"] = "true" 
  elif mode=="dos": # default mode
    task["dos"] = "true" 
  elif mode=="spismj": # default mode
    task["spismj"] = "true" 
  else: raise
  self.task = task # initialize
