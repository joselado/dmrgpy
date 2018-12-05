import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

def get(name):
  return np.genfromtxt(name).transpose()

##################
# create figures #
##################
# energy
def get_energy():
  fig_energy = plt.figure()
  fig_energy.set_facecolor("white")
  ax_energy = fig_energy.add_subplot(1,1,1)
  ax_energy.set_title("Energy")
  
  m = get("ENERGY.OUT") 
  plt.scatter(m[0],m[1]) # scatter data
  return fig_energy

# entropy
def get_entropy():
  fig_entropy = plt.figure()
  plt.title("Entropy")
  fig_entropy.set_facecolor("white")
  m = get("ENTROPY.OUT")  # get entropies
  ns = len(m)-1
  for i in range(1,ns+1):  
  #  fig_entropy.add_subplot(ns,1,i)  # list of plots
    plt.plot(m[0],m[i],marker="o") # scatter data
  return fig_entropy


def get_map(name,cmap="seismic"):
  fig = plt.figure() # create figure
  fig.set_facecolor("white")
  ax = fig.add_subplot(111)
  m = get(name) # get data
  ns  = len(m) -1 # number of spins
  y = m[0] # x points
  x = list(range(1,ns+1)) # y points
  z = np.array([m[i] for i in range(1,ns+1)]).transpose() # z points
  minz = np.min(z)
  maxz = np.max(z)
  levels = np.linspace(minz,maxz, 80)
  plt.contourf(x,y,z,cmap=plt.get_cmap(cmap),levels=levels)
  plt.colorbar(label=name.split(".")[0])
  plt.ylabel("Time")
  plt.xlabel("Index")


def get_multiplot(name):
  fig = plt.figure() # create figure
  fig.set_facecolor("white")
  m = get(name) # get data
  ns  = len(m) -1 # number of spins
  x = m[0] # x points
  for i in range(1,len(m)):
    fig.add_subplot(len(m)-1,1,i)
    plt.plot(x,m[i],marker="o")
    plt.ylabel(name.split(".")[0])
  plt.xlabel("Time")


def get_sz():
  # sz expectation value
  fig_sz = plt.figure()
  fig_sz.set_facecolor("white")
  m = get("SZ.OUT")  # get entropies
  ns = len(m)-1 # number of spins
  for i in range(1,ns+1):  
    plt.plot(m[0],m[i],marker="o") # scatter data

def get_xy(name,ylabel="",xlabel="Time",s=10):
  # sz expectation value
  fig = plt.figure()
  fig.set_facecolor("white")
  m = get(name)  # get entropies
  for i in range(1,len(m)):
    plt.plot(m[0],m[i],c="red") # scatter data
    plt.scatter(m[0],m[i],c="red",s=s) # scatter data
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  return fig




if __name__=="__main__":
  #get_sz()
  get_map("SZ.OUT")
  get_map("ENTROPY.OUT")
  
  #get_entropy() # entropy plot
  
  plt.show()


