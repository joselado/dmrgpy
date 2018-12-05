from __future__ import print_function
import build
import read
import pyximport; pyximport.install()
import matplotlib.pyplot as plt
import numpy as np
import spectrum
import os
import examples
import entanglement



n = 5  # number of spins
#spins = [.5,.5,1,1.5,2.5,1.5,1,.5,.5]
spins = [.5 for i in range(n)]
#spin:w
#s = [1.5,.5]
#couplings = []
#for i in range(n):
#  for j in range(n):
#    couplings.append([i,j])
#couplings = [(0,1),(1,2),(2,3),(0,3)]
#js = [1.0,0.1,1.0,0.1] 
couplings = [(i,i+1) for i in range(len(spins)-1)] 
#couplings += [(n-1,0)]

js = [-1.0 for i in couplings]
# generate the input files

sc = build.Spin_chain()

sc.build(spins,couplings=couplings)
xs = [c[0] for c in couplings]
ys = [c[1] for c in couplings]
h0 = sc.generate_hamiltonian(xs,ys,js,mode="operator")

# kitaev trimer
ts = examples.kitaev_coupling() # kitaev coupling
#h0 = sc.add_tensor_interaction([0,1,2],[1,2,0],ts)

# uniaxial anisotropy
ds = [.1 for s in spins] # add anisotropy 

# correlation function

if True:
  ds = [-.1 for s in spins] # add anisotropy 
#  h0 = h0 + sc.add_uniaxial_anisotropy(ds)
#  ts = examples.kitaev_coupling() # kitaev coupling
#  h0 = sc.add_tensor_interaction([0,1,2],[1,2,0],ts)
  eig,evec = spectrum.eigenstates(h0,k=10,evals=True)
  v0 = evec[0] # ground state
  e0,waves = spectrum.ground_states(h0) # wave functions of the ground state
  print("GS energy and degeneracy",e0,len(waves))
  print("Total spin of the GS",[build.exp_val(w,sc.s) for w in waves])
  spectrum.write_splitting(waves,sc.szi) # write matrix elements
  # entanglement
  for wave in waves: # loop over waves
    dm = entanglement.reduced_density_matrix(wave,sc.basis)
    entropy = entanglement.entropy(wave,sc.basis)
#    print("Wavefunction",wave)
#    print("Density matrix",dm)
    print("Entropy",entropy,"\n")
  mins = entanglement.minimize_entropy(waves,sc.basis)
  print("Minimum entropy",mins)
  exit()
  os.system("rm *.OUT")
  spectrum.write_correlation(sc,v0) # write correlation matrix
  exit()


import scipy.sparse.linalg as lg
import scipy.linalg as dlg

n = 5 # number os eigenvalues to plot
zs = np.linspace(0.,.2,40)

for z in zs: # loop over fields 
  zee = [[z,z,z] for s in spins] # create zeeman fields
#  zee = [[0.,0.,0.] for s in spins] # create zeeman fields
#  zee[0] = [0.,0.,z] # create zeeman fields
  h = h0 + sc.add_exchange(zee) # add exchange field
  ds = [-z for s in spins] # add anisotropy 
#  h = h0 + sc.add_uniaxial_anisotropy(ds)
  h = build.traceless(h) # traceless hamiltonian
  eig = spectrum.eigenstates(h,k=100)
  print("Done",z)
  plt.scatter([z for e in eig],eig,marker="o")




plt.show()
