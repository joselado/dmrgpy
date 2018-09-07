from __future__ import print_function
import numpy as np

def get_moments_dmrg(self,n=1000):
  """Get the moments with DMRG"""
  self.setup_sweep("accurate")
  self.setup_task("dos",task={"nkpm":str(n)})
  self.write_hamiltonian() # write the Hamiltonian to a file
  self.run() # perform the calculation
  return np.genfromtxt("KPM_MOMENTS.OUT").transpose()[0]



def get_moments_spismj_dmrg(self,n=1000,i=0,j=0,smart=True):
  """Get the moments with DMRG"""

  self.setup_sweep()
  task= {"nkpm":str(n),"kpmmaxm":str(self.kpmmaxm),
                "site_i_kpm":str(i),"site_j_kpm":str(j),
                "kpm_scale":str(self.kpmscale)}
  if smart: task["smart_kpm_window"] = "true" 
  self.setup_task("spismj",task=task) 
  self.write_hamiltonian() # write the Hamiltonian to a file
  self.run() # perform the calculation
  return np.genfromtxt("KPM_MOMENTS.OUT").transpose()[0]



def get_moments_dynamical_correlator_dmrg(self,n=1000,i=0,j=0,name="XX"):
  """Get the moments with DMRG"""
  self.setup_sweep("accurate")
  fcorr = ["cdc","cdcup","ccd","cdcdn","deltadelta"] # correlators for fermions
  if len(name)==2:
    if name[0]=="X": namei="Sx"
    elif name[0]=="Y": namei="Sy"
    elif name[0]=="Z": namei="Sz"
    if name[1]=="X": namej="Sx"
    elif name[1]=="Y": namej="Sy"
    elif name[1]=="Z": namej="Sz"
    else: raise
    if self.sites[i] !=1 or self.sites[j]!=1:
        if name!="ZZ": raise  # fermions only accept ZZ
  elif name in fcorr: # fermionic correlator
    if self.sites[i] !=1 or self.sites[j]!=1: raise # only for fermions
    if name=="cdc": namei = "Cdag" ; namej = "Cdag"
    elif name=="cdcup": namei = "Cdagup" ; namej = "Cdagup"
    elif name=="cdcdn": namei = "Cdagdn" ; namej = "Cdagdn"
    elif name=="ccd": namei = "C" ; namej = "C"
    elif name=="deltadelta" or name=="delta":
        namei = "delta" ; namej = "delta"
    else: raise
  task= {"nkpm":str(n),"kpmmaxm":str(self.kpmmaxm),
                "site_i_kpm":str(i),"site_j_kpm":str(j),
                "kpm_scale":str(self.kpmscale),
                "kpm_cutoff":str(self.kpmcutoff),
                "kpm_operator_i":namei,"kpm_operator_j":namej}
  self.setup_task("dynamical_correlator",task=task) 
  self.write_hamiltonian() # write the Hamiltonian to a file
  self.run() # perform the calculation
  m = np.genfromtxt("KPM_MOMENTS.OUT").transpose()
#  return m[1]
  return m[0]+1j*m[1]





import pychain.kpm
from pychain.kpm import generate_profile

def get_dos(self,n=1000,mode="DMRG",ntries=10):
  if mode=="DMRG": 
#  if False: 
    mus = [get_moments_dmrg(self,n=n) for i in range(ntries)] # get the moments
    scale = np.genfromtxt("DOS_KPM_SCALE.OUT") # scale of the dos
  else:
    m  = self.get_full_hamiltonian()
    mus = [pychain.kpm.random_trace(m/15.0,ntries=1,n=1000)
                 for i in range(ntries)]
    scale = 1./15.
  mus = np.mean(np.array(mus),axis=0)
  mus = mus[0:100]
  xs = 0.99*np.linspace(-1.0,1.0,2000,endpoint=True) # energies
  ys = generate_profile(mus,xs,use_fortran=False).real # generate the DOS
  xs /= scale
  ys *= scale
  np.savetxt("DOS.OUT",np.matrix([xs,ys]).T)
  return (xs,ys)


def get_spismj(self,n=1000,mode="DMRG",i=0,j=0,smart=True,window=[-1,10]):
  self.to_folder() # go to temporal folder
  if mode=="DMRG": 
# get the moments
    mus = get_moments_spismj_dmrg(self,n=n,i=i,j=j,smart=smart) 
    if smart: # smart energy window
      scale,shift = np.genfromtxt("DOS_KPM_SCALE.OUT") # scale of the dos
    else: # convertional brute force window
      scale = np.genfromtxt("DOS_KPM_SCALE.OUT") # scale of the dos
      shift = 0.0 # energy shift
    # check that the moments do not take absur values
    mmu = np.max(np.abs(mus[1:n//10])) # average
    xs = 0.99*np.linspace(-1.0,1.0,n*10,endpoint=True) # energies
    ys = generate_profile(mus,xs,use_fortran=False,kernel="lorentz").real # generate the DOS
    xs -= shift # energy shift
    xs /= scale
    ys *= scale
#    e0 = self.gs_energy() # ground state energy
    e0 = np.genfromtxt("GS_ENERGY.OUT") # ground state energy
    xs -= e0 # substract GS energy
    # now retain only an energy window
  else: 
    h = self.get_full_hamiltonian()
    sc = self.get_pychain()
    import pychain.correlator
    if mode=="fullKPM":
      (xs,ys) = pychain.correlator.spismj_kpm(sc,h,n=n,i=i,j=j)
    elif mode=="full" or mode=="ED":
      delta = float(self.ns)/n*1.5
      (xs,ys) = pychain.correlator.spismj(sc,h,delta=delta,i=i,j=j)
    else: raise
  self.to_origin() # go to origin folder
  (xs,ys) = restrict_interval(xs,ys,window) # restrict the interval
  np.savetxt("SPISPJ.OUT",np.matrix([xs.real,ys.real]).T)
  return (xs.real,ys.real)


def restrict_interval(x,y,window):
  """Restrict the result to a certain energy window"""
  if window is None: return (x,y)
  i = np.argwhere(x<window[0]) # last one
  j = np.argwhere(x>window[1]) # last one
  if len(i)==0: i = 0
  else: i = i[0][-1]
  if len(j)==0: j = len(x)
  else: j = j[0][0]
  return x[i:j].real,y[i:j]






def get_dynamical_correlator(self,n=1000,mode="DMRG",i=0,j=0,
             window=[-1,10],name="XX",delta=None):
  if delta is not None: # estimate the number of polynomials
    scale = 0.
    for s in self.sites:
        if s>1: scale += s**2
        if s==1: scale += 4.
    n = int(scale/(2.*delta))
  self.to_folder() # go to temporal folder
  if mode=="DMRG": 
# get the moments
    mus = get_moments_dynamical_correlator_dmrg(self,n=n,i=i,j=j,name=name) 
    scale = np.genfromtxt("DOS_KPM_SCALE.OUT") # scale of the dos
    e0 = np.genfromtxt("GS_ENERGY.OUT") # ground state energy
    minx = -1.0 + (e0*scale) 
    maxx = -1.0 
    xs = 0.99*np.linspace(-1.0,1.0,n*10,endpoint=True) # energies
    ys = generate_profile(mus,xs,use_fortran=False,kernel="lorentz") # generate the DOS
    xs /= scale
    ys *= scale
#    e0 = self.gs_energy() # ground state energy
    xs -= e0 # substract GS energy
    # now retain only an energy window
  else: 
    h = self.get_full_hamiltonian()
    sc = self.get_pychain()
    import pychain.correlator
    if delta is None: delta = float(self.ns)/n*1.5
    if mode=="fullKPM":
      (xs,ys) = pychain.correlator.dynamical_correlator_kpm(sc,h,n=n,i=i,j=j,
                         namei=name[0],namej=name[1])
    elif mode=="ED":
      (xs,ys) = pychain.correlator.dynamical_correlator(sc,h,delta=delta,i=i,
                        j=j,namei=name[0],namej=name[1])
    else: raise
  self.to_origin() # go to origin folder
  (xs,ys) = restrict_interval(xs,ys,window) # restrict the interval
  from scipy.interpolate import interp1d
  fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
  fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
  ne = int(10*(window[1] - window[0])/delta) # number of energies
  xs = np.linspace(window[0],window[1],ne)
  ys = fr(xs) + 1j*fi(xs) # evaluate the interpolator
  np.savetxt("DYNAMICAL_CORRELATOR.OUT",np.matrix([xs.real,ys.real,ys.imag]).T)
  return (xs,ys)



