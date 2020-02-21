import numpy as np
from .kpmdmrg import restrict_interval
from .algebra import algebra
from .edtk import finitetemperature
from . import operatornames
from . import multioperator

def old2ampo(self):
    """Transform an old Hamiltonian into the AMPO form"""
    if self.use_ampo_hamiltonian: return self.hamiltonian
    else:
        h = 0
        Si = [self.Sx,self.Sy,self.Sz]
        for c in self.exchange:
            for i in range(3):
              for j in range(3): h = h + c.g[i,j]*Si[i][c.i]*Si[j][c.j]
        if len(self.fields)==len(self.Sz):
            for i in range(len(fields)):
                b = fields[i]
                for j in range(3): h = h + b[j]*Si[j][i]
        return h



# wrapper function for pychain

def get_full_hamiltonian(self):
    sc = get_pychain(self) # get pychain object
    if self.use_ampo_hamiltonian:
        return multioperator.MO2matrix(self.hamiltonian,sc)
    else: # conventional way
        h = old2ampo(self)
        return multioperator.MO2matrix(h,sc)
#      def get_coupling(i,j):
#        """Return the coupling between two sites"""
#        for c in self.exchange:
#          if i==c.i and j==c.j: return c.g
#        return np.zeros((3,3))
#      h = sc.add_tensor_interaction(get_coupling) # add interaction
#      h = h + sc.add_exchange(self.fields) # add magnetic fields
#      return h


def get_pychain(self):
  if self.pychain_object is None:
    from .pychain import build
    sc = build.Spin_chain()
    # the pychain library assumes that s=1/2 is for spin one-half
    # whereas in DMRG s = 2 is for S=1/2
    sc.build((np.array(self.sites)-1.)/2.) 
    sc.hamiltonian = self.hamiltonian
    return sc
  else: return self.pychain_object # return stored object


def get_dynamical_correlator(self,T=0.0,i=0,j=0,name="XX",**kwargs):
    """Return the dynamical correlator"""
    if T==0.0: # zero temperature
        return get_dynamical_correlator_T0(self,name=name,i=i,j=j,**kwargs)
    else: # finite temperature
        h = get_full_hamiltonian(self) # get the Hamiltonian
        sc = get_pychain(self) # get the object
        n1,n2 = operatornames.recognize(name) # return the names
        a = sc.get_operator(n1,i) # get operator
        b = sc.get_operator(n2,j) # get operator
        return finitetemperature.dynamical_correlator(h,a,b,T=T,**kwargs)




def get_dynamical_correlator_T0(self,submode="ED",
             window=[-1,10],name="XX",delta=2e-2,es=None,**kwargs):
  """
  Compute a dynamical correlator using the KPM-DMRG method
  """
  if delta is not None: # estimate the number of polynomials
    scale = 0.
    for s in self.sites:
        if s>1: scale += s**2 # spins
        else: scale += 4. # anything else
    npol = int(scale/(4.*delta))
  self.to_folder() # go to temporal folder
  h = self.get_full_hamiltonian()
  sc = self.get_pychain()
  from .pychain import correlator as pychain_correlator
  if delta is None: delta = float(self.ns)/n*1.5
  if submode=="KPM":
    (xs,ys) = pychain_correlator.dynamical_correlator_kpm(sc,h,n=npol,
                       namei=name[0],namej=name[1],es=es,**kwargs)
  elif submode=="ED" or submode=="CVM":
    if submode=="ED": mode = "full"
    elif submode=="CVM": mode = "cv"
    else: raise
    (xs,ys) = pychain_correlator.dynamical_correlator(sc,h,delta=delta,
                      namei=name[0],namej=name[1],mode=mode,es=es,
                      **kwargs)
  else: raise
  self.to_origin() # go to origin folder
  if es is None:
    (xs,ys) = restrict_interval(xs,ys,window) # restrict the interval
  else:
    (xs,ys) = restrict_interval(xs,ys,[min(es),max(es)]) 
  from scipy.interpolate import interp1d
  fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
  fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
  if es is None:
      ne = int(100*(window[1] - window[0])/delta) # number of energies
      xs = np.linspace(window[0],window[1],ne)
  else: xs = np.array(es).copy() # copy input array
  ys = fr(xs) + 1j*fi(xs) # evaluate the interpolator
  np.savetxt("DYNAMICAL_CORRELATOR.OUT",np.matrix([xs.real,ys.real,ys.imag]).T)
  return (xs,ys)


def gs_energy(self,T=0.0):
    """Return the ground state energy"""
    pyc = self.get_pychain()
    pyc.hamiltonian = self.get_hamiltonian()
    return pyc.gs_energy()
    h = get_full_hamiltonian(self)
    if T==0.0: # zero temperature
        return algebra.ground_state(h)[0] # return energy
    else: # non zero temperature
        return finitetemperature.gs_energy(h,beta=1./T) # return energy



def get_magnetization(sc,T=0.0):
    """Compute a static correlator"""
    scp = sc.get_pychain() # get pychain spinchain object
    h = get_full_hamiltonian(sc) # get Hamiltonian
    ns = sc.ns
    opx = [scp.sxi[i] for i in range(ns)]
    opy = [scp.syi[i] for i in range(ns)]
    opz = [scp.szi[i] for i in range(ns)]
    if T==0.0: # zero temperature
        wf = algebra.ground_state(h)[1] # get GS wavefunction
        mx = np.array([algebra.braket_wAw(wf,op).real for op in opx])
        my = np.array([algebra.braket_wAw(wf,op).real for op in opy])
        mz = np.array([algebra.braket_wAw(wf,op).real for op in opz])
    else: # non zero temperature
        mx = finitetemperature.measure(h,opx,beta=1./T).real
        my = finitetemperature.measure(h,opy,beta=1./T).real
        mz = finitetemperature.measure(h,opz,beta=1./T).real
    return (mx,my,mz)
