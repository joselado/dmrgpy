from .manybodychain import Many_Body_Hamiltonian
import numpy as np
import scipy.linalg as lg


class Fermionic_Hamiltonian(Many_Body_Hamiltonian):
    """Class for fermionic Hamiltonians"""
    def __init__(self,n,spinful=True):
        Many_Body_Hamiltonian.__init__(self,[1 for i in range(n)])
        self.spinful = spinful
    def get_density(self):
        """Return the electronic density"""
        m = self.get_file("MEASURE_N.OUT") # get the file
        return m.transpose()[1] # return density
    def get_density_fluctuation(self):
        """Return the electronic density"""
        d = self.get_file("MEASURE_N.OUT").transpose()[1] # get the file
        d2 = self.get_file("MEASURE_N2.OUT").transpose()[1] # get the file
        return d2-d**2 # return density fluctuations
    def get_delta(self):
        """Return the electronic density"""
        m = self.get_file("MEASURE_DELTA.OUT") # get the file
        return m.transpose()[1] # return delta
    def hamiltonian_free(self,pairs=[[]]):
        """Compute the free correlator"""
        if len(self.hubbard)!=0: raise
        else:
          m = np.zeros((self.ns*2,self.ns*2)) # matrix
          for key in self.hoppings:
              t = self.hoppings[key]
              m[2*t.i,2*t.j] = t.g
              m[2*t.i+1,2*t.j+1] = t.g
        if type(self.spinful_hoppings)!=type(dict()):
          m = m + self.spinful_hoppings
        return m
    def correlator_free(self,pairs=[[]]):
          """Get the correlator for free fermions"""
          m = self.hamiltonian_free() # get the single body matrix
          (es,vs) = lg.eigh(m) # diagonalize
          vs = vs.transpose()
          out = []
          for p in pairs:
              o = 0.0 # initialize
              for (e,v) in zip(es,vs):
                  if e<=0.0: 
                      if self.spinful: # spinful Hamiltonian
                          for i in range(2):
                            o += v[2*p[0]+i]*np.conjugate(v[2*p[1]+i]) # add
                      else: raise # not implemented
              out.append(o)
          return np.array(out) # return
    def gs_energy_free(self):
        """Get the energy for free fermions"""
        m = self.hamiltonian_free() # get the single body matrix
        es = lg.eigvalsh(m) # get the energies
        return np.sum(es[es<0.0]) # return energies
    def get_gr(self,**kwargs):
        return get_gr(self,**kwargs)
    def get_gr_free(self,**kwargs):
        return get_gr_free(self,**kwargs)





def get_gr_free(self,es=np.linspace(-10.,10.,800),delta=0.1,i=0,j=0):
    m = self.hamiltonian_free() # get the single body matrix
    print(m)
    y = np.zeros(es.shape[0],dtype=np.complex) # output
    iden = np.identity(m.shape[0])
    for ii in range(len(es)):
        yi = np.matrix(m-(es[ii]+1j*delta)*iden).I[i,j]
        y[ii] = yi # store
    return es,y













def get_gr(self,delta=0.002,es=np.linspace(-10.0,10.0,800),i=0,j=0):
    """Compute the advanced Green's function"""
    from . import kpmdmrg
    (x1,y1) = kpmdmrg.get_dynamical_correlator(self,es=es,i=i,j=j,
            name="cdc",delta=delta)
    (x2,y2) = kpmdmrg.get_dynamical_correlator(self,es=es,i=i,j=j,
            name="ccd",delta=delta)
    x1 = x1 + self.e0 # shift by the fermi energy
    x2 = x2 + self.e0 # shift by the fermi energy
    # define interpolating function
    from scipy.interpolate import interp1d
    f1r = interp1d(x1, y1.real,fill_value=0.0,bounds_error=False)
    f1i = interp1d(x1, y1.imag,fill_value=0.0,bounds_error=False)
    f2r = interp1d(x2, y2.real,fill_value=0.0,bounds_error=False)
    f2i = interp1d(x2, y2.imag,fill_value=0.0,bounds_error=False)
    # compute the result
    yr = f1r(es) #+ f2r(-es) # real part
    yi = f1i(es) #- f2i(-es) # imaginary part
    # now add the imaginary part
    from scipy.signal import hilbert
    y = yr + 1j*hilbert(yr) + 1j*yi - hilbert(yi)
#    y = 1j*y
#    y = 1j*yr
    return (es,y)









