from .manybodychain import Many_Body_Hamiltonian
import numpy as np
import scipy.linalg as lg
from .pyfermion import mbfermion
from .algebra import algebra


class Fermionic_Hamiltonian(Many_Body_Hamiltonian):
    """Class for fermionic Hamiltonians"""
    def __init__(self,n,spinful=False):
        if spinful:
          Many_Body_Hamiltonian.__init__(self,[1 for i in range(n)])
        else:
          Many_Body_Hamiltonian.__init__(self,[0 for i in range(n)])
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
        """
        Return the superfluid density
        """
        m = self.get_file("MEASURE_DELTA.OUT") # get the file
        return m.transpose()[1] # return delta
    def hamiltonian_free(self,pairs=[[]]):
        """
        Return the free part of the fermionic Hamiltonian
        """
#        if len(self.hubbard)!=0: raise # not implemented
#        else: # everythin ok so far
        if self.spinful: # spinful
          m = np.zeros((self.ns*2,self.ns*2),dtype=np.complex) # matrix
          for key in self.hoppings:
              t = self.hoppings[key]
              m[2*t.i,2*t.j] = t.g
              m[2*t.i+1,2*t.j+1] = t.g
          if type(self.spinful_hoppings)!=type(dict()):
            m = m + self.spinful_hoppings
        else: # spinless Hamiltonian
          m = np.zeros((self.ns,self.ns),dtype=np.complex) # matrix
          for key in self.hoppings:
              t = self.hoppings[key]
              m[t.i,t.j] = t.g
        return m
    def get_excited(self,mode="DMRG",**kwargs):
          """
          Wrapper for static correlator
          """
          if mode=="DMRG": # using DMRG
            return Many_Body_Hamiltonian.get_excited(self,**kwargs)
          elif mode=="ED":
            MBF = self.get_MBF() # get the object
            return algebra.lowest_eigenvalues(MBF.h,**kwargs)
    def get_correlator(self,name="cdc",mode="DMRG",**kwargs):
          """
          Wrapper for static correlator
          """
          if mode=="DMRG": # using DMRG
            return Many_Body_Hamiltonian.get_correlator(self,name=name,**kwargs)
          elif mode=="ED": # using ED
            MBF = self.get_MBF() # get the object
            return MBF.get_correlator(name=name,**kwargs)
          else: raise
    def get_dynamical_correlator(self,name="densitydensity",
            mode="DMRG",**kwargs):
        """
        Compute a dynamical correlator
        """
        if name is not "density" and self.spinful: raise
        if mode=="DMRG":
#            if name is not "densitydensity": raise # this does not work
            return Many_Body_Hamiltonian.get_dynamical_correlator(self,
                    name=name,**kwargs)
        elif mode=="ED":
            MBF = self.get_MBF() # get the object
            return MBF.get_dynamical_correlator(name=name,**kwargs)
        else: raise
    def get_correlator_free(self,pairs=[[]]):
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
                      else: 
                            o += v[p[0]]*np.conjugate(v[p[1]]) # add
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
    def gs_energy(self,mode="DMRG",**kwargs):
        """Compute ground state energy, overrriding the method"""
        if mode=="DMRG": return Many_Body_Hamiltonian.gs_energy(self,**kwargs)
        elif mode=="ED": 
            if np.max(np.abs(self.hubbard_matrix))<1e-6 and self.vijkl is None:
                return self.gs_energy_free()
            else:
                MBF = self.get_MBF()
                return algebra.lowest_eigenvalues(MBF.h,n=1)[0]
        else: raise # unrecognised
    def get_MBF(self):
        """
        Return the many body fermion object
        """
        if not self.spinful: # spinless
          m0 = self.hamiltonian_free() # free Hamiltonian
          MBf = mbfermion.MBFermion(m0.shape[0]) # create object
          MBf.add_hopping(m0)
          MBf.add_hubbard(self.hubbard_matrix)
          MBf.add_vijkl(self.vijkl)
          return MBf # return the object
        else: raise 
    def get_kpm_scale(self):
        """Energy scale for KPM method"""
        return 4*self.ns*(2.+10*np.mean(np.abs(self.hubbard_matrix)))





def get_gr_free(self,es=np.linspace(-10.,10.,800),delta=0.1,i=0,j=0):
    m = self.hamiltonian_free() # get the single body matrix
#    print(m)
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




class Spinful_Fermionic_Hamiltonian(Fermionic_Hamiltonian):
    """
    Class to deal with fermionic Hamiltonians with
    spin degree of freedom
    """
    pass





