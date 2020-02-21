from .manybodychain import Many_Body_Hamiltonian
import numpy as np
import scipy.linalg as lg
from .pyfermion import mbfermion
from .algebra import algebra
from .fermionchaintk import dynamicalcorrelator
from .fermionchaintk import staticcorrelator
from .fermionchaintk import hamiltonian
from . import funtk

class Fermionic_Hamiltonian(Many_Body_Hamiltonian):
    """Class for fermionic Hamiltonians"""
    def __init__(self,n):
        self.C = [self.get_operator("C",i) for i in range(n)]
        self.N = [self.get_operator("N",i) for i in range(n)]
        self.Cdag = [self.get_operator("Cdag",i) for i in range(n)]
        Many_Body_Hamiltonian.__init__(self,[0 for i in range(n)])
        self.fermionic = True
        self.use_ampo_hamiltonian = True # use ampo
    def set_hoppings(self,fun):
        """Add the spin independent hoppings"""
        self.set_hoppings_MB(fun)
    def get_density_spinless(self,**kwargs):
        """Return the electronic density"""
        return staticcorrelator.get_density_spinless(self,**kwargs)
    def get_density(self,**kwargs):
        """Return the electronic density"""
        return staticcorrelator.get_density_spinless(self,**kwargs)
    def set_hubbard_spinless(self,fun):
        """ Hubbard term """
        hamiltonian.set_hubbard_spinless(self,fun)
    def set_hubbard(self,fun):
        """ Hubbard term """
        hamiltonian.set_hubbard_spinless(self,fun)
    def get_density_fluctuation_spinless(self,**kwargs):
        """Return the electronic density fluctuations"""
        return staticcorrelator.get_density_fluctuation_spinless(self,**kwargs)
    def get_density_fluctuation(self,**kwargs):
        """Return the electronic density fluctuations"""
        return staticcorrelator.get_density_fluctuation_spinless(self,**kwargs)
    def get_pairing(self):
        """
        Return the superfluid density
        """
        return staticcorrelator.get_pairing_spinless(self,**kwargs)
    def vev_spinless(self,MO,mode="DMRG",**kwargs):
        """ Return a vaccum expectation value"""
        if mode=="DMRG":
            return self.vev_MB(MO,**kwargs)
        elif mode=="ED": 
            MBF = self.get_ED_obj() # get the object
            return MBF.vev(MO,**kwargs) # return overlap
        else: raise # unrecognized
    def vev(self,MO,**kwargs): return self.vev_spinless(MO,**kwargs)
    def excited_vev_spinless(self,MO,mode="DMRG",**kwargs):
        """ Return a vaccum expectation value"""
        if mode=="DMRG": return self.excited_vev_MB(MO,**kwargs)
        elif mode=="ED": return self.get_ED_obj().excited_vev(MO,**kwargs) 
    def excited_vev(self,MO,**kwargs): 
        return self.excited_vev_spinless(MO,**kwargs)
    def hamiltonian_free(self,pairs=[[]]):
        """
        Return the free part of the fermionic Hamiltonian
        """
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
            MBF = self.get_ED_obj() # get the object
            return algebra.lowest_eigenvalues(MBF.h,**kwargs)
    def get_correlator_spinless(self,**kwargs):
          """
          Compute static correlator
          """
          return staticcorrelator.get_correlator_spinless(self,**kwargs)
    def get_correlator(self,**kwargs):
          """
          Compute static correlator
          """
          return staticcorrelator.get_correlator_spinless(self,**kwargs)
    def get_dynamical_correlator(self,**kwargs):
        """
        Compute a dynamical correlator, standard name
        """
        return dynamicalcorrelator.get_dynamical_correlator_spinless(self,
                **kwargs)
    def get_dynamical_correlator_spinless(self,**kwargs):
        """
        Compute a dynamical correlator, specific function for spinless
        """
        return dynamicalcorrelator.get_dynamical_correlator_spinless(self,
                **kwargs)
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
        if mode=="DMRG": 
            return Many_Body_Hamiltonian.gs_energy(self,**kwargs)
        elif mode=="ED": 
#            if np.max(np.abs(self.hubbard_matrix))<1e-6 and self.vijkl is None:
#                return self.gs_energy_free()
#            else:
                MBF = self.get_ED_obj()
                return algebra.lowest_eigenvalues(MBF.h,n=1)[0]
        else: raise # unrecognised
    def get_ED_obj(self):
        """
        Return the many body fermion object
        """
        MBf = mbfermion.MBFermion(self.ns) # create object
        MBf.add_multioperator(self.hamiltonian) 
        return MBf # return the object
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
    def __init__(self,n):
        """ Rewrite the init method to take twice as many sites"""
        super().__init__(2*n) # initialize the Hamiltonian
        self.Sx = [0.5*self.Cdag[2*i]*self.C[2*i+1] +
                0.5*self.Cdag[2*i+1]*self.C[2*i] for i in range(n)]
        self.Sy = [-0.5*1j*self.Cdag[2*i]*self.C[2*i+1] +
                1j*0.5*self.Cdag[2*i+1]*self.C[2*i] for i in range(n)]
        self.Sz = [0.5*self.Cdag[2*i]*self.C[2*i] +
                (-1)*0.5*self.Cdag[2*i+1]*self.C[2*i+1] for i in range(n)]
        self.Delta = [0.5*self.C[2*i]*self.C[2*i+1] for i in range(n)]
        self.Cup = [self.C[2*i] for i in range(n)]
        self.Cdagup = [self.Cdag[2*i] for i in range(n)]
        self.Cdn = [self.C[2*i+1] for i in range(n)]
        self.Cdagdn = [self.Cdag[2*i+1] for i in range(n)]
        self.Nup = [self.get_operator("N",2*i) for i in range(n)]
        self.Ndn = [self.get_operator("N",2*i+1) for i in range(n)]
        self.Ntot = [self.Nup[i]+self.Ndn[i] for i in range(n)]
    def get_density_spinful(self,**kwargs):
        """
        Return the density in each site, summing over spin channels
        """
        return staticcorrelator.get_density_spinful(self,**kwargs)
    def get_density(self,**kwargs):
        """
        Return the density in each site, summing over spin channels
        """
        return staticcorrelator.get_density_spinful(self,**kwargs)
    def get_magnetization(self,**kwargs):
        """Return magnetization"""
        return staticcorrelator.get_magnetization_spinful(self,**kwargs)
    def get_onsite_pairing(self,**kwargs):
        """
        Return the expectation value of the onsite pairing
        """
        return staticcorrelator.get_onsite_pairing_spinful(self,**kwargs)
    def get_dynamical_correlator_spinful(self,**kwargs):
        """Return the dynamical correlator of an spinful system"""
        return dynamicalcorrelator.get_dynamical_correlator_spinful(self,
                **kwargs)
    def get_dynamical_correlator(self,**kwargs):
        """Return the dynamical correlator of an spinful system"""
        return self.get_dynamical_correlator_spinful(**kwargs)
    def set_hubbard_spinful(self,fun):
        """
        Add Hubbard interation in a spinful manner
        The Hubbard term will be defined as
        n_i n_j, with n_i = n_{i,up} + n_{i,,down}
        """
        hamiltonian.set_hubbard_spinful(self,fun)
    def set_hubbard(self,fun):
        """
        Add Hubbard interation in a spinful manner
        The Hubbard term will be defined as
        n_i n_j, with n_i = n_{i,up} + n_{i,,down}
        """
        hamiltonian.set_hubbard_spinful(self,fun)
    def set_swave_pairing(self,fun):
        """
        Add onsite swave pairing to a spinful Hamiltonian
        The pairing term is of the form
        Delta_i c_{i,up} c_{i,down} + h.c.
        """
        hamiltonian.set_swave_pairing_spinful(self,fun)
    def get_density_fluctuation_spinful(self,**kwargs):
        """Return the electronic density"""
        return staticcorrelator.get_density_fluctuation_spinful(self,**kwargs)
    def get_density_fluctuation(self,**kwargs):
        """Return the electronic density"""
        return staticcorrelator.get_density_fluctuation_spinful(self,**kwargs)
    def get_correlator_spinful(self,**kwargs):
        """
        Get a static correlator
        """
        return staticcorrelator.get_correlator_spinful(self,**kwargs)
    def get_correlator(self,**kwargs):
        """
        Get a static correlator
        """
        return staticcorrelator.get_correlator_spinful(self,**kwargs)
    def set_exchange(self,fun):
        """
        Add exchange coupling betwwen the spinful fermionic sites
        """
        hamiltonian.set_exchange_spinful(self,fun) # set the exchange
    def set_hoppings_spinful(self,fun):
        """
        Function to Add hopping in a spinful manner
        """
        def fun2(i,j):
            if i%2==j%2: return fun(i//2,j//2)
            return 0.0
        self.set_hoppings(fun2)








