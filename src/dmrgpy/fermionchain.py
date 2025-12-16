from .manybodychain import Many_Body_Chain
import numpy as np
import scipy.linalg as lg
from .pyfermion import mbfermion
from .algebra import algebra
from .fermionchaintk import staticcorrelator
from .fermionchaintk import hamiltonian
from . import funtk
from . import gap

class Fermionic_Chain(Many_Body_Chain):
    """Class for fermionic Hamiltonians"""
    def __init__(self,n,**kwargs):
        self.C = [self.get_operator("C",i) for i in range(n)]
        self.Cdag = [self.get_operator("Cdag",i) for i in range(n)]
        self.A = [self.get_operator("A",i) for i in range(n)]
        self.Adag = [self.get_operator("Adag",i) for i in range(n)]
        self.N = [self.get_operator("N",i) for i in range(n)]
#        self.N = [self.Cdag[i]*self.C[i] for i in range(n)]
        self.Id = self.get_operator("Id",1)
        Many_Body_Chain.__init__(self,[0 for i in range(n)],**kwargs)
        self.fermionic = True
        self.use_ampo_hamiltonian = True # use ampo
    def get_charge_gap(self,**kwargs):
        """Return the charge gap"""
        return gap.sector_gap(self,sum(self.N),**kwargs)
    def set_hoppings(self,fun):
        """Add the spin independent hoppings"""
        self.set_hoppings_MB(fun)
    def get_logdimension(self):
        return len(self.C)*np.log(2) # log dimension
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
    def vev_spinless(self,MO,**kwargs):
        """ Return a vacuum expectation value"""
        return self.vev(MO,**kwargs)
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
        m = np.zeros((self.ns,self.ns),dtype=np.complex128) # matrix
        for key in self.hoppings:
              t = self.hoppings[key]
              m[t.i,t.j] = t.g
        return m
#    def get_excited(self,mode="DMRG",**kwargs):
#          """
#          Wrapper for static correlator
#          """
#          if mode=="DMRG": # using DMRG
#            return Many_Body_Chain.get_excited(self,**kwargs)
#          elif mode=="ED":
#            MBF = self.get_ED_obj() # get the object
#            return algebra.lowest_eigenvalues(MBF.h,**kwargs)
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
            return Many_Body_Chain.gs_energy(self,**kwargs)
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
        if self.has_ED_obj: # if the ED object has been computed
            return self.ED_obj # return the stored object
        else:
            MBf = mbfermion.MBFermion(self.ns) # create object
            MBf.add_multioperator(self.hamiltonian) # add the Hamiltonian
            self.ED_obj = MBf # store the object
            self.has_ED_obj = True # set to True
            return self.ED_obj # return the object
    def execute(self,f):
        """
        This is a temporal fix to use the C operators in Julia ITensor
        """
        if self.itensor_version=="julia": # use the fermionic representation
            from . import multioperator
            multioperator.use_jordan_wigner = False
        return Many_Body_Chain.execute(self,f)






def get_gr_free(self,es=np.linspace(-10.,10.,800),delta=0.1,i=0,j=0):
    m = self.hamiltonian_free() # get the single body matrix
#    print(m)
    y = np.zeros(es.shape[0],dtype=np.complex128) # output
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



class Majorana_Chain(Fermionic_Chain):
    def __init__(self,n):
        """ Rewrite the init method to take twice as many sites"""
        nf = (n+1)//2 # number of conventional fermions
        nf = max([2,nf]) # fix
        super().__init__(nf) # initialize the Hamiltonian
        # define the Majorana operators
        G = [0 for i in range(nf*2)] # empty list
        for jf in range(nf): # loop over fermions
            G[2*jf] = (self.C[jf] + self.Cdag[jf])/np.sqrt(2)
            G[2*jf+1] = 1j*(self.C[jf] - self.Cdag[jf])/np.sqrt(2)
        self.G = [G[i] for i in range(n)] # sotre those operators
        del self.C  # clean
        del self.Cdag  # clean
        del self.N # clean





class Spinful_Fermionic_Chain(Fermionic_Chain):
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
        self.Sz = [0.5*self.N[2*i] +
                (-1)*0.5*self.N[2*i+1] for i in range(n)]
        self.Delta = [0.5*self.C[2*i]*self.C[2*i+1] for i in range(n)]
        self.Cup = [self.C[2*i] for i in range(n)]
        self.Cdagup = [self.Cdag[2*i] for i in range(n)]
        self.Cdn = [self.C[2*i+1] for i in range(n)]
        self.Cdagdn = [self.Cdag[2*i+1] for i in range(n)]
        self.Nup = [self.N[2*i] for i in range(n)]
        self.Ndn = [self.N[2*i+1] for i in range(n)]
        self.Ntot = [self.Nup[i]+self.Ndn[i] for i in range(n)]
        self.use_ampo_hamiltonian = True # use ampo
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



class Spinon_Chain(Spinful_Fermionic_Chain):
    """Class for spinon chains"""
    def get_gs(self,**kwargs):
        """Redefine the ground state method"""
        if self.computed_gs: return self.wf0 # return the wavefunction
        P = 1. # parton projector
        for i in range(len(self.Sx)): 
            P = P*(-2*self.Nup[i]*self.Ndn[i] + self.Ndn[i] + self.Nup[i])
            #P = P*(1.- self.Nup[i]*self.Ndn[i]) #  Gutzwiller projection
        from .mpsalgebra import mpsarnoldi
        super().gs_energy(**kwargs) # get the GS
        wf = self.wf0
#        wf = None # no initial guess
#        print("Projection",wf.dot(P*wf).real)
        wf = mpsarnoldi(self,self.hamiltonian,mode="GS",P=P,wf=wf)
        print("Projection",wf.dot(P*wf).real)
        self.computed_gs = True # computed GS
        self.wf0 = wf # store ground state
        return wf
    def gs_energy(self,**kwargs):
        wf = self.get_gs(**kwargs) # ground state
        return wf.dot(self.hamiltonian*wf).real # return energy










Fermionic_Hamiltonian = Fermionic_Chain
Spinful_Fermionic_Hamiltonian = Spinful_Fermionic_Chain


def isfermion(self):
    """Function to determine if an object is a valid fermionic object"""
    from .pyfermion.mbfermion import MBFermion
    if type(self)==Fermionic_Chain: return True
    if type(self)==Spinful_Fermionic_Chain: return True
    if type(self)==Spinful_F_Fermionic_Chain: return True
    if type(self)==MBFermion: return True
    else: return False
    



class Spinful_F_Fermionic_Chain(Fermionic_Chain):
    """
    Class to deal with fermionic Hamiltonians with
    spin degree of freedom
    """
    def __init__(self,n):
        """ Rewrite the init method to take twice as many sites"""
        super().__init__(3*n) # initialize the Hamiltonian
        self.Sx = [0.5*self.Cdag[3*i]*self.C[3*i+1] +
                0.5*self.Cdag[3*i+1]*self.C[3*i] for i in range(n)]
        self.Sy = [-0.5*1j*self.Cdag[3*i]*self.C[3*i+1] +
                1j*0.5*self.Cdag[3*i+1]*self.C[3*i] for i in range(n)]
        self.Sz = [0.5*self.N[3*i] +
                (-1)*0.5*self.N[3*i+1] for i in range(n)]
        self.Cup = [self.C[3*i] for i in range(n)]
        self.Cdagup = [self.Cdag[3*i] for i in range(n)]
        self.Cdn = [self.C[3*i+1] for i in range(n)]
        self.Cdagdn = [self.Cdag[3*i+1] for i in range(n)]
        self.F = [self.C[3*i+2] for i in range(n)]
        self.Fdag = [self.Cdag[3*i+2] for i in range(n)]
        self.Nup = [self.N[3*i] for i in range(n)]
        self.Ndn = [self.N[3*i+1] for i in range(n)]
        self.NF = [self.N[3*i+2] for i in range(n)]










