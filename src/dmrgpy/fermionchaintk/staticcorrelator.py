import numpy as np
from .. import correlator


def get_correlator_spinless(self,name="cdc",mode="DMRG",**kwargs):
      """
      Wrapper for static correlator

      The DMRG branch goes through correlator.get_correlator (which turns
      the name into a pair of per-site MultiOperators and takes their
      vev) rather than a self.get_correlator_MB method: that method has
      never existed -- it sits commented out in manybodychain.py -- so
      every DMRG-mode caller of this function raised AttributeError
      instead of returning a correlator. That took the whole spinful
      static-correlator path down with it (get_magnetization_spinful and
      get_pairing_spinless both route through here), which is how
      Spinful_Fermionic_Chain.get_magnetization() came to be dead on
      arrival.
      """
      if mode=="DMRG": # using DMRG
        return correlator.get_correlator(self,name=name,**kwargs)
      elif mode=="ED": # using ED
        MBF = self.get_ED_obj() # get the object
        return MBF.get_correlator(name=name,**kwargs)
      else: raise




def get_density_spinless(self,**kwargs):
    """Return the electronic density"""
    out = [self.vev(self.N[i],**kwargs) for i in range(self.ns)]
    return np.array(out).real








def get_density_fluctuation_spinless(self,**kwargs):
    """Return the electronic density"""
    d = self.get_density_spinless(**kwargs) # get the density
    d2 = np.array([self.vev(self.N[i]*self.N[i],**kwargs) for i in range(self.ns)])
    return d2-d**2 # return density fluctuations





def get_density_spinful(self,**kwargs):
    """
    Return the density in each site, summing over spin channels
    """
    ds = self.get_density_spinless(**kwargs) # get density
    return np.array([ds[2*i]+ds[2*i+1] for i in range(len(ds)//2)])






def get_magnetization_spinful(self,**kwargs):
    """Return the magnetization (mx,my,mz) of every orbital of a spinful
    fermionic chain, as an array of shape (norbitals,3).

    Orbital i occupies the two interleaved spinless sites 2i (up) and
    2i+1 (down), so with z = <Cdag_{2i} C_{2i+1}> = <cdag_up c_dn>,

        <Sx_i> = <(cdag_up c_dn + cdag_dn c_up)/2> = Re z
        <Sy_i> = <(-i cdag_up c_dn + i cdag_dn c_up)/2> = Im z
        <Sz_i> = (n_up - n_dn)/2

    (<cdag_dn c_up> is the conjugate of z, the two operators being each
    other's dagger.) This matches Spinful_Fermionic_Chain's own
    Sx/Sy/Sz MultiOperators exactly.
    """
    pairsxc = [(2*i,2*i+1) for i in range(self.ns//2)]
    xc = np.array(self.get_correlator_spinless(pairs=pairsxc,
            name="cdc",**kwargs))
    mx = xc.real # get mx
    my = xc.imag # get my
    ds = self.get_density_spinless(**kwargs) # get density
    mz = (np.array([ds[2*i]-ds[2*i+1] for i in range(len(ds)//2)]))/2.
    return np.array([mx,my,mz]).T # return magnetization








def get_onsite_pairing_spinful(self,**kwargs):
    """
    Return the expectation value of the onsite pairing
    """
    pairs = [(2*i,2*i+1) for i in range(self.ns//2)]
    cs = np.array([self.vev(self.C[i]*self.C[j],**kwargs) for (i,j) in pairs])
    return cs


def get_pairing_spinless(self,pairs=[[]],**kwargs):
    """
    Return the expectation value of the onsite pairing
    """
    cs = self.get_correlator_spinless(pairs=pairs,
            name="cc",**kwargs)
    return cs





def get_density_fluctuation_spinful(self,**kwargs):
    """Return the electronic density"""
    d = self.get_density(**kwargs) # total density
    d2 = [self.vev(Ni*Ni,**kwargs) for Ni in self.Ntot] # density fluctuation
    d2 = np.array(d2).real
    return d2-d**2 # return density fluctuations

