import numpy as np


def get_correlator_spinless(self,name="cdc",mode="DMRG",**kwargs):
      """
      Wrapper for static correlator
      """
      if mode=="DMRG": # using DMRG
        return self.get_correlator_MB(name=name,**kwargs)
      elif mode=="ED": # using ED
        MBF = self.get_MBF() # get the object
        return MBF.get_correlator(name=name,**kwargs)
      else: raise




def get_correlator_spinful(self,name="ZZ",pairs=[[]],**kwargs):
    """
    Get a static correlator
    """
    # overwrite the method for spinful fermions
    pu = [[p1*2,p2*2] for (p1,p2) in pairs] # up pairs
    pd = [[p1*2+1,p2*2+1] for (p1,p2) in pairs] # down pairs
    pud = [[p1*2,p2*2+1] for (p1,p2) in pairs] # up-down pairs
    pdu = [[p1*2+1,p2*2] for (p1,p2) in pairs] # down-up pairs
    def f(pp):
      return self.get_correlator_spinless(
              name="densitydensity",pairs=pp,**kwargs)
    if name=="ZZ": # ZZ correlator
        return (f(pu) + f(pd) - f(pud) - f(pdu))/4.
    elif name=="densitydensity": # density-density correlator
        return f(pu) + f(pd) + f(pud) + f(pdu)
    else: raise











def get_density_spinless(self,**kwargs):
    """Return the electronic density"""
    pairs = [(i,i) for i in range(self.ns)]
    return self.get_correlator_spinless(pairs=pairs,
            name="cdc",**kwargs).real








def get_density_fluctuation_spinless(self,**kwargs):
    """Return the electronic density"""
    d = self.get_density_spinless(**kwargs) # get the density
    pairs = [(i,i) for i in range(self.ns)]
    d2 = self.get_correlator_spinless(pairs=pairs,name="densitydensity",
            **kwargs).real
    return d2-d**2 # return density fluctuations





def get_density_spinful(self,**kwargs):
    """
    Return the density in each site, summing over spin channels
    """
    ds = self.get_density_spinless(**kwargs) # get density
    return np.array([ds[2*i]+ds[2*i+1] for i in range(len(ds)//2)])






def get_magnetization_spinful(self,**kwargs):
    """Return magnetization"""
    pairsxc = [(2*i,2*i+1) for i in range(self.ns//2)]
    xc = self.get_correlator_spinless(pairs=pairsxc,
            name="cdc",**kwargs)
    mx = xc.real # get mx
    my = xc.imag # get my
    pairs = [(i,i) for i in range(self.ns)]
    ds = self.get_density_spinless(**kwargs) # get density
    mz = (np.array([ds[2*i]-ds[2*i+1] for i in range(len(ds)//2)]))/2.
    return np.array([mx,my,mz]).T # return magnetization








def get_onsite_pairing_spinful(self,**kwargs):
    """
    Return the expectation value of the onsite pairing
    """
    pairs = [(2*i,2*i+1) for i in range(self.ns//2)]
    cs = self.get_correlator_spinless(pairs=pairs,
            name="cc",**kwargs)
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
    d = self.get_density_spinful(**kwargs) # total density
    d2 = self.get_correlator_spinful(name="densitydensity", # fluctuation
            pairs=[(i,i) for i in range(self.ns//2)]).real
    return d2-d**2 # return density fluctuations

