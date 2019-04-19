import numpy as np


def get_dynamical_correlator_spinless(self,name="densitydensity",
        mode="DMRG",**kwargs):
    """
    Compute a dynamical correlator for a spinless chain
    """
    if name is not "density" and self.spinful: raise
    if mode=="DMRG":
        return self.get_dynamical_correlator_MB(name=name,**kwargs)
    elif mode=="ED":
        MBF = self.get_MBF() # get the object
        return MBF.get_dynamical_correlator(name=name,**kwargs)
    else: raise




def get_dynamical_correlator_spinful(self,name="densitydensity",
        i=0,j=0,**kwargs):
    """Return the dynamical correlator of an spinful system"""
    if name=="densitydensity":
        (es,uu) = self.get_dynamical_correlator_spinless(
                name="densitydensity",i=2*i,j=2*j,**kwargs)
        (es,dd) = self.get_dynamical_correlator_spinless(
                name="densitydensity",i=2*i+1,j=2*j+1,**kwargs)
        return (es,uu+dd) # return the contributions
    elif name=="cdc":
        (es,uu) = self.get_dynamical_correlator_spinless(
                name="cdc",i=2*i,j=2*j,**kwargs)
        (es,dd) = self.get_dynamical_correlator_spinless(
                name="cdc",i=2*i+1,j=2*j+1,**kwargs)
        return (es,uu+dd) # return the contributions
    elif name=="ZZ":
        def getd(ii,jj):
            return self.get_dynamical_correlator_spinless(
                name="densitydensity",i=2*i+ii,j=2*j+jj,**kwargs)
        (es,uu) = getd(0,0) # up up
        (es,dd) = getd(1,1) # down down
        (es,ud) = getd(0,1) # up down
        (es,du) = getd(1,0) # down up
        return (es,(uu+dd-ud-du)/4.0) # return the contributions
    else: raise # not implemented




