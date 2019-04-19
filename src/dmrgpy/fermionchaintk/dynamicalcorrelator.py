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
    def getd(ii,jj):
        return self.get_dynamical_correlator_spinless(
            name="densitydensity",i=2*i+ii,j=2*j+jj,**kwargs)
    def getcdc(ii,jj):
        return = self.get_dynamical_correlator_spinless(
                name="cdc",i=2*i+ii,j=2*j+jj,**kwargs)
    def getcc(ii,jj):
        return = self.get_dynamical_correlator_spinless(
                name="cc",i=2*i+ii,j=2*j+jj,**kwargs)
    if name=="densitydensity":
        (es,uu) = getd(0,0) # up up
        (es,ud) = getd(0,1) # up up
        (es,dd) = getd(1,1) # down down
        (es,du) = getd(1,0) # down down
        return (es,uu+dd+ud+du) # return the contributions
    elif name=="cdc":
        (es,uu) = getcdc(0,0)
        (es,dd) = getcdc(1,1)
        return (es,uu+dd) # return the contributions
    elif name=="cc":
        (es,uu) = getcc(0,0)
        (es,dd) = getcc(1,1)
    elif name=="cdcup": return getcdc(0,0)
    elif name=="cdcdn": return getcdc(1,1)
    elif name=="ccup": return getcc(0,0)
    elif name=="ccdn": return getcc(1,1)
    elif name=="ZZ":
        (es,uu) = getd(0,0) # up up
        (es,dd) = getd(1,1) # down down
        (es,ud) = getd(0,1) # up down
        (es,du) = getd(1,0) # down up
        return (es,(uu+dd-ud-du)/4.0) # return the contributions
    elif name=="densityZ":
        (es,uu) = getd(0,0) # up up
        (es,dd) = getd(1,1) # down down
        (es,ud) = getd(0,1) # up down
        (es,du) = getd(1,0) # down up
        return (es,(uu-dd-ud+du)/2.0) # return the contributions
    else: raise # not implemented




