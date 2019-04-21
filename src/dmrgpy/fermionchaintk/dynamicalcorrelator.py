import numpy as np
from .. import multioperator
from . import mop


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
        return self.get_dynamical_correlator_spinless(
                name="cdc",i=2*i+ii,j=2*j+jj,**kwargs)
    def getccd(ii,jj):
        return self.get_dynamical_correlator_spinless(
                name="ccd",i=2*i+ii,j=2*j+jj,**kwargs)
    def getcc(ii,jj):
        return self.get_dynamical_correlator_spinless(
                name="cc",i=2*i+ii,j=2*j+jj,**kwargs)
    ### Worksround for four field operators
    def caca(ii,jj,kk,ll): # four field operators
        mi = multioperator.obj2MO([["Cdag",2*i+jj],["C",2*i+ii]]
            ,name="kpm_multioperator_i")
        mj = multioperator.obj2MO([["Cdag",2*j+kk],["C",2*j+ll]]
            ,name="kpm_multioperator_j")
        return self.get_dynamical_correlator_spinless(name=(mi,mj),
                **kwargs)
    def aaaa(ii,jj,kk,ll): # four field operators
        mi = multioperator.obj2MO([["Cdag",2*i+jj],["Cdag",2*i+ii]]
            ,name="kpm_multioperator_i")
        mj = multioperator.obj2MO([["C",2*j+kk],["C",2*j+ll]]
            ,name="kpm_multioperator_j")
        return self.get_dynamical_correlator_spinless(name=(mi,mj),
                **kwargs)
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
    elif name=="ccd":
        (es,uu) = getccd(0,0)
        (es,dd) = getccd(1,1)
        return (es,uu+dd) # return the contributions
    elif name=="cc":
        (es,uu) = getcc(0,0)
        (es,dd) = getcc(1,1)
        return (es,uu+dd) # return the contributions
    elif name=="cdcup": return getcdc(0,0)
    elif name=="cdcdn": return getcdc(1,1)
    elif name=="ccdup": return getccd(0,0)
    elif name=="ccddn": return getccd(1,1)
    elif name=="ccup": return getcc(0,0)
    elif name=="ccdn": return getcc(1,1)
    elif name=="ZZ":
        mi = mop.get_sz(i=i,name="kpm_multioperator_i") # multioprator for i
        mj = mop.get_sz(i=j,name="kpm_multioperator_j") # multioprator for i
        return self.get_dynamical_correlator_spinless(name=(mi,mj),
                **kwargs)
#        (es,uu) = getd(0,0) # up up
#        (es,dd) = getd(1,1) # down down
#        (es,ud) = getd(0,1) # up down
#        (es,du) = getd(1,0) # down up
#        return (es,(uu+dd-ud-du)/4.0) # return the contributions
    elif name=="XX":
        mi = mop.get_sx(i=i,name="kpm_multioperator_i") # multioprator for i
        mj = mop.get_sx(i=j,name="kpm_multioperator_j") # multioprator for i
        return self.get_dynamical_correlator_spinless(name=(mi,mj),
                **kwargs)
#        (es,out) = caca(0,1,0,1)
#        out = out + caca(0,1,1,0)[1] 
#        out = out + caca(1,0,1,0)[1] + caca(1,0,1,0)[1] 
#        return (es,out/4.0) # return the contributions
    elif name=="YY":
#        mi = mop.get_sy(i=i,name="kpm_multioperator_i") # multioprator for i
#        mj = mop.get_sy(i=j,name="kpm_multioperator_j") # multioprator for i
#        return self.get_dynamical_correlator_spinless(name=(mi,mj),
#                **kwargs)
        (es,out) = caca(0,1,0,1)
        out = out - caca(0,1,1,0)[1] 
        out = out + caca(1,0,1,0)[1] - caca(1,0,1,0)[1] 
        return (es,-out/4.0) # return the contributions
    elif name=="deltadelta": # swave pairing
        return aaaa(0,1,0,1) # return the swave pairing amplitude
    elif name=="densityZ":
        (es,uu) = getd(0,0) # up up
        (es,dd) = getd(1,1) # down down
        (es,ud) = getd(0,1) # up down
        (es,du) = getd(1,0) # down up
        return (es,(uu-dd-ud+du)/2.0) # return the contributions
    else: raise # not implemented




