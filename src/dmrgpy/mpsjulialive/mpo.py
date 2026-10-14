from ..multioperator import MO2list


def text_mpo(MO):
    """Transform a multioperator into the text input for Julia"""
    out = MO2list(MO)
    ls = []
    ls += [str(len(out))] # number of lines
    for o in out:
        n = (len(o)-2)//2 # number of terms
        if n>=200: raise ValueError("Too many factors in one term (%d), the limit is 199" % n) # C++ code needs to be recompiled
        ls += [str((len(o)-2)//2)] # number of terms
        ils = ""
        for io in o:
            ils += str(io)+"  "
        ls += [ils]
    return ls



def get_MPO(MO,MBO):
    """Transform a multioperator into a matrix product operator"""
    ls = text_mpo(MO)
    from .juliasession import Main, to_julia_strvec
    MPO = Main.toMPO(to_julia_strvec(ls),MBO.jlsites)
    return MPO

import numpy as np


class MPO():
    """Class for a Julia matrix product operator"""
    def __init__(self,MO,MBO=None):
        """Initialize, transforming a multioperator into a matrix
        product operator"""
        self.jlmpo = get_MPO(MO,MBO) # store
        self.MBO = MBO # store
        self.jlsites = MBO.jlsites # store the sites
    def __mul__(self,v):
        from ..multioperator import MultiOperator
        from .mps import MPS
        if type(v)==MPS: # input is an MPS
            from .juliasession import Main
            jlmps = Main.applyoperator(self.jlsites,self.jlmpo,
                    v.jlmps,self.MBO.maxm,self.MBO.cutoff)
            return MPS(jlmps,MBO=self.MBO)
        elif type(v)==MPO: # input is an MPO
            raise NotImplementedError("MPO*MPO is not implemented on julia_live") # not implemented
        elif type(v)==MultiOperator: # input is a multioperator
            self*MPO(v,MBO=self.MBO)
        else:
            raise TypeError("Unrecognized type in __mul__ %s" % (type(v),))
    def __rmul__(self,v):
        from ..multioperator import MultiOperator
        if type(v)==MultiOperator: # input is a multioperator
            return MPO(v,MBO=self.MBO)*self
        else: raise TypeError("Unrecognized type in __rmul__ %s" % (type(v),))
    def get_dagger(self):
        raise NotImplementedError("get_dagger is not implemented for julia_live MPOs")
    def copy(self):
        from copy import deepcopy
        return deepcopy(self)
    def trace(self):
        raise NotImplementedError("trace is not implemented for julia_live MPOs")
    def aMb(self,wf1,wf2):
        raise NotImplementedError("aMb is not implemented for julia_live MPOs")


