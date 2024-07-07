from ..multioperator import MO2list


def text_mpo(MO):
    """Transform a multioperator into the text input for Julia"""
    out = MO2list(MO)
    ls = []
    ls += [str(len(out))] # number of lines
    for o in out:
        n = (len(o)-2)//2 # number of terms
        if n>=200: raise # C++ code needs to be recompiled
        ls += [str((len(o)-2)//2)] # number of terms
        ils = ""
        for io in o:
            ils += str(io)+"  "
        ls += [ils]
    return ls



def get_MPO(MO,MBO):
    """Transform a multioperator into a matrix product operator"""
    ls = text_mpo(MO)
#    print(ls) ; exit()
    from .juliasession import Main
    Main.mpo_ls = ls
#    print(type(Main.mpo_ls),type(MBO.jlsites))
#    Main.mpo_ls = Main.eval('convert(Vector{String}, mpo_ls)')
#    MPO = Main.eval('toMPO(mpo_ls)')
    MPO = Main.toMPO(Main.mpo_ls,MBO.jlsites)
#    print(MPO)
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
            raise # not implemented
        elif type(v)==MultiOperator: # input is a multioperator
            self*MPO(v,MBO=self.MBO)
        else:
            print("Unrecognized type in __mul__",type(v))
            raise
    def __rmul__(self,v):
        from ..multioperator import MultiOperator
        if type(v)==MultiOperator: # input is a multioperator
            return MPO(v,MBO=self.MBO)*self
        else: raise
    def get_dagger(self):
        raise
    def copy(self):
        from copy import deepcopy
        return deepcopy(self)
    def trace(self):
        raise
    def aMb(self,wf1,wf2):
        raise


