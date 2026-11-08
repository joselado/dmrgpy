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



def is_parity_odd(MO):
    """True if some term of MO has an odd number of fermionic (C-type)
    factors, the operators ITensorMPS's OpSum cannot turn into an MPO
    from the raw names"""
    from ..multioperatortk.canonical import _ODD
    for term in MO.op:
        if sum(1 for f in term[1:] if f[0] in _ODD) % 2 == 1: return True
    return False


def text_terms(terms):
    """Text input for Julia from [(coeff,[(name,1-based site),...]),...],
    the Jordan-Wigner-dressed form MultiOperator.to_terms() gives"""
    ls = [str(len(terms))] # number of lines
    for (c,factors) in terms:
        ls += [str(len(factors))] # number of factors
        ils = str(c.real)+"  "+str(c.imag)+"  "
        for (name,i) in factors: ils += str(name)+"  "+str(i)+"  "
        ls += [ils]
    return ls


def operator_text(MO):
    """Text input for Julia. A parity-even operator keeps its raw names,
    and ITensorMPS threads the Jordan-Wigner string itself. A parity-odd
    one (a single C or Cdag) is sent pre-dressed with A/Adag/F names,
    bosonic for OpSum (fermion_ops.jl), since OpSum refuses parity-odd
    fermionic terms."""
    if is_parity_odd(MO): return text_terms(MO.to_terms())
    return text_mpo(MO)


def get_MPO(MO,MBO):
    """Transform a multioperator into a matrix product operator"""
    ls = operator_text(MO)
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
    @classmethod
    def _from_jlmpo(cls,jlmpo,MBO):
        """An MPO object around an already-built Julia MPO"""
        out = cls.__new__(cls)
        out.jlmpo = jlmpo
        out.MBO = MBO
        out.jlsites = MBO.jlsites
        return out
    def __mul__(self,v):
        from ..multioperator import MultiOperator, isnumber
        from .mps import MPS
        from .juliasession import Main
        if type(v)==MPS: # input is an MPS
            jlmps = Main.applyoperator(self.jlsites,self.jlmpo,
                    v.jlmps,self.MBO.maxm,self.MBO.cutoff)
            return MPS(jlmps,MBO=self.MBO)
        elif type(v)==MPO: # operator product, v acts first
            # cut at the operator bond dimension mpomaxm, as v3 and
            # python do; only the sum uses the wavefunction maxm
            jlmpo = Main.mult_mpo(self.jlmpo,v.jlmpo,self.MBO.mpomaxm,
                    self.MBO.cutoff)
            return MPO._from_jlmpo(jlmpo,self.MBO)
        elif type(v)==MultiOperator: # input is a multioperator
            return self*MPO(v,MBO=self.MBO)
        elif isnumber(v): # scalar rescale, no contraction
            return MPO._from_jlmpo(Main.scale_mpo(self.jlmpo,complex(v)),
                    self.MBO)
        else:
            raise TypeError("Unrecognized type in __mul__ %s" % (type(v),))
    def __rmul__(self,v):
        from ..multioperator import MultiOperator, isnumber
        if type(v)==MultiOperator: # input is a multioperator
            return MPO(v,MBO=self.MBO)*self
        elif isnumber(v): return self*v
        else: raise TypeError("Unrecognized type in __rmul__ %s" % (type(v),))
    def __truediv__(self,v): return self*(1./v)
    def __neg__(self): return (-1)*self
    def __radd__(self,v): return self + v
    def __add__(self,v):
        """Sum of two already-built MPOs, compressed at the chain's maxm
        and cutoff"""
        if type(v)!=MPO:
            if type(v) in (int,float,complex) and v==0: return self # sum([...])
            raise TypeError("Can only add an MPO to another MPO (got "
                    +str(type(v))+")")
        from .juliasession import Main
        jlmpo = Main.sum_mpo(self.jlmpo,v.jlmpo,self.MBO.maxm,self.MBO.cutoff)
        return MPO._from_jlmpo(jlmpo,self.MBO)
    def __sub__(self,v): return self + (-1)*v
    def get_dagger(self):
        from .juliasession import Main
        return MPO._from_jlmpo(Main.dagger_mpo(self.jlmpo),self.MBO)
    def copy(self):
        # shallow, as StaticOperator.copy(): no method changes jlmpo in
        # place, and the chain must stay shared rather than deepcopied
        from copy import copy as shallow_copy
        return shallow_copy(self)
    def __deepcopy__(self,memo): return self.copy()
    def trace(self):
        """Full trace of the operator over the Hilbert space"""
        from .juliasession import Main
        return complex(Main.trace_mpo(self.jlmpo))
    def aMb(self,wf1,wf2):
        """<wf1|M|wf2>, exact, without truncating M|wf2>"""
        from .juliasession import Main
        return complex(Main.aMb_mpo(wf1.jlmps,self.jlmpo,wf2.jlmps))
