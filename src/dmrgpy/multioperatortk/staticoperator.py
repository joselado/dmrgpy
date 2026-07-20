# library for immutable operators: they can act over a wavefunction,
# and (see __add__/__sub__/__neg__/__mul__/__truediv__ below) support
# direct-sum/scalar algebra with each other on itensor_version 2, 3
# and "python" (not "julia_live" yet)

from ..mps import MPS
from copy import copy as _shallow_copy

class StaticOperator():
    def __init__(self,MO,MBO):
        """Init, takes as input a multioperator and the MBO"""
        self.MBO = MBO # store the many-body object
        self.cpp_handle = MBO._session.build_operator(MO.to_terms())
    def __mul__(self,v):
        from ..multioperator import MultiOperator, isnumber
        if type(v)==MPS: # input is an MPS
            handle = self.MBO._session.apply_pure_operator(self.cpp_handle,v.cpp_handle)
            return MPS(self.MBO,cpp_handle=handle).copy()
        elif type(v)==StaticOperator: # input is an MPO
            out = self.copy()
            out.cpp_handle = self.MBO._session.multiply_operators(self.cpp_handle,v.cpp_handle)
            return out
        elif type(v)==MultiOperator: # input is a multioperator
            return self*StaticOperator(v,self.MBO)
        elif isnumber(v): # scalar rescale, no contraction/bond growth
            out = self.copy()
            out.cpp_handle = self.MBO._session.scale_operator(self.cpp_handle,complex(v))
            return out
        else:
            print("Unrecognized type in __mul__",type(v))
            raise
    def __rmul__(self,v):
        from ..multioperator import MultiOperator, isnumber
        if type(v)==MultiOperator: # input is a multioperator
            return StaticOperator(v,self.MBO)*self
        elif isnumber(v): return self*v
        else: raise
    def __truediv__(self,v): return self*(1./v)
    def __neg__(self): return (-1)*self
    def __radd__(self,v): return self + v
    def __add__(self,v):
        """Sum of two already-built MPOs, mirroring ITensorMPS.jl's
        `+(::MPO, ::MPO)` (abstractmps.jl): a compressed direct sum at the
        tensor-network level, without going back through the symbolic
        MultiOperator representation. MultiOperator already supports `+`
        on its own (see multioperator.py) for the common case of combining
        Hamiltonians before ever building an MPO; this is for combining
        operators that only exist as already-built StaticOperators (e.g.
        two independently constructed products/exponentials)."""
        if v==0: return self # allows sum([...]) starting from 0
        if type(v)!=StaticOperator:
            raise TypeError("Can only add a StaticOperator to another "
                    "StaticOperator (got "+str(type(v))+")")
        if not hasattr(self.MBO._session,"sum_operators"):
            raise NotImplementedError("MPO sum is not implemented for "
                    "this backend yet (itensor_version 2, 3 and "
                    "'python' support it, 'julia_live' doesn't)")
        out = self.copy()
        out.cpp_handle = self.MBO._session.sum_operators(self.cpp_handle,v.cpp_handle)
        return out
    def __sub__(self,v): return self + (-1)*v
    def get_dagger(self):
        out = self.copy()
        out.cpp_handle = self.MBO._session.hermitian_operator(self.cpp_handle)
        return out
    def copy(self):
        # A shallow copy is enough, for the same reasons as MPS.copy()
        # (see mps.py): the opaque cpp_handle is never mutated in place
        # (every Chain operator method takes A by const reference and
        # returns a new handle), and self.MBO should stay shared. A real
        # deepcopy would also choke on cpp_handle, which has no
        # pickle/deepcopy support.
        return _shallow_copy(self)
    def __deepcopy__(self,memo):
        return self.copy()
    def trace(self):
        return self.MBO._session.trace_operator(self.cpp_handle)
    def aMb(self,wf1,wf2):
        return self.MBO._session.overlap_aMb_operator(
                wf1.cpp_handle,self.cpp_handle,wf2.cpp_handle)
