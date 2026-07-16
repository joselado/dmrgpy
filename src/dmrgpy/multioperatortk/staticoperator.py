# library for immutable operators,
# they can act over a wavefunction, but they do not have
# algebra

from ..mps import MPS
from copy import copy as _shallow_copy

class StaticOperator():
    def __init__(self,MO,MBO):
        """Init, takes as input a multioperator and the MBO"""
        self.MBO = MBO # store the many-body object
        self.cpp_handle = MBO._session.build_operator(MO.to_terms())
    def __mul__(self,v):
        from ..multioperator import MultiOperator
        if type(v)==MPS: # input is an MPS
            handle = self.MBO._session.apply_pure_operator(self.cpp_handle,v.cpp_handle)
            return MPS(self.MBO,cpp_handle=handle).copy()
        elif type(v)==StaticOperator: # input is an MPO
            out = self.copy()
            out.cpp_handle = self.MBO._session.multiply_operators(self.cpp_handle,v.cpp_handle)
            return out
        elif type(v)==MultiOperator: # input is a multioperator
            return self*StaticOperator(v,self.MBO)
        else:
            print("Unrecognized type in __mul__",type(v))
            raise
    def __rmul__(self,v):
        from ..multioperator import MultiOperator
        if type(v)==MultiOperator: # input is a multioperator
            return StaticOperator(v,self.MBO)*self
        else: raise
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
