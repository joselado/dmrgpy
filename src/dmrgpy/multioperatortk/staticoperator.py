# library for immutable operators,
# they can act over a wavefunction, but they do not have
# algebra

from ..mps import MPS
from copy import copy as _shallow_copy
import numpy as np

class StaticOperator():
    def __init__(self,MO,MBO):
        """Init, takes as input a multioperator and the MBO"""
        self.MBO = MBO # store the many-body object
        if getattr(MBO,"use_cpp_extension",False) and MBO._session is not None:
            self.cpp_handle = MBO._session.build_operator(MO.to_terms())
            self.SO = None
        else:
            self.cpp_handle = None
            self.SO = generate_SO(MO,MBO) # generate the static operator
    def __mul__(self,v):
        from ..multioperator import MultiOperator
        if type(v)==MPS: # input is an MPS
            if self.cpp_handle is not None and v.cpp_handle is not None:
                handle = self.MBO._session.apply_pure_operator(self.cpp_handle,v.cpp_handle)
                return MPS(self.MBO,cpp_handle=handle).copy()
            return pure_applyoperator_dmrg(self.MBO,self.SO,v)
        elif type(v)==StaticOperator: # input is an MPO
            out = self.copy()
            if self.cpp_handle is not None and v.cpp_handle is not None:
                out.cpp_handle = self.MBO._session.multiply_operators(self.cpp_handle,v.cpp_handle)
            else:
                out.SO = mult_pureoperator(self.MBO,self.SO,v.SO)
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
        if self.cpp_handle is not None:
            out.cpp_handle = self.MBO._session.hermitian_operator(self.cpp_handle)
        else:
            out.SO = hermitian_mpo_operator(self.MBO,self.SO)
        return out
    def copy(self):
        # A shallow copy is enough, for the same reasons as MPS.copy()
        # (see mps.py): self.SO is immutable bytes (or None, in
        # cpp_handle mode), the opaque cpp_handle is never mutated in
        # place (every Chain operator method takes A by const reference
        # and returns a new handle), and self.MBO should stay shared. A
        # real deepcopy would also choke on cpp_handle, which has no
        # pickle/deepcopy support.
        return _shallow_copy(self)
    def __deepcopy__(self,memo):
        return self.copy()
    def trace(self):
        if self.cpp_handle is not None:
            return self.MBO._session.trace_operator(self.cpp_handle)
        return trace_pureoperator(self.MBO,self.SO)
    def aMb(self,wf1,wf2):
        if (self.cpp_handle is not None and wf1.cpp_handle is not None
                and wf2.cpp_handle is not None):
            return self.MBO._session.overlap_aMb_operator(
                    wf1.cpp_handle,self.cpp_handle,wf2.cpp_handle)
        return overlap_aMb_static(self.MBO,wf1,self.SO,wf2)



def pure_applyoperator_dmrg(self,A,wf):
    """Apply a pure operator to a many body wavefunction"""
    self.execute(lambda: wf.write(name="pureapplyoperator_wf1.mps"))
    task = {"pureapplyoperator":"true",
            "pureapplyoperator_wf0":"pureapplyoperator_wf1.mps",
            "pureapplyoperator_operator":"pureapplyoperator_operator.mpo",
            "pureapplyoperator_wf1":"pureapplyoperator_wf1.mps",
            }
    self.execute(lambda: open(self.path+"/pureapplyoperator_operator.mpo","wb").write(A))
    self.task = task
    self.execute(lambda : self.run()) # run calculation
    return MPS(self,name="pureapplyoperator_wf1.mps").copy() # copy


def generate_SO(A,MBO):
    """Generate a static many-body object"""
    task = {"gen_pureoperator":"true",
            "gen_pureoperator_operator_in":"gen_pureoperator_operator.in",
            "gen_pureoperator_operator_out":"gen_pureoperator_operator.mpo",
            }
    MBO.execute(lambda: A.write(name="gen_pureoperator_operator.in"))
    MBO.task = task
    MBO.execute(lambda : MBO.run()) # run calculation
    return open(MBO.path+"/gen_pureoperator_operator.mpo","rb").read()




def mult_pureoperator(self,A,B):
    """Apply a pure operator to a many body wavefunction"""
    task = {"multmpo_operator":"true",
            "multmpo_pureoperator_A":"multmpo_pureoperator_A.mpo",
            "multmpo_pureoperator_B":"multmpo_pureoperator_B.mpo",
            "multmpo_pureoperator_C":"multmpo_pureoperator_C.mpo",
            }
    self.execute(lambda: open(self.path+"/multmpo_pureoperator_A.mpo","wb").write(A))
    self.execute(lambda: open(self.path+"/multmpo_pureoperator_B.mpo","wb").write(B))
    self.task = task
    self.execute(lambda : self.run()) # run calculation
    return open(self.path+"/multmpo_pureoperator_C.mpo","rb").read()




def trace_pureoperator(self,A):
    """Compute the trace of a pure operator"""
    task = {"trace_mpo_operator":"true",
            "trace_pureoperator":"trace_pureoperator.mpo",
            }
    self.execute(lambda: open(self.path+"/trace_pureoperator.mpo","wb").write(A))
    self.task = task
    self.execute(lambda : self.run()) # run calculation
    m = self.execute( lambda : np.genfromtxt("TRACE.OUT")) # run calculation
    return m[0] + 1j*m[1]


def overlap_aMb_static(self,wf1,A,wf2):
    """Compute the overlap between wavefunctions, with A a StaticOperator"""
    task = {"overlap_aMb_static":"true",
            }
    self.task = task # override tasks
    self.execute(lambda: wf1.write(name="overlap_aMb_wf1.mps"))
    self.execute(lambda: wf2.write(name="overlap_aMb_wf2.mps"))
    self.execute(lambda: open(self.path+"/overlap_aMb_pureoperator.mpo","wb").write(A))
    self.execute(lambda: self.run()) # run calculation
    m = self.execute(lambda : np.genfromtxt("OVERLAP_aMb.OUT")) # read
    return m[0] + 1j*m[1]




def hermitian_mpo_operator(self,A):
    """Get the dagger of an operator"""
    task = {"hermitian_mpo_operator":"true",
            "mpo_pureoperator_A":"mpo_pureoperator_A.mpo",
            "mpo_pureoperator_B":"mpo_pureoperator_B.mpo",
            }
    self.execute(lambda: open(self.path+"/mpo_pureoperator_A.mpo","wb").write(A))
    self.task = task
    self.execute(lambda : self.run()) # run calculation
    return open(self.path+"/mpo_pureoperator_B.mpo","rb").read()
